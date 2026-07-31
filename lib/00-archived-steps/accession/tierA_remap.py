"""Tier A -- sequence-driven re-mapping with an author tiebreaker.

The tiers before this one *flag* a bad accession; none of them says what the row
should point at instead. Tier A closes that gap for the rows that are stuck:

  * obsolete_record whose Identical Protein Group is itself entirely suppressed, so
    Tier S recovered no live replacement (`recovered_accession` empty);
  * organism_mismatch / title_unmapped -- live accessions where the record's
    organism or annotation disagrees with PAZy, i.e. a candidate wrong-protein.

For each, it turns the row's sequence into a ranked set of *live* candidate
accessions and picks a best re-map target:

  1. candidates come cheap-first -- Identical Protein Group members (exact), then
     BLAST near-identical neighbours (blast.blastp, refseq_protein then nr);
  2. each candidate is scored on organism concordance (Tier O's comparator) and on
     cited-author overlap -- whether the candidate's own linked paper shares an
     author surname with the row's cited DOI, the tiebreaker that distinguishes the
     right protein from a look-alike;
  3. it also measures the *existing* accession's author overlap, which settles the
     organism_mismatch / title_unmapped rows the other way: if the flagged record's
     own literature shares an author with the cited paper, the flag is benign.

Tier A never re-labels a row (an obsolete accession must still be re-pointed by a
human). It attaches a ranked target + the evidence and leaves the row reviewable.
Absence of a candidate is silence, not proof -- same discipline as Tier S.
"""

from __future__ import annotations

import collections
import json
import os
import sys

import blast
import dataset
import ncbi
import normalize
import tier3_report
import tierO_organism
import tierP_provenance

# A candidate must be this close to be a plausible re-map: near-identical and
# near-full-length. Below this it is a homologue, not the same protein.
MIN_IDENTITY = 98.0
MIN_COVERAGE = 90.0
# Cap BLAST neighbours carried forward per row, best-scoring first, to bound the
# downstream docsum/elink/taxonomy volume. The author tiebreaker only needs the
# closest few; carrying more just multiplies elink->pubmed calls (mostly RefSeq
# entries with no linked paper) for no ranking gain.
MAX_BLAST_CANDS = 5

ORG_RANK = {"exact": 4, "same_species": 3, "same_genus": 2,
            "uncomparable": 1, "mismatch": 0}

FIELDS = [
    "row", "pazy_id", "enzyme_name", "accession", "source_label",
    "n_candidates", "remap_accession", "remap_source", "remap_identity",
    "remap_coverage", "remap_is_refseq", "remap_organism_flag",
    "remap_organism_detail", "remap_author_overlap", "remap_author_shared",
    "self_author_overlap", "tierA_flag", "candidates",
]

REFSEQ_PREFIXES = ("WP_", "NP_", "XP_", "YP_", "AP_")


def _key(acc: str) -> str:
    return normalize.strip_version(acc or "").upper()


def _is_refseq(acc: str, sourcedb: str) -> bool:
    return acc.upper().startswith(REFSEQ_PREFIXES) or "refseq" in (sourcedb or "").lower()


def _select_targets(t0s, t1s, t2s, tSs, tOs) -> list[dict]:
    """The rows Tier A works on, using tier3's own labeller as the source of truth.

    label_row is the single definition of what makes a row obsolete_record /
    organism_mismatch / title_unmapped, so reusing it keeps Tier A's scope from
    drifting from the report. Obsolete rows are taken only when Tier S found no IPG
    replacement -- the ones that still need an accession re-found.
    """
    targets = []
    for row, t0 in t0s.items():
        ts = tSs[row]
        label, _ = tier3_report.label_row(t0, t1s[row], t2s[row], ts, tOs[row])
        if label == "obsolete_record" and not ts["recovered_accession"]:
            targets.append({"t0": t0, "label": "obsolete_no_ipg"})
        elif label == "organism_mismatch":
            targets.append({"t0": t0, "label": "organism_mismatch"})
        elif label == "title_unmapped":
            targets.append({"t0": t0, "label": "title_unmapped"})
    return targets


def main() -> int:
    rows = {r.index: r for r in dataset.load_rows()}
    t0s = {r["row"]: r for r in dataset.read_json("tier0.json")}
    t1s = {r["row"]: r for r in dataset.read_json("tier1.json")}
    t2s = {r["row"]: r for r in dataset.read_json("tier2.json")}
    tSs = {r["row"]: r for r in dataset.read_json("tierS.json")}
    tOs = {r["row"]: r for r in dataset.read_json("tierO.json")}
    seqs = tierP_provenance.load_enriched_sequences()

    targets = _select_targets(t0s, t1s, t2s, tSs, tOs)
    limit = int(os.environ.get("TIERA_LIMIT", "0") or "0")
    if limit:
        targets = targets[:limit]
        print(f"tierA: TIERA_LIMIT={limit} -> processing {len(targets)} of the stuck rows")
    print(f"tierA: {len(targets)} stuck rows to re-map "
          f"({collections.Counter(t['label'] for t in targets)})")

    # --- pass 1a: IPG candidates (exact, cached) + collect each row's query -----
    row_cands: dict[int, dict[str, dict]] = {}
    own_key: dict[int, str] = {}
    queries: dict[int, str] = {}
    for t in targets:
        acc = t["t0"]["accession"]
        row = t["t0"]["row"]
        own = _key(acc)
        own_key[row] = own
        cands: dict[str, dict] = {}
        try:
            for m in ncbi.efetch_ipg(acc):
                k = _key(m["protein"])
                if k and k != own:
                    cands[k] = {"accession": m["protein"], "identity": 100.0,
                                "coverage": 100.0, "source": "ipg"}
        except Exception as exc:
            print(f"  ipg {acc}: {exc}", file=sys.stderr)
        row_cands[row] = cands
        seq = seqs.get(own, "")
        if seq:
            queries[row] = seq
        else:
            print(f"  {acc}: no query sequence available", file=sys.stderr)

    # --- pass 1b: BLAST all queries at once (refseq_protein, then nr fallback) ---
    print(f"tierA: BLASTing {len(queries)} sequences against refseq_protein")
    refseq_hits = blast.blast_many({str(r): s for r, s in queries.items()}, db="refseq_protein")
    row_blast: dict[int, list[dict]] = {}
    nr_needed: dict[int, str] = {}
    for row, seq in queries.items():
        keep = _passing(refseq_hits.get(str(row), []))
        if keep:
            row_blast[row] = keep
        else:
            nr_needed[row] = seq
    if nr_needed:
        print(f"tierA: {len(nr_needed)} rows had no refseq_protein hit; BLASTing nr")
        nr_hits = blast.blast_many({str(r): s for r, s in nr_needed.items()}, db="nr")
        for row in nr_needed:
            row_blast[row] = _passing(nr_hits.get(str(row), []))

    # merge BLAST neighbours into each row's candidate set (IPG members win on key)
    all_cand_accs: set[str] = set()
    for row, cands in row_cands.items():
        for h in row_blast.get(row, [])[:MAX_BLAST_CANDS]:
            k = _key(h["accession"])
            if not k or k == own_key[row] or k in cands:
                continue
            cands[k] = {"accession": h["accession"], "identity": h["pident"],
                        "coverage": h["coverage"], "source": "blast"}
        all_cand_accs |= {c["accession"] for c in cands.values()}

    # --- resolve every candidate once, drop non-current / non-protein ----------
    print(f"tierA: resolving {len(all_cand_accs)} distinct candidate accessions")
    cand_doc = _resolve_current(sorted(all_cand_accs))

    # --- taxonomy for organism concordance (candidates + PAZy) -----------------
    tax_ids = {cand_doc[k]["taxid"] for k in cand_doc if cand_doc[k]["taxid"]}
    for t in targets:
        tax_ids.add(tierO_organism._clean_taxid(rows[t["t0"]["row"]].ncbi_taxonomy_id))
    tax = ncbi.efetch_taxonomy(sorted(t for t in tax_ids if t))

    # --- author sets: candidates + each row's own accession + cited DOIs -------
    cand_uids = [cand_doc[k]["uid"] for k in cand_doc if cand_doc[k]["uid"]]
    own_uids = [t["t0"]["uid"] for t in targets if t["t0"].get("uid")]
    uid_pmids = _elink_pubmed(sorted(set(cand_uids) | set(own_uids)))

    cited_pmids_by_row: dict[int, list[str]] = {}
    all_dois = sorted({d for t in targets for d in normalize.split_dois(rows[t["t0"]["row"]].doi)})
    doi_pmid = ncbi.doi_to_pmid(all_dois)
    for t in targets:
        row = t["t0"]["row"]
        cited_pmids_by_row[row] = [doi_pmid[d] for d in normalize.split_dois(rows[row].doi)
                                   if d in doi_pmid]

    all_pmids = {p for ps in uid_pmids.values() for p in ps}
    all_pmids |= {p for ps in cited_pmids_by_row.values() for p in ps}
    authors_by_pmid = ncbi.pubmed_authors(sorted(all_pmids))

    def authors_of_uid(uid: str) -> set[str]:
        return {a for p in uid_pmids.get(str(uid), []) for a in authors_by_pmid.get(p, set())}

    def cited_authors(row: int) -> set[str]:
        return {a for p in cited_pmids_by_row.get(row, []) for a in authors_by_pmid.get(p, set())}

    # --- pass 2: score, rank, decide ------------------------------------------
    out = []
    for t in targets:
        t0 = t["t0"]
        row = t0["row"]
        pazy_tax = tierO_organism._clean_taxid(rows[row].ncbi_taxonomy_id)
        cited = cited_authors(row)

        scored = []
        for k, c in row_cands[row].items():
            doc = cand_doc.get(k)
            if not doc:  # candidate did not resolve to a live protein
                continue
            oflag, odetail = tierO_organism.compare(pazy_tax, doc["taxid"], tax)
            shared = authors_of_uid(doc["uid"]) & cited
            scored.append({
                "accession": doc["accession"], "source": c["source"],
                "identity": c["identity"], "coverage": c["coverage"],
                "is_refseq": _is_refseq(doc["accession"], doc["sourcedb"]),
                "organism_flag": oflag, "organism_detail": odetail,
                "author_overlap": len(shared), "author_shared": sorted(shared),
            })
        scored.sort(key=lambda s: (s["author_overlap"], ORG_RANK.get(s["organism_flag"], 0),
                                   s["identity"], s["is_refseq"]), reverse=True)

        self_overlap = len(authors_of_uid(t0.get("uid", "")) & cited)
        best = scored[0] if scored else None
        flag = _decide(best, self_overlap)

        rec = {
            "row": row, "pazy_id": t0["pazy_id"], "enzyme_name": t0["enzyme_name"],
            "accession": t0["accession"], "source_label": t["label"],
            "n_candidates": len(scored), "self_author_overlap": self_overlap,
            "tierA_flag": flag,
            "candidates": json.dumps(scored[:5]) if scored else "",
        }
        if best:
            rec.update({
                "remap_accession": best["accession"], "remap_source": best["source"],
                "remap_identity": best["identity"], "remap_coverage": best["coverage"],
                "remap_is_refseq": best["is_refseq"],
                "remap_organism_flag": best["organism_flag"],
                "remap_organism_detail": best["organism_detail"],
                "remap_author_overlap": best["author_overlap"],
                "remap_author_shared": "; ".join(best["author_shared"]),
            })
        else:
            rec.update({k: "" for k in (
                "remap_accession", "remap_source", "remap_identity", "remap_coverage",
                "remap_is_refseq", "remap_organism_flag", "remap_organism_detail",
                "remap_author_overlap", "remap_author_shared")})
        out.append(rec)

    out.sort(key=lambda r: r["row"])
    dataset.write_json("tierA.json", out)
    dataset.write_csv("tierA.csv", FIELDS, out)

    counts = collections.Counter(r["tierA_flag"] for r in out)
    print(f"\ntierA: {len(out)} rows")
    for k, v in counts.most_common():
        print(f"  {k:<18} {v}")
    strong = [r for r in out if r["tierA_flag"] == "remap_strong"]
    if strong:
        print("\n  strong re-maps (live, ≥99% id, organism-concordant, author-corroborated):")
        for r in strong[:12]:
            print(f"    {r['accession']:<18} -> {r['remap_accession']:<18} "
                  f"id={r['remap_identity']} {r['remap_organism_flag']}")
    return 0


def _passing(hits: list[dict]) -> list[dict]:
    """BLAST neighbours close enough to be the same protein, not just a homologue."""
    return [h for h in hits if h["pident"] >= MIN_IDENTITY and h["coverage"] >= MIN_COVERAGE]


def _resolve_current(accs: list[str]) -> dict[str, dict]:
    """Batch-resolve candidate accessions; keep only live protein records.

    Returns dict[versionless-upper] -> {accession, uid, taxid, sourcedb, title}.
    A suppressed/dead candidate is no better a re-map target than the row's own
    dead accession, so it is dropped here rather than ranked.
    """
    out: dict[str, dict] = {}
    # Small batches: docsum records are bulky and NCBI truncates big id lists.
    for batch in ncbi.chunked(sorted(set(accs)), 50):
        for cap, rec in ncbi.efetch_docsum(batch).items():
            status = (rec.get("status") or "").strip()
            replaced = (rec.get("replacedby") or "").strip()
            moltype = (rec.get("moltype") or "").lower()
            biomol = (rec.get("biomol") or "").lower()
            is_protein = moltype == "aa" or biomol == "peptide"
            is_current = status in ("", "live") and not replaced
            if not (is_protein and is_current):
                continue
            out[cap] = {
                "accession": rec.get("accessionversion") or cap,
                "uid": str(rec.get("gi") or rec.get("uid") or ""),
                "taxid": str(rec.get("taxid") or ""),
                "sourcedb": rec.get("sourcedb", ""),
                "title": rec.get("title", ""),
            }
    return out


def _elink_pubmed(uids: list[str]) -> dict[str, list[str]]:
    """protein uid -> linked PMIDs, bisecting bad id sets (never invents absence)."""
    out: dict[str, list[str]] = {}
    for batch in ncbi.chunked([u for u in uids if u], 50):
        out.update(ncbi.elink_pubmed_bisect(batch, db_from="protein"))
    return out


def _decide(best: dict | None, self_overlap: int) -> str:
    if best and best["identity"] >= 99.0 and ORG_RANK.get(best["organism_flag"], 0) >= 2 \
            and best["author_overlap"] >= 1:
        return "remap_strong"
    if best:
        return "remap_plausible"
    if self_overlap >= 1:
        return "self_corroborated"
    return "remap_none"


if __name__ == "__main__":
    raise SystemExit(main())
