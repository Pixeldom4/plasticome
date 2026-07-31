"""Tier S -- sequence identity to the cited paper.

This is the check that breaks Tier 2's ceiling. Tier 2 can only confirm a row when
NCBI's citation graph links the accession to the cited paper -- which never happens
for the mining/screening papers that merely *analysed* a sequence someone else
deposited. Those papers instead list the accessions and sequences they studied in
their supplementary files. Tier S:

  1. pulls every resolvable protein's amino-acid sequence from NCBI (the sheet's
     aa_sequence column is empty), filling it in;
  2. matches each row's accession and sequence against the cited papers' parsed
     supplementary material (supp.py) -- confirming rows Tier 2 structurally cannot;
  3. recovers a live replacement accession for suppressed records via the Identical
     Protein Group (their `replacedby` pointer is empty);
  4. flags rows that share an identical sequence (duplicate curation).

Absence of a supplementary match is silence, not error: most cited papers are not
open-access or deposit no parseable sequence table.
"""

from __future__ import annotations

import collections
import sys

import dataset
import extdb
import ncbi
import normalize
import supp

BATCH = 150
# Shortest window we will accept for a substring (not exact) sequence match, to
# absorb signal-peptide / expression-tag trimming without matching by chance.
MIN_SUBSTR = 60

FIELDS = [
    "row", "pazy_id", "enzyme_name", "accession", "aa_len", "sequence_source",
    "tierS_flag", "supp_doi", "match_kind", "recovered_accession", "recovered_source",
    "n_ipg_members", "duplicate_group", "duplicate_partners",
]


def _seq_match(row_seq: str, si_seqs: set[str]) -> bool:
    """Exact membership, then a length-guarded substring either direction."""
    if not row_seq:
        return False
    if row_seq in si_seqs:
        return True
    if len(row_seq) >= MIN_SUBSTR:
        for s in si_seqs:
            if len(s) >= MIN_SUBSTR and (row_seq in s or s in row_seq):
                return True
    return False


def main() -> int:
    rows = {r.index: r for r in dataset.load_rows()}
    tier0 = dataset.read_json("tier0.json")

    # --- 1. fetch aa_sequence for every row we can -------------------------
    # NCBI first (includes suppressed records -- is_protein, not is_current --
    # whose sequence we still want for IPG recovery and duplicate detection),
    # then fill the gaps from the databases the accession actually lives in.
    seq_by_key: dict[str, str] = {}
    src_by_key: dict[str, str] = {}

    def key_of(accession: str) -> str:
        return normalize.strip_version(accession).upper()

    def record(key: str, seq: str, src: str) -> None:
        if seq and key not in seq_by_key:
            seq_by_key[key], src_by_key[key] = seq, src

    protein_accs = sorted({
        t["accession"] for t in tier0 if t["resolves"] and t["is_protein"]
    })
    print(f"tierS: fetching aa_sequence from NCBI for {len(protein_accs)} protein accessions")
    for i, batch in enumerate(ncbi.chunked(protein_accs, BATCH), 1):
        print(f"  efetch fasta batch {i} ({len(batch)} ids)")
        for k, v in ncbi.efetch_fasta(batch).items():
            record(k, v, "ncbi")

    # PDB entries/chains efetch will not take by accession, but Tier 0 resolved a
    # uid for them -- fetch by uid and map the sequence back onto the accession.
    pdb = [t for t in tier0 if t["resolves"] and t["is_protein"]
           and key_of(t["accession"]) not in seq_by_key and t["uid"] not in ("", "0")]
    if pdb:
        print(f"tierS: fetching {len(pdb)} sequences by uid (PDB etc.)")
        for t in pdb:
            got = ncbi.efetch_fasta_uids([t["uid"]])
            record(key_of(t["accession"]), got.get(t["uid"], ""), "ncbi_uid")

    # Accessions valid outside NCBI: UniProt/TrEMBL (active or UniParc-archived)
    # and MGnify. This is the "check the other databases" step for no_sequence rows.
    external = [t for t in tier0 if key_of(t["accession"]) not in seq_by_key
                and t["accession_class"] in ("uniprot", "mgnify")]
    if external:
        print(f"tierS: querying UniProt/MGnify for {len(external)} non-NCBI accessions")
        for t in external:
            if t["accession_class"] == "uniprot":
                seq, src = extdb.uniprot_sequence(t["accession"])
            else:
                seq, src = extdb.mgnify_sequence(t["accession"])
            record(key_of(t["accession"]), seq, src)
    ext_hits = sum(1 for s in src_by_key.values() if s in ("uniprot", "uniparc", "mgnify"))
    print(f"tierS: {len(seq_by_key)} sequences total ({ext_hits} from non-NCBI databases)")

    def seq_of(accession: str) -> str:
        return seq_by_key.get(key_of(accession), "")

    def src_of(accession: str) -> str:
        return src_by_key.get(key_of(accession), "")

    # --- 2. supplementary index for every cited DOI -------------------------
    all_dois = sorted({d for t in tier0 for d in normalize.split_dois(rows[t["row"]].doi)})
    print(f"tierS: building supplementary index for {len(all_dois)} cited DOIs")
    supp_index: dict[str, dict] = {}
    for doi in all_dois:
        try:
            supp_index[doi] = supp.build_index(doi)
        except Exception as exc:  # never let one paper sink the run
            print(f"  supp {doi}: {exc}", file=sys.stderr)
            supp_index[doi] = {"status": "unavailable", "accessions": [],
                               "sequences": [], "main_text_accessions": []}
    parsed = sum(1 for v in supp_index.values() if v["status"] == "parsed")
    print(f"tierS: {parsed}/{len(all_dois)} DOIs yielded parseable supplementary data")

    # --- 3. duplicate sequences across the sheet ----------------------------
    by_seq: dict[str, list[int]] = collections.defaultdict(list)
    for t in tier0:
        s = seq_of(t["accession"])
        if s:
            by_seq[s].append(t["row"])
    dup_group: dict[int, int] = {}
    dup_partners: dict[int, list[str]] = {}
    gid = 0
    for members in by_seq.values():
        if len(members) > 1:
            gid += 1
            for m in members:
                dup_group[m] = gid
                dup_partners[m] = [rows[o].accession for o in members if o != m]

    # --- 4. per-row adjudication --------------------------------------------
    out = []
    ipg_cache: dict[str, list[dict]] = {}
    for t in tier0:
        r = rows[t["row"]]
        row_seq = seq_of(t["accession"])
        acc_key = normalize.strip_version(t["accession"]).upper()
        dois = normalize.split_dois(r.doi)

        rec = {
            "row": r.index, "pazy_id": r.pazy_id, "enzyme_name": r.enzyme_name,
            "accession": t["accession"], "aa_len": len(row_seq),
            "sequence_source": src_of(t["accession"]),
            "supp_doi": "", "match_kind": "", "recovered_accession": "",
            "recovered_source": "", "n_ipg_members": "",
            "duplicate_group": dup_group.get(r.index, ""),
            "duplicate_partners": "; ".join(dup_partners.get(r.index, [])),
        }

        # supplementary match, strongest evidence first across all cited DOIs
        match_kind, match_doi = "", ""
        for doi in dois:
            idx = supp_index.get(doi, {})
            if row_seq and _seq_match(row_seq, set(idx.get("sequences", []))):
                match_kind, match_doi = "sequence", doi
                break
            if acc_key in set(idx.get("accessions", [])):
                match_kind, match_doi = "accession", doi
            elif not match_kind and acc_key in set(idx.get("main_text_accessions", [])):
                match_kind, match_doi = "accession_maintext", doi

        # IPG recovery for suppressed / dead protein records
        if t["resolves"] and t["is_protein"] and not t["is_current"]:
            if t["accession"] not in ipg_cache:
                try:
                    ipg_cache[t["accession"]] = ncbi.efetch_ipg(t["accession"])
                except Exception as exc:
                    print(f"  ipg {t['accession']}: {exc}", file=sys.stderr)
                    ipg_cache[t["accession"]] = []
            members = ipg_cache[t["accession"]]
            rec["n_ipg_members"] = len(members)
            preferred = _pick_replacement(members, acc_key)
            if preferred:
                rec["recovered_accession"] = preferred["protein"]
                rec["recovered_source"] = preferred["source"]

        if match_kind == "sequence":
            flag = "seq_in_supp"
        elif match_kind == "accession":
            flag = "acc_in_supp"
        elif match_kind == "accession_maintext":
            flag = "acc_in_maintext"
        elif not row_seq:
            flag = "no_sequence"
        elif any(supp_index.get(d, {}).get("status") == "parsed" for d in dois):
            flag = "supp_no_match"
        else:
            flag = "supp_unavailable"

        rec["tierS_flag"] = flag
        rec["supp_doi"] = match_doi
        rec["match_kind"] = match_kind
        out.append(rec)

    # --- 5. write outputs ---------------------------------------------------
    dataset.write_json("tierS.json", out)
    dataset.write_csv("tierS.csv", FIELDS, out)
    _write_enriched(rows, seq_of)

    counts = collections.Counter(r["tierS_flag"] for r in out)
    src_counts = collections.Counter(r["sequence_source"] for r in out if r["sequence_source"])
    confirmed = sum(counts[k] for k in ("seq_in_supp", "acc_in_supp"))
    recovered = sum(1 for r in out if r["recovered_accession"])
    dups = sum(1 for r in out if r["duplicate_group"])
    print(f"\ntierS: {len(out)} rows")
    for k, v in counts.most_common():
        print(f"  {k:<18} {v}")
    print(f"\ntierS: aa_sequence obtained for {sum(src_counts.values())} rows by source:")
    for k, v in src_counts.most_common():
        print(f"  {k:<10} {v}")
    print(f"\ntierS: {confirmed} rows confirmed via supplementary material "
          f"(breaks Tier 2's citation-graph ceiling)")
    print(f"tierS: {recovered} suppressed records got a recovered accession")
    print(f"tierS: {dups} rows share an identical sequence with another row")
    return 0


def _pick_replacement(members: list[dict], query_key: str) -> dict | None:
    """Choose a live replacement from an Identical Protein Group.

    Prefer a RefSeq accession (curated, stable), then any accession that differs
    from the suppressed one. The IPG guarantees sequence identity, so any member
    is a valid re-pointer; RefSeq is simply the best citation target.
    """
    candidates = [m for m in members
                  if normalize.strip_version(m["protein"]).upper() != query_key]
    if not candidates:
        return None
    refseq = [m for m in candidates if "refseq" in m["source"].lower()]
    return (refseq or candidates)[0]


def _write_enriched(rows: dict, seq_of) -> None:
    """Re-emit the input sheet with aa_sequence filled from NCBI."""
    records = []
    for r in sorted(rows.values(), key=lambda x: x.index):
        records.append({
            "enzyme_name": r.enzyme_name, "accession": r.accession,
            "direct_or_indirect_activity": r.direct_or_indirect_activity,
            "pazy_id": r.pazy_id, "doi": r.doi,
            "aa_sequence": r.aa_sequence or seq_of(r.accession),
            "organism": r.organism, "ncbi_taxonomy_id": r.ncbi_taxonomy_id,
        })
    dataset.write_csv(
        "accession-validation-enriched.csv",
        ["enzyme_name", "accession", "direct_or_indirect_activity", "pazy_id",
         "doi", "aa_sequence", "organism", "ncbi_taxonomy_id"],
        records,
    )


if __name__ == "__main__":
    raise SystemExit(main())
