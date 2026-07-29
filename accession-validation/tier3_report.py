"""Tier 3 -- adjudication queue.

Merges Tiers 0-2 into one confidence label per row and emits a review sheet that
puts the record title, organism, linked DOIs, enzyme_name and cited DOI side by
side, ordered so the most actionable rows sit at the top.

On the plan's open decision (single score vs independent flags): we keep the tier
flags independent and *derive* a label from them. A weighted numeric score would
blend a hard fact (the accession does not resolve) with a soft heuristic (the
title is generically named) into one number that no curator could act on. The
label says what to do; the flags say why.
"""

from __future__ import annotations

import collections
import csv

import dataset
import normalize

# Ordered most- to least-urgent. First matching rule wins. Each predicate takes
# the five per-row tier records (t0, t1, t2, ts, to).
CONFIDENCE_RULES = [
    # (label, predicate, what the curator should do)
    # A UniProt/TrEMBL or MGnify id is a *valid* identifier in its own database,
    # just not an NCBI one. Telling a curator to re-source it would be wrong.
    ("external_database",
     lambda t0, t1, t2, ts, to: t0["accession_class"] in ("uniprot", "mgnify") and not t0["resolves"],
     "Valid identifier, wrong database. Fetch the sequence from UniProt/MGnify, "
     "or map it to an NCBI protein accession."),
    ("broken_accession",
     lambda t0, t1, t2, ts, to: not t0["resolves"],
     "Accession resolves to nothing in NCBI. Re-source from the paper."),
    ("wrong_molecule",
     lambda t0, t1, t2, ts, to: t0["resolves"] and not t0["is_protein"],
     "Accession names a nucleotide record. Replace with its protein product."),
    # Suppressed/dead stays a review item even when Tier S confirms the sequence:
    # the biology may be right but the *accession* must still be re-pointed. The
    # recovered_accession field carries the IPG replacement when one was found.
    ("obsolete_record",
     lambda t0, t1, t2, ts, to: t0["resolves"] and t0["is_protein"] and not t0["is_current"],
     "NCBI record is suppressed/dead. Re-map to the current accession "
     "(see recovered_accession)."),
    # A cross-clan family conflict is only trusted when Tier S does NOT find the
    # exact sequence in the cited paper. If the paper's SI reports this very
    # sequence, the 'conflict' is an annotation artefact, not a wrong protein.
    ("family_conflict",
     lambda t0, t1, t2, ts, to: t1["tier1_flag"] == "family_conflict"
     and ts["tierS_flag"] not in ("seq_in_supp", "acc_in_supp"),
     "Record encodes a different enzyme class than enzyme_name implies. Likely wrong protein."),
    ("confirmed",
     lambda t0, t1, t2, ts, to: t2["tier2_flag"] == "doi_match",
     "Cited paper is linked to this accession in PubMed. High confidence."),
    # Tier S -- the citation-graph ceiling breaker. The cited paper's own
    # supplementary material reports this accession/sequence, which PubMed's
    # deposit-only link graph can never show for a mining/screening paper.
    ("confirmed_supp",
     lambda t0, t1, t2, ts, to: ts["tierS_flag"] in ("seq_in_supp", "acc_in_supp"),
     "Cited paper's supplementary material reports this accession/sequence. High confidence."),
    ("title_unmapped",
     lambda t0, t1, t2, ts, to: t1["tier1_flag"] == "title_unmapped",
     "Record is annotated as a function unrelated to enzyme_name "
     "(e.g. amidase vs glutamyl-tRNA amidotransferase). Eyeball it."),
    # Tier O -- the two sources disagree on the source organism at genus level.
    # Only reached by rows Tier 2/S could not confirm (confirmations are matched
    # above), so it surfaces a wrong-protein / mis-deposit candidate that would
    # otherwise have sat silently in `unconfirmed`. The PAZy taxid is itself a
    # claim, so this is "sources conflict", not proof.
    ("organism_mismatch",
     lambda t0, t1, t2, ts, to: to["tierO_flag"] == "mismatch",
     "NCBI record's source organism disagrees with PAZy at genus level "
     "(see organism_detail). Likely wrong protein or mis-deposit; verify."),
    ("name_supported",
     lambda t0, t1, t2, ts, to: t1["tier1_flag"] == "concordant",
     "No DOI link, but record title matches enzyme family. Medium confidence."),
    ("unverifiable_title",
     lambda t0, t1, t2, ts, to: t1["tier1_flag"] == "uninformative_title",
     "Record has no functional annotation (hypothetical protein). Cannot confirm by name."),
    ("fetch_error",
     lambda t0, t1, t2, ts, to: t2["tier2_flag"] == "elink_failed",
     "PubMed link lookup failed. Re-run before drawing any conclusion."),
]

FALLBACK = ("unconfirmed",
            "No DOI link and enzyme_name implies no family. Needs manual adjudication.")

# Rows a curator must look at, most urgent first.
REVIEW_LABELS = [
    "broken_accession", "wrong_molecule", "provenance_mismatch", "obsolete_record",
    "family_conflict", "title_unmapped", "organism_mismatch", "external_database",
    "fetch_error", "unverifiable_title", "unconfirmed",
]

REVIEW_FIELDS = [
    "priority", "confidence", "auto_confidence", "action",
    "curator_verdict", "curator_note",
    "pazy_id", "enzyme_name", "accession",
    "record_title", "organism", "cited_doi", "linked_dois",
    "accession_class", "tier0_pass", "tier1_flag", "tier2_flag", "tierS_flag",
    "tierO_flag", "organism_detail",
    "sequence_source", "supp_doi", "recovered_accession", "duplicate_partners",
    "tierA_flag", "remap_accession", "remap_identity", "remap_organism_flag",
    "remap_author_overlap", "remap_author_shared",
    "name_families", "title_families",
]

# Human curation overrides. Each row a curator has adjudicated is keyed by its
# versionless accession. verdict `correct` clears the row out of the review queue
# (relabelled `verified`, its automated label preserved in auto_confidence);
# verdict `note` leaves the row in the queue but attaches the caveat so the next
# reviewer sees it. This is how a human's decision beats the heuristic without
# loosening any tier's logic (which would suppress genuine future flags).
VERDICTS_PATH = dataset.BASE / "curation" / "verdicts.csv"


def load_verdicts() -> dict[str, dict]:
    if not VERDICTS_PATH.exists():
        return {}
    out: dict[str, dict] = {}
    with open(VERDICTS_PATH, newline="", encoding="utf-8-sig") as fh:
        for rec in csv.DictReader(fh):
            acc = (rec.get("accession") or "").strip()
            if not acc:
                continue
            out[normalize.strip_version(acc).upper()] = {
                "verdict": (rec.get("verdict") or "").strip().lower(),
                "evidence_accession": (rec.get("evidence_accession") or "").strip(),
                "note": (rec.get("note") or "").strip(),
            }
    return out


def label_row(t0: dict, t1: dict, t2: dict, ts: dict, to: dict) -> tuple[str, str]:
    for name, pred, action in CONFIDENCE_RULES:
        if pred(t0, t1, t2, ts, to):
            return name, action
    return FALLBACK


def main() -> int:
    t0s = {r["row"]: r for r in dataset.read_json("tier0.json")}
    t1s = {r["row"]: r for r in dataset.read_json("tier1.json")}
    t2s = {r["row"]: r for r in dataset.read_json("tier2.json")}
    tSs = {r["row"]: r for r in dataset.read_json("tierS.json")}
    tOs = {r["row"]: r for r in dataset.read_json("tierO.json")}
    tPs = ({normalize.strip_version(r["accession"]).upper(): r
            for r in dataset.read_json("tierP.json")}
           if (dataset.OUT / "tierP.json").exists() else {})
    tAs = ({r["row"]: r for r in dataset.read_json("tierA.json")}
           if (dataset.OUT / "tierA.json").exists() else {})
    verdicts = load_verdicts()
    applied = collections.Counter()

    merged = []
    for row, t0 in t0s.items():
        t1, t2, ts, to = t1s[row], t2s[row], tSs[row], tOs[row]
        confidence, action = label_row(t0, t1, t2, ts, to)
        auto_confidence = confidence

        # Apply a human verdict, if one exists for this accession.
        key = normalize.strip_version(t0["accession"]).upper()
        v = verdicts.get(key)
        curator_verdict = v["verdict"] if v else ""
        curator_note = v["note"] if v else ""
        if v:
            applied[v["verdict"]] += 1
            if v["verdict"] == "correct":
                confidence = "verified"
                action = "Human-verified correct. " + curator_note
            elif v["verdict"] == "verify_cds":
                # Tier P checked the protein against the cited nucleotide record.
                p = tPs.get(key)
                pflag = p["tierP_flag"] if p else ""
                pdetail = p["provenance_detail"] if p else "Tier P not run"
                curator_note = (curator_note + " | " + pdetail).strip(" |")
                if pflag in ("cds_id_match", "cds_seq_match"):
                    confidence = "verified"
                    action = "Provenance auto-verified (NCBI): " + pdetail
                elif pflag == "cds_mismatch":
                    confidence = "provenance_mismatch"
                    action = "Cited nucleotide record does NOT encode this protein: " + pdetail

        # Tier A -- sequence-driven re-mapping. It never re-labels a row; it
        # attaches a ranked replacement target + evidence. For an obsolete record
        # with no IPG replacement, surface Tier A's candidate as the target the
        # curator should re-point to.
        tA = tAs.get(row, {})
        recovered = ts["recovered_accession"] or tA.get("remap_accession", "")

        merged.append({
            "row": row,
            "priority": REVIEW_LABELS.index(confidence) if confidence in REVIEW_LABELS else 99,
            "remap_strong": tA.get("tierA_flag") == "remap_strong",
            "confidence": confidence,
            "auto_confidence": auto_confidence,
            "action": action,
            "curator_verdict": curator_verdict,
            "curator_note": curator_note,
            "pazy_id": t0["pazy_id"],
            "enzyme_name": t0["enzyme_name"],
            "accession": t0["accession"],
            "record_title": t0["title"],
            "organism": t0["organism"],
            "cited_doi": t2["row_dois"],
            "linked_dois": t2["linked_dois"],
            "accession_class": t0["accession_class"],
            "resolution_note": t0["resolution_note"],
            "tier0_pass": t0["tier0_pass"],
            "tier1_flag": t1["tier1_flag"],
            "tier2_flag": t2["tier2_flag"],
            "tierS_flag": ts["tierS_flag"],
            "tierO_flag": to["tierO_flag"],
            "organism_detail": to["organism_detail"],
            "sequence_source": ts["sequence_source"],
            "supp_doi": ts["supp_doi"],
            "recovered_accession": recovered,
            "duplicate_partners": ts["duplicate_partners"],
            "tierA_flag": tA.get("tierA_flag", ""),
            "remap_accession": tA.get("remap_accession", ""),
            "remap_identity": tA.get("remap_identity", ""),
            "remap_organism_flag": tA.get("remap_organism_flag", ""),
            "remap_author_overlap": tA.get("remap_author_overlap", ""),
            "remap_author_shared": tA.get("remap_author_shared", ""),
            "name_families": t1["name_families"],
            "title_families": t1["title_families"],
        })

    # Within a label, a Tier A `remap_strong` row is the most actionable (a live,
    # organism- and author-corroborated target is ready), so float it up.
    merged.sort(key=lambda r: (r["priority"], not r["remap_strong"], r["enzyme_name"].lower()))
    dataset.write_csv("validation_full.csv", ["row"] + REVIEW_FIELDS + ["resolution_note"], merged)

    queue = [r for r in merged if r["confidence"] in REVIEW_LABELS]
    dataset.write_csv("review_queue.csv", REVIEW_FIELDS, queue)
    dataset.write_json("summary.json", {
        "total_rows": len(merged),
        "confidence": dict(collections.Counter(r["confidence"] for r in merged)),
        "tier0": dict(collections.Counter(r["tier0_pass"] for r in merged)),
        "tier1": dict(collections.Counter(r["tier1_flag"] for r in merged)),
        "tier2": dict(collections.Counter(r["tier2_flag"] for r in merged)),
        "tierS": dict(collections.Counter(r["tierS_flag"] for r in merged)),
        "tierO": dict(collections.Counter(r["tierO_flag"] for r in merged)),
        "tierA": dict(collections.Counter(r["tierA_flag"] for r in merged if r["tierA_flag"])),
        "curator_verdicts_applied": dict(applied),
        "review_queue_size": len(queue),
    })

    counts = collections.Counter(r["confidence"] for r in merged)
    print(f"tier3: {len(merged)} rows -> out/validation_full.csv")
    if applied:
        print(f"tier3: applied {sum(applied.values())} curator verdict(s): "
              f"{dict(applied)} ({counts.get('verified', 0)} cleared from queue)")
    print(f"tier3: {len(queue)} rows need review -> out/review_queue.csv\n")
    order = ["verified"] + [n for n, _, _ in CONFIDENCE_RULES] + [FALLBACK[0]]
    for name in order:
        if counts.get(name):
            mark = "!" if name in REVIEW_LABELS[:4] else " "
            print(f" {mark} {name:<20} {counts[name]:>4}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
