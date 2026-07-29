"""Tier P -- protein<->nucleotide provenance.

Some screening papers report a source enzyme by its **nucleotide** accession (a
gene / cds entry) while the sheet carries the **protein** accession. That is not a
mismatch on its own -- the protein is usually the CDS product of that nucleotide
record -- but it cannot be confirmed from the paper text. This tier confirms the
link against NCBI: it fetches the cited nucleotide record's coding sequences and
checks whether the sheet's protein is among them, by accession or by sequence.

Driven by `curation/verdicts.csv` rows with `verdict=verify_cds` and the cited
nucleotide in `evidence_accession`. Flags (`tierP_flag`):
  cds_id_match    the protein accession is a `protein_id` of the nucleotide record
  cds_seq_match   the protein's sequence equals one of the record's CDS translations
  cds_mismatch    the record has CDS(s) but none is this protein -- a real problem
  cds_no_cds      the nucleotide record yielded no CDS translations
  cds_unavailable the nucleotide record could not be fetched
"""

from __future__ import annotations

import csv

import dataset
import ncbi
import normalize


def _norm(acc: str) -> str:
    return normalize.strip_version(acc or "").upper()


def load_enriched_sequences() -> dict[str, str]:
    """accession (versionless, upper) -> aa_sequence, from the Tier S enriched sheet."""
    path = dataset.OUT / "accession-validation-enriched.csv"
    out: dict[str, str] = {}
    if not path.exists():
        return out
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for rec in csv.DictReader(fh):
            acc = _norm(rec.get("accession", ""))
            seq = (rec.get("aa_sequence") or "").strip().upper()
            if acc and seq:
                out[acc] = seq
    return out


def check(protein_acc: str, nuc_acc: str, seq: str) -> tuple[str, str, str]:
    """Return (flag, matched_protein_id, detail) for one protein<->nucleotide pair."""
    try:
        cds = ncbi.efetch_cds_aa(nuc_acc)
    except RuntimeError as exc:
        return "cds_unavailable", "", f"could not fetch {nuc_acc}: {exc}"
    if not cds:
        return "cds_no_cds", "", f"{nuc_acc} returned no CDS translations"

    want = _norm(protein_acc)
    for c in cds:
        if _norm(c["protein_id"]) == want:
            return "cds_id_match", c["protein_id"], f"{protein_acc} is a CDS of {nuc_acc}"
    if seq:
        for c in cds:
            if c["seq"] == seq:
                return ("cds_seq_match", c["protein_id"],
                        f"sequence matches CDS {c['protein_id']} of {nuc_acc}")
    ids = ", ".join(c["protein_id"] for c in cds if c["protein_id"]) or "(unnamed CDS)"
    return "cds_mismatch", "", f"{protein_acc} not among {nuc_acc} CDS products: {ids}"


def main() -> int:
    # Reuse Tier 3's verdict loader so the schema stays in one place.
    import tier3_report
    verdicts = tier3_report.load_verdicts()
    pairs = [(acc, v["evidence_accession"]) for acc, v in verdicts.items()
             if v["verdict"] == "verify_cds" and v.get("evidence_accession")]
    print(f"tierP: {len(pairs)} protein<->nucleotide pair(s) to verify")

    seqs = load_enriched_sequences()
    t0 = {_norm(r["accession"]): r for r in dataset.read_json("tier0.json")}

    out = []
    for prot_key, nuc_acc in pairs:
        rec = t0.get(prot_key, {})
        # prot_key is versionless; recover the sheet's exact accession for display.
        prot_acc = rec.get("accession", prot_key)
        flag, matched, detail = check(prot_acc, nuc_acc, seqs.get(prot_key, ""))
        print(f"  {prot_acc} <- {nuc_acc}: {flag} ({detail})")
        out.append({
            "row": rec.get("row", ""),
            "accession": prot_acc,
            "enzyme_name": rec.get("enzyme_name", ""),
            "nuccore_accession": nuc_acc,
            "tierP_flag": flag,
            "matched_protein_id": matched,
            "provenance_detail": detail,
        })

    out.sort(key=lambda r: str(r["accession"]))
    dataset.write_json("tierP.json", out)
    dataset.write_csv(
        "tierP.csv",
        ["row", "accession", "enzyme_name", "nuccore_accession",
         "tierP_flag", "matched_protein_id", "provenance_detail"],
        out,
    )
    confirmed = sum(1 for r in out if r["tierP_flag"] in ("cds_id_match", "cds_seq_match"))
    mismatch = sum(1 for r in out if r["tierP_flag"] == "cds_mismatch")
    print(f"\ntierP: {confirmed} confirmed, {mismatch} mismatch, "
          f"{len(out) - confirmed - mismatch} inconclusive")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
