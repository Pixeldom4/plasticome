"""Tier O -- organism concordance (taxid, rank-aware).

Compares what PAZy says the enzyme's source organism is (the `organism` /
`ncbi_taxonomy_id` columns) against what the NCBI accession record *itself*
reports (`taxid`, captured by Tier 0). This is the check doc `02` specified for
the rows Tier 2 and Tier S leave `unconfirmed`: the citation graph is silent and
no supplementary file was reachable, but the two sources' organism claims can
still be cross-examined.

Why taxid, not free-text name:
  Names diverge on reclassification (Ideonella -> Pseudideonella / Piscinibacter),
  strain suffixes, abbreviations and typos; the taxid collapses all of those. And
  because PAZy may hold a species-level taxid while the record holds a strain-level
  one (or vice versa), both are walked to a common rank (species, then genus)
  before a verdict is called.

Flags (`tierO_flag`):
  exact        same taxid (after aka/merge resolution)
  same_species same species-rank ancestor
  same_genus   same genus-rank ancestor
  mismatch     genus differs on both sides -- a wrong-protein / mis-deposit signal
  uncomparable a taxid is missing or unresolvable on either side

Caveat: the PAZy taxid is itself a claim, so `mismatch` means "the two sources
disagree", not ground truth. It is a review prompt, never a verdict.
"""

from __future__ import annotations

import collections

import dataset
import ncbi


def _clean_taxid(value: str) -> str:
    """PAZy taxids arrive as '1547922', sometimes '1547922.0' or with stray text."""
    value = (value or "").strip()
    if value.endswith(".0"):
        value = value[:-2]
    return value if value.isdigit() else ""


def compare(pazy_taxid: str, rec_taxid: str, tax: dict[str, dict]) -> tuple[str, str]:
    """Return (flag, human-readable detail) for one row's two taxids."""
    if not pazy_taxid and not rec_taxid:
        return "uncomparable", "no taxid on either side"
    if not pazy_taxid:
        return "uncomparable", "no PAZy taxid"
    if not rec_taxid:
        return "uncomparable", "no record taxid"

    p, r = tax.get(pazy_taxid), tax.get(rec_taxid)
    if not p or not r:
        which = "PAZy" if not p else "record"
        return "uncomparable", f"{which} taxid did not resolve in NCBI taxonomy"

    # Canonical taxids (aka/merged ids fold onto the current record).
    if p["taxid"] == r["taxid"]:
        return "exact", f"{p['name']} (taxid {p['taxid']})"

    ps, rs = p["ranks"].get("species"), r["ranks"].get("species")
    if ps and ps == rs:
        return "same_species", f"same species (taxid {ps}); {p['name']} vs {r['name']}"

    pg, rg = p["ranks"].get("genus"), r["ranks"].get("genus")
    if pg and pg == rg:
        return "same_genus", f"same genus (taxid {pg}); {p['name']} vs {r['name']}"
    if pg and rg and pg != rg:
        return "mismatch", f"different genus; PAZy {p['name']} vs record {r['name']}"

    return "uncomparable", "no shared species/genus rank to compare"


def main() -> int:
    rows = {r.index: r for r in dataset.load_rows()}
    t0s = {r["row"]: r for r in dataset.read_json("tier0.json")}
    print(f"tierO: {len(t0s)} rows")

    # Gather every taxid needed on both sides, then resolve lineages in one pass.
    pazy_taxids = {i: _clean_taxid(rows[i].ncbi_taxonomy_id) for i in t0s}
    rec_taxids = {row: str(t0s[row].get("taxid") or "").strip() for row in t0s}
    wanted = sorted({t for t in (*pazy_taxids.values(), *rec_taxids.values()) if t})
    print(f"tierO: resolving {len(wanted)} distinct taxids in NCBI taxonomy")
    tax = ncbi.efetch_taxonomy(wanted)
    n_resolved = sum(1 for t in wanted if t in tax)
    print(f"tierO: {n_resolved}/{len(wanted)} taxids resolved to a lineage")

    out = []
    for row, t0 in t0s.items():
        ptx, rtx = pazy_taxids[row], rec_taxids[row]
        flag, detail = compare(ptx, rtx, tax)
        out.append({
            "row": row,
            "pazy_id": t0["pazy_id"],
            "enzyme_name": t0["enzyme_name"],
            "accession": t0["accession"],
            "pazy_organism": rows[row].organism,
            "pazy_taxid": ptx,
            "record_organism": t0.get("organism", ""),
            "record_taxid": rtx,
            "tierO_flag": flag,
            "organism_detail": detail,
        })

    out.sort(key=lambda r: r["row"])
    dataset.write_json("tierO.json", out)
    dataset.write_csv(
        "tierO.csv",
        ["row", "pazy_id", "enzyme_name", "accession", "pazy_organism", "pazy_taxid",
         "record_organism", "record_taxid", "tierO_flag", "organism_detail"],
        out,
    )

    counts = collections.Counter(r["tierO_flag"] for r in out)
    print(f"\ntierO summary: {len(out)} rows")
    for flag in ("exact", "same_species", "same_genus", "mismatch", "uncomparable"):
        if counts.get(flag):
            print(f"  {flag:<13} {counts[flag]:>4}")
    mism = [r for r in out if r["tierO_flag"] == "mismatch"]
    if mism:
        print("\n  genus-level mismatches (wrong-protein candidates):")
        for r in mism[:12]:
            print(f"    {r['accession']:<18} {r['enzyme_name']:<16} {r['organism_detail'][:60]}")
        if len(mism) > 12:
            print(f"    ... and {len(mism) - 12} more")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
