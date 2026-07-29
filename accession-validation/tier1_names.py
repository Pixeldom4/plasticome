"""Tier 1 -- name concordance.

Compares each resolved record's title against `enzyme_name`. Surfaces candidates
for review; it does not issue verdicts.

The dominant complication in this dataset: 430 of 446 enzyme_names are opaque lab
codes (PET5, PpCutA, LCC) that contain no family word at all. A literal fuzzy
match of enzyme_name against a record title is therefore near-useless. We instead
extract a canonical enzyme family from each side -- using affix conventions on the
code side (`...ase` stems, Cut->cutinase, Est->esterase, Lip->lipase) -- and
compare families. When a name implies no family, Tier 1 abstains and Tier 2 decides.
"""

from __future__ import annotations

import collections

import dataset
import normalize

FIELDS = [
    "row", "pazy_id", "enzyme_name", "accession", "tier0_pass",
    "clean_title", "organism", "name_families", "title_families", "shared_families",
    "family_agreement", "name_informative", "uninformative_title", "title_unmapped",
    "token_jaccard", "fuzzy_ratio", "tier1_flag",
]


def classify(sc: dict, tier0_pass: bool) -> str:
    """Reduce the concordance signals to one triage label."""
    if not tier0_pass:
        return "skipped_tier0"
    agreement = sc["family_agreement"]
    if agreement == "match":
        return "concordant"
    if agreement == "conflict":
        return "family_conflict"          # cross-clan: strongest name-side signal
    if agreement == "related":
        return "related_family"           # within-clan drift: weak, usually benign
    if agreement == "generic":
        return "generic_title"            # superfamily label: consistent, unconfirmed
    if sc["uninformative_title"]:
        # e.g. "hypothetical protein BSQ33_03270" -- resolving proves nothing.
        return "uninformative_title"
    if sc["title_unmapped"] and sc["name_informative"]:
        # Name implies a family; the title names some *other* function we do not
        # recognise (e.g. Amidase -> "glutamyl-tRNA amidotransferase"). Needs eyes.
        return "title_unmapped"
    if not sc["name_informative"]:
        return "name_uninformative"
    return "indeterminate"


def main() -> int:
    rows = {r.index: r for r in dataset.load_rows()}
    tier0 = dataset.read_json("tier0.json")

    out = []
    for t0 in tier0:
        r = rows[t0["row"]]
        sc = normalize.score_concordance(r.enzyme_name, t0["title"])
        flag = classify(sc, t0["tier0_pass"])
        out.append({
            "row": r.index, "pazy_id": r.pazy_id, "enzyme_name": r.enzyme_name,
            "accession": r.accession, "tier0_pass": t0["tier0_pass"],
            **{k: (", ".join(v) if isinstance(v, list) else v) for k, v in sc.items()},
            "tier1_flag": flag,
        })

    dataset.write_json("tier1.json", out)
    dataset.write_csv("tier1.csv", FIELDS, out)

    counts = collections.Counter(r["tier1_flag"] for r in out)
    print(f"tier1: scored {len(out)} rows")
    for k, v in counts.most_common():
        print(f"  {k:<22} {v}")

    conflicts = [r for r in out if r["tier1_flag"] == "family_conflict"]
    if conflicts:
        print(f"\ntier1: {len(conflicts)} family conflicts (strongest wrong-protein signal)")
        for r in conflicts[:12]:
            print(f"  {r['enzyme_name']:<16} {r['accession']:<18} "
                  f"name={r['name_families']:<18} title={r['title_families']:<22} "
                  f"| {r['clean_title'][:38]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
