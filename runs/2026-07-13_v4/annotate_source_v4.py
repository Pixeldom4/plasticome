#!/usr/bin/env python3
"""Annotate the original plasticome final CSV in place with its v4 component.

Emits the 476 source rows in ORIGINAL order and with the original columns
(enzyme_name, accession, organism, aa_sequence) preserved, plus the v4 partition
columns joined 1:1 on sequence md5 (the input is md5-unique, so 1 row = 1 node).
"""
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent           # v4/
CSV_IN = ROOT / "plasticome.v2.cleaned_plasticome_final.csv"
ASSIGN = ROOT / "component_assignment_v4_2026-07-13.csv"
OUT = ROOT / "plasticome_final_v4_components_2026-07-13.csv"


def md5_of(seq):
    return hashlib.md5(re.sub(r"[^A-Za-z]", "", seq).upper().encode()).hexdigest()


def main():
    by_md5 = {r["sequence_md5"]: r for r in csv.DictReader(ASSIGN.open())}
    rows = list(csv.DictReader(CSV_IN.open()))
    add = ["component_id", "size_rank", "v1_component", "cath",
           "sequence_length", "node_id", "sequence_md5"]
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["enzyme_name", "accession", "organism", "aa_sequence"] + add)
        for r in rows:
            a = by_md5[md5_of(r["aa_sequence"])]
            w.writerow([r["enzyme_name"], r["accession"], r["organism"],
                        r["aa_sequence"]] + [a[c] for c in add])
    print(f"wrote {OUT.name}: {len(rows)} rows, "
          f"{len({by_md5[md5_of(r['aa_sequence'])]['component_id'] for r in rows})} components")


if __name__ == "__main__":
    main()
