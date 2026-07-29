#!/usr/bin/env python3
"""Annotate the updated cleaned PAZy final CSV in place with its v5 component.

Emits every source row in ORIGINAL order with the original columns preserved, plus
the v5 partition columns joined on sequence md5. Unlike v4 the input is NOT
md5-unique and has blank-sequence rows:
  * md5-duplicate rows all receive the SAME node's component (many rows -> 1 node).
  * blank-sequence rows (dropped from the graph) get empty partition columns and
    node_status = "no_sequence"; graphed rows get node_status = "assigned".
"""
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_IN = ROOT / "cleaned_pazy_final.csv"
ASSIGN = ROOT / "component_assignment_v5_2026-07-16.csv"
OUT = ROOT / "cleaned_pazy_final_v5_components_2026-07-16.csv"


def md5_of(seq):
    s = re.sub(r"[^A-Za-z]", "", seq or "").upper()
    return hashlib.md5(s.encode()).hexdigest() if s else ""


def main():
    by_md5 = {r["sequence_md5"]: r for r in csv.DictReader(ASSIGN.open())}
    rows = list(csv.DictReader(CSV_IN.open()))
    add = ["component_id", "size_rank", "v1_component", "cath",
           "sequence_length", "node_id", "sequence_md5"]
    n_assigned = n_blank = 0
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["enzyme_name", "accession", "organism", "aa_sequence",
                    "node_status"] + add)
        for r in rows:
            h = md5_of(r["aa_sequence"])
            a = by_md5.get(h)
            if a:
                n_assigned += 1
                w.writerow([r["enzyme_name"], r["accession"], r["organism"],
                            r["aa_sequence"], "assigned"] + [a[c] for c in add])
            else:
                n_blank += 1
                w.writerow([r["enzyme_name"], r["accession"], r["organism"],
                            r["aa_sequence"], "no_sequence"] + [""] * len(add))
    comps = {by_md5[md5_of(r["aa_sequence"])]["component_id"]
             for r in rows if md5_of(r["aa_sequence"]) in by_md5}
    print(f"wrote {OUT.name}: {len(rows)} rows "
          f"({n_assigned} assigned, {n_blank} no_sequence), {len(comps)} components")


if __name__ == "__main__":
    main()
