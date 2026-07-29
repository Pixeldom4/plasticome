#!/usr/bin/env python3
"""Annotate the updated cleaned PAZy final TSV in place with its v6 component.

Emits every source row in ORIGINAL order with the original columns preserved, plus
the v6 partition columns joined on sequence md5. md5-duplicate rows all receive the
SAME node's component; any blank-sequence rows (none in this input) would get empty
partition columns and node_status = "no_sequence". Graphed rows get "assigned".
"""
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TSV_IN = ROOT / "cleaned_pazy_final.tsv"
ASSIGN = ROOT / "component_assignment_v6_2026-07-17.csv"
OUT = ROOT / "cleaned_pazy_final_v6_components_2026-07-17.csv"

SRC_COLS = ["identifier", "enzyme_name", "pazy_id", "accession", "organism", "aa_sequence"]


def md5_of(seq):
    s = re.sub(r"[^A-Za-z]", "", seq or "").upper()
    return hashlib.md5(s.encode()).hexdigest() if s else ""


def main():
    by_md5 = {r["sequence_md5"]: r for r in csv.DictReader(ASSIGN.open())}
    rows = list(csv.DictReader(TSV_IN.open(), delimiter="\t"))
    add = ["component_id", "size_rank", "v1_component", "cath",
           "sequence_length", "node_id", "sequence_md5"]
    n_assigned = n_blank = 0
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(SRC_COLS + ["node_status"] + add)
        for r in rows:
            h = md5_of(r["aa_sequence"])
            a = by_md5.get(h)
            if a:
                n_assigned += 1
                w.writerow([r[c] for c in SRC_COLS] + ["assigned"] + [a[c] for c in add])
            else:
                n_blank += 1
                w.writerow([r[c] for c in SRC_COLS] + ["no_sequence"] + [""] * len(add))
    comps = {by_md5[md5_of(r["aa_sequence"])]["component_id"]
             for r in rows if md5_of(r["aa_sequence"]) in by_md5}
    print(f"wrote {OUT.name}: {len(rows)} rows "
          f"({n_assigned} assigned, {n_blank} no_sequence), {len(comps)} components")


if __name__ == "__main__":
    main()
