#!/usr/bin/env python3
"""Annotate the curated PAZy final TSV in place with its recomputed component.

Emits every source row in ORIGINAL order with the original columns preserved, plus
the partition columns joined on sequence md5. md5-duplicate rows all receive the
SAME node's component; any blank-sequence rows (none in this input) would get empty
partition columns and node_status = "no_sequence"; graphed rows get "assigned".

The input carries a pre-existing `component_id` column (a prior assignment); it is
preserved as `prior_component_id` so the freshly computed partition can be compared
against it. The de-novo partition never reads that prior column.
"""
import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402

SRC_COLS = ["identifier", "enzyme_name", "pazy_id", "accession", "organism", "aa_sequence"]
ADD_COLS = ["component_id", "size_rank", "v1_component", "cath",
            "sequence_length", "node_id", "sequence_md5"]


def md5_of(seq):
    s = re.sub(r"[^A-Za-z]", "", seq or "").upper()
    return hashlib.md5(s.encode()).hexdigest() if s else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", type=Path, default=config.INPUT_TSV)
    ap.add_argument("--assignment", type=Path, required=True,
                    help="component_assignment_*.csv from step23_graph.py")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    by_md5 = {r["sequence_md5"]: r for r in csv.DictReader(args.assignment.open())}
    rows = list(csv.DictReader(args.tsv.open(), delimiter="\t"))
    n_assigned = n_blank = 0
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(SRC_COLS + ["prior_component_id", "node_status"] + ADD_COLS)
        for r in rows:
            a = by_md5.get(md5_of(r["aa_sequence"]))
            prior = (r.get("component_id") or "").strip()
            if a:
                n_assigned += 1
                w.writerow([r[c] for c in SRC_COLS] + [prior, "assigned"] + [a[c] for c in ADD_COLS])
            else:
                n_blank += 1
                w.writerow([r[c] for c in SRC_COLS] + [prior, "no_sequence"] + [""] * len(ADD_COLS))
    comps = {by_md5[md5_of(r["aa_sequence"])]["component_id"]
             for r in rows if md5_of(r["aa_sequence"]) in by_md5}
    print(f"[annotate] wrote {args.out.name}: {len(rows)} rows "
          f"({n_assigned} assigned, {n_blank} no_sequence), {len(comps)} components")


if __name__ == "__main__":
    main()
