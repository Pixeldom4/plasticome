#!/usr/bin/env python3
"""Join the computed partition back onto the cluster-centroids TSV.

Emits every cluster row in ORIGINAL order with all of its original columns
preserved, plus the partition columns joined on the md5 of `rep_aa_sequence`
(the same node identity step1_nodes.py uses). This is the cluster-table twin of
`annotate_source.py`, which does the same job for the curated PAZy TSV.

Two files come out of one pass:
  --out       the full join: every original column + node_status + the partition
              columns (component_id, size_rank, v1_component, cath, ...)
  --slim-out  the run deliverable ("03 alignment.tsv"): every original column
              plus exactly one new one, component_id

Clusters whose representatives are byte-identical after normalize() collapse to
one node upstream and therefore land in the same component -- that is reported,
not an error. A cluster with no sequence gets empty partition columns and
`node_status = no_sequence`.
"""
import argparse
import csv
import hashlib
import re
import sys
from collections import Counter
from pathlib import Path

ADD_COLS = ["component_id", "size_rank", "v1_component", "cath",
            "sequence_length", "node_id", "sequence_md5"]


def md5_of(seq: str) -> str:
    s = re.sub(r"[^A-Za-z]", "", seq or "").upper()
    return hashlib.md5(s.encode()).hexdigest() if s else ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clusters", type=Path, required=True, help="the original clusters TSV")
    ap.add_argument("--assignment", type=Path, required=True,
                    help="component_assignment_*.csv from step23_graph.py")
    ap.add_argument("--out", type=Path, required=True,
                    help="full join: every original column + node_status + the partition columns")
    ap.add_argument("--slim-out", type=Path, default=None,
                    help="also write the run deliverable: every original column + component_id only")
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.slim_out:
        args.slim_out.parent.mkdir(parents=True, exist_ok=True)

    with args.assignment.open(newline="") as fh:
        by_md5 = {r["sequence_md5"]: r for r in csv.DictReader(fh)}
    with args.clusters.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        src_cols, rows = rdr.fieldnames or [], list(rdr)
    if "rep_aa_sequence" not in src_cols:
        sys.exit(f"[annotate_clusters] {args.clusters} has no rep_aa_sequence column")

    n_assigned = n_missing = 0
    # Tab-separated, like the input -- these rows carry `; `-joined member lists
    # and free-text enzyme names, which read badly as CSV.
    slim_fh = args.slim_out.open("w", newline="") if args.slim_out else None
    try:
        with args.out.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(src_cols + ["node_status"] + ADD_COLS)
            slim = csv.writer(slim_fh, delimiter="\t") if slim_fh else None
            if slim:
                slim.writerow(src_cols + ["component_id"])
            for r in rows:
                a = by_md5.get(md5_of(r["rep_aa_sequence"]))
                src = [r.get(c, "") for c in src_cols]
                if a:
                    n_assigned += 1
                    w.writerow(src + ["assigned"] + [a[c] for c in ADD_COLS])
                else:
                    n_missing += 1
                    status = "no_sequence" if not md5_of(r["rep_aa_sequence"]) else "unassigned"
                    w.writerow(src + [status] + [""] * len(ADD_COLS))
                if slim:
                    slim.writerow(src + [a["component_id"] if a else ""])
    finally:
        if slim_fh:
            slim_fh.close()

    comps, nodes = set(), Counter()
    for r in rows:
        a = by_md5.get(md5_of(r["rep_aa_sequence"]))
        if a:
            comps.add(a["component_id"])
            nodes[a["node_id"]] += 1
    shared = sum(1 for n, c in nodes.items() if c > 1)
    note = f", {shared} node(s) shared by >1 cluster (identical representatives)" if shared else ""
    if n_missing:
        note += f", {n_missing} row(s) not assigned"
    print(f"[annotate_clusters] wrote {args.out.name}: {len(rows)} clusters "
          f"({n_assigned} assigned) -> {len(comps)} components{note}")
    if args.slim_out:
        print(f"[annotate_clusters] wrote {args.slim_out.name}: same rows + component_id")


if __name__ == "__main__":
    main()
