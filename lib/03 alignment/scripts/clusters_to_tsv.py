#!/usr/bin/env python3
"""Adapter: cluster-centroids OR union TSV -> curated-TSV shape step1_nodes.py consumes.

The usual input is the one-row-per-cluster table written by
`lib/02 clustering/cluster_reference_seeded.py` (columns `cluster_id, size,
rep_*, member_*, rep_aa_sequence`). Each cluster is already the de-duplicated
unit we want to align, so this script does no biology -- it only re-expresses
each row's REPRESENTATIVE as one row with the six columns step1_nodes.py reads,
then the unchanged alignment pipeline runs on top.

    cluster_id + rep_label -> identifier      (see --id-mode)
    rep_enzyme_name        -> enzyme_name
    rep_pazy_id            -> pazy_id
    rep_accession          -> accession
    (not carried)          -> organism        (always blank)
    rep_aa_sequence        -> aa_sequence

Only the representative is carried forward. The cluster's other members are
recoverable from the original TSV via `annotate_clusters.py`, which joins the
partition back onto every cluster row.

The UNION table from step 1 of the pipeline (`plasticome_id, enzyme_name,
accession, pazy_id, aa_sequence, source`) is also accepted, and the shape is
detected from the header -- that is the align-BEFORE-cluster arrangement, where
every union row is its own node and clustering never runs. Same six output
columns, so steps 1-3 are untouched:

    plasticome_id + accession -> identifier   (PL0001|WP_054022242.1)
    enzyme_name / pazy_id / accession / aa_sequence -> themselves
    (not carried)             -> organism     (always blank)

`identifier` is not cosmetic: step1_nodes.py and step23_graph.py order nodes by
the FIRST NUMBER in it (representative pick on md5-collapse, canonical label,
C### numbering). The clustering script's `rep_label` is `S<nnn>|<acc>` for seeds
and `U<nnnn>|<acc>` for union rows, so its numbers COLLIDE (S001 and U0001 both
read as 1). Hence the default `--id-mode cluster`, which prefixes the globally
unique `cluster_id` and keeps the original label visible behind it:

    cluster (default)  CL0001|S001|WP_054022242.1     unique, ordered by cluster_id
    label              S001|WP_054022242.1            verbatim; collisions reported

Component/edge COUNTS are label-independent either way -- the mode only affects
which node represents a tie and how components are numbered and named.
`--id-mode` applies to cluster input only; `plasticome_id` is already globally
unique, so union rows always take the numbered form.
"""
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

IN_COLS = ["cluster_id", "rep_label", "rep_enzyme_name", "rep_accession",
           "rep_pazy_id", "rep_aa_sequence"]
UNION_COLS = ["plasticome_id", "enzyme_name", "accession", "pazy_id", "aa_sequence"]
OUT_COLS = ["identifier", "enzyme_name", "pazy_id", "accession", "organism", "aa_sequence"]


def pl_num(ident: str) -> int:
    """The number step1_nodes.py/step23_graph.py will read out of an identifier."""
    m = re.search(r"(\d+)", ident or "")
    return int(m.group(1)) if m else 10**9


def clean(s: str) -> str:
    return (s or "").strip()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clusters", type=Path, required=True,
                    help="clusters TSV (one row per cluster) or union TSV (one row per sequence)")
    ap.add_argument("--out", type=Path, required=True, help="curated-shape TSV for step1_nodes.py")
    ap.add_argument("--id-mode", choices=("cluster", "label"), default="cluster",
                    help="how to build `identifier` from cluster input (default: cluster)")
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    with args.clusters.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"[clusters_to_tsv] no data rows in {args.clusters}")

    # Shape from the header: a `plasticome_id` column means the union table (align
    # before cluster), anything else is read as the clusters table.
    shape = "union" if "plasticome_id" in rows[0] else "clusters"
    want = UNION_COLS if shape == "union" else IN_COLS
    missing = [c for c in want if c not in rows[0]]
    if missing:
        sys.exit(f"[clusters_to_tsv] {args.clusters} looks like a {shape} TSV but is "
                 f"missing column(s): {', '.join(missing)}")

    out_rows, n_blank_seq = [], 0
    for i, r in enumerate(rows, start=1):
        if shape == "union":
            pid = clean(r["plasticome_id"])
            n = int(pid) if pid.isdigit() else i
            ident = f"PL{n:04d}|{clean(r['accession'])}"
            src = {"enzyme_name": "enzyme_name", "pazy_id": "pazy_id",
                   "accession": "accession", "aa_sequence": "aa_sequence"}
        else:
            label = clean(r["rep_label"])
            cid = clean(r["cluster_id"])
            cnum = int(cid) if cid.isdigit() else i
            ident = label if args.id_mode == "label" else f"CL{cnum:04d}|{label}"
            src = {"enzyme_name": "rep_enzyme_name", "pazy_id": "rep_pazy_id",
                   "accession": "rep_accession", "aa_sequence": "rep_aa_sequence"}
        seq = re.sub(r"\s", "", r[src["aa_sequence"]] or "")
        if not seq:
            n_blank_seq += 1
        out_rows.append({
            "identifier": ident,
            "enzyme_name": clean(r[src["enzyme_name"]]),
            "pazy_id": clean(r[src["pazy_id"]]),
            "accession": clean(r[src["accession"]]),
            "organism": "",
            "aa_sequence": seq,
        })

    # The ordering key must be unique or the C### numbering stops being stable.
    dupes = [n for n, c in Counter(pl_num(r["identifier"]) for r in out_rows).items() if c > 1]
    if dupes:
        print(f"[clusters_to_tsv] WARNING: {len(dupes)} identifier number(s) shared by >1 "
              f"row (e.g. {sorted(dupes)[:5]}) -- component numbering/representative "
              f"picks fall back to label order. Use --id-mode cluster to avoid this.")

    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    blank = f", {n_blank_seq} blank sequence(s) -- step1 will drop these" if n_blank_seq else ""
    mode = f"id-mode={args.id_mode}" if shape == "clusters" else "id=plasticome_id"
    unit = "clusters" if shape == "clusters" else "union rows"
    print(f"[clusters_to_tsv] {len(out_rows)} {unit} -> {args.out} ({mode}{blank})")


if __name__ == "__main__":
    main()
