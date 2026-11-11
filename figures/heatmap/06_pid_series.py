#!/usr/bin/env python3
"""
Step 6: the pairwise heatmap of step 3 on centroids at a series of identities.

For each identity the 609 records are clustered with the project's own
lib/02 clustering/cluster_reference_seeded.py (DIAMOND engine, its defaults:
member covered 90 percent by its representative, the 13 curated seeds kept as
representatives). The representatives are then drawn against each other by
step 3.

Everything goes to pid-series/:
  idNN/clusters.tsv                      the cluster table at NN percent
  idNN/figure_heatmap_sequences_idNN.*   the heatmap
  idNN/sequence_order_idNN.tsv, sequence_matrix_idNN.tsv
  summary.tsv                            one row per identity
  overview.png                           all heatmaps on one sheet
"""

import argparse
import csv
import os
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tree"))
from common import SURFACE, load_rows, read_tsv, write_tsv  # noqa: E402

REPO = os.path.join(HERE, "..", "..")
CLUSTER = os.path.join(REPO, "lib", "02 clustering", "cluster_reference_seeded.py")
SEEDS = os.path.join(REPO, "source-data", "plasticome-curated-seeds",
                     "plasticome-curated-seeds.tsv")
OUT = os.path.join(HERE, "pid-series")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", default="90,80,70,60,50,40,30")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    ids = [int(x) for x in a.ids.split(",")]

    os.makedirs(OUT, exist_ok=True)
    rows = load_rows()
    union = os.path.join(OUT, "union.tsv")
    write_tsv(union, ["plasticome_id", "enzyme_name", "accession", "pazy_id",
                      "source", "aa_sequence"],
              [[p, r["enzyme_name"], r["accession"], r["pazy_id"], r["source"],
                r["aa_sequence"].strip()] for p, r in rows.items()])

    summary = []
    for i in ids:
        d = os.path.join(OUT, f"id{i}")
        os.makedirs(d, exist_ok=True)
        clusters = os.path.join(d, "clusters.tsv")
        log = os.path.join(d, "clustering.log")
        with open(log, "w") as fh:
            subprocess.run([sys.executable, CLUSTER, union, SEEDS, "--id",
                            f"{i / 100:.2f}", "-o", clusters, "--workdir",
                            os.path.join(d, "work"), "--threads", str(a.threads)],
                           check=True, stdout=fh, stderr=subprocess.STDOUT)
        k = len(read_tsv(clusters))
        subprocess.run([sys.executable, os.path.join(HERE, "03_sequence_heatmap.py"),
                        "--centroids", clusters, "--out-dir", d, "--tag", f"_id{i}",
                        "--title", f"{k} centroids at {i}% identity"],
                       check=True, stdout=subprocess.DEVNULL)

        order = read_tsv(os.path.join(d, f"sequence_order_id{i}.tsv"))
        with open(os.path.join(d, f"sequence_matrix_id{i}.tsv")) as fh:
            mat = [r[1:] for r in list(csv.reader(fh, delimiter="\t"))[1:]]
        n = len(order)
        vals = [mat[r][c] for r in range(n) for c in range(r)]
        al = np.array([float(v) for v in vals if v != "NA"])
        fold = [o["fold"] for o in order]
        top = max(set(fold), key=fold.count)
        summary.append([i, k, n, fold.count(top), f"{100 * fold.count(top) / n:.0f}",
                        len(set(fold)), len(set(o["component"] for o in order)),
                        len(vals), len(al), f"{100 * len(al) / len(vals):.1f}",
                        f"{np.median(al):.1f}" if len(al) else ""])
        print(f"id {i}: {k} clusters, {n} drawn, {fold.count(top)} {top} "
              f"({100 * fold.count(top) / n:.0f}%), "
              f"{100 * len(al) / len(vals):.1f}% of pairs aligned")
    write_tsv(os.path.join(OUT, "summary.tsv"),
              ["identity", "clusters", "centroids_drawn", "alpha_beta_hydrolase",
               "hydrolase_percent", "folds", "ssn_components", "pairs",
               "pairs_aligned", "aligned_percent", "median_identity_aligned"],
              summary)

    # one sheet with every heatmap
    cols = 2
    nrow = -(-len(ids) // cols)
    fig, axes = plt.subplots(nrow, cols, figsize=(cols * 9, nrow * 6.6))
    for ax in axes.ravel():
        ax.axis("off")
    for ax, i in zip(axes.ravel(), ids):
        ax.imshow(plt.imread(os.path.join(
            OUT, f"id{i}", f"figure_heatmap_sequences_id{i}.png")))
    fig.patch.set_facecolor(SURFACE)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01,
                        wspace=0.02, hspace=0.02)
    fig.savefig(os.path.join(OUT, "overview.png"), dpi=200, facecolor=SURFACE)


if __name__ == "__main__":
    main()
