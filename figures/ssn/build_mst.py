#!/usr/bin/env python3
"""
Panel a-i as a minimum spanning tree, after the reviewer feedback of
2026-10-27.

Nodes are the 492 distinct sequences of the tree (../tree, all-sequences run),
so a node stands for every record with that sequence and its area is the
number of records. Encodings follow the tree: fill is substrate, shape is
organism, the background blob is the fold in its colour from
../tree/folds.py.

Edges come from two kinds of evidence and are drawn differently so they are
never confused:

  solid grey   sequence. Inside each SSN component (30 percent family), the
               minimum spanning tree on 100 minus pairwise identity, from the
               exhaustive Smith-Waterman pairs of ../heatmap (a pair
               counts when it covers half the shorter sequence). Width and
               darkness follow identity
  dashed       structure. Inside a fold, the components are joined by the
               minimum spanning tree on 1 minus TM-score between their
               representative structures (../tree step 4), drawn between the
               representative sequences in the fold colour. Edge number: TM
  none         between folds. Folds share no detectable homology (tree step
               4), so they are separate trees and the figure is a forest of 16

A component that the coverage rule leaves in pieces (single linkage in the SSN
ran through DIAMOND hits) is joined by its best remaining pair, drawn dotted.

Layout is Graphviz neato with edge length set by distance, each fold packed as
its own connected piece.

Writes out/figure_mst.{svg,pdf,png} and out/mst_edges.tsv.
"""

import collections
import os
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["svg.hashsalt"] = "plasticome-mst"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(HERE, "..", "tree")
sys.path.insert(0, TREE)
os.environ["PLASTICOME_TREE_SET"] = "all"
from common import (INK, MUTED, NODE_STROKE, OUT as TREE_OUT, SUBSTRATE_COLOR, SURFACE,  # noqa: E402
                    WORK as TREE_WORK, load_rows, read_tsv, seed_nodes, write_tsv)
import folds  # noqa: E402

PAIRS = os.path.join(HERE, "..", "heatmap", "work", "pairs.tsv")
OUT = os.path.join(HERE, "out")
COVERAGE = 0.5
SHAPE = {"Bacteria": "D", "Archaea": "^", "Eukaryota": "h"}


def main():
    rows = load_rows()
    tips = {t["tip"]: t for t in read_tsv(os.path.join(TREE_OUT, "tree_tips.tsv"))}
    ids = sorted(tips, key=int)
    ix = {t: i for i, t in enumerate(ids)}
    n = len(ids)
    ident = np.zeros((n, n))
    cov = np.zeros((n, n))
    with open(PAIRS) as fh:
        next(fh)
        for line in fh:
            a, b, idv, _, _, c, *_ = line.split("\t")
            if a in ix and b in ix:
                ident[ix[a], ix[b]] = ident[ix[b], ix[a]] = float(idv)
                cov[ix[a], ix[b]] = cov[ix[b], ix[a]] = float(c)

    edges = []          # (a, b, kind, value)
    comps = collections.defaultdict(list)
    for t in ids:
        comps[tips[t]["component"]].append(t)
    for c, mem in comps.items():
        if len(mem) < 2:
            continue
        idx = [ix[t] for t in mem]
        sub_i, sub_c = ident[np.ix_(idx, idx)], cov[np.ix_(idx, idx)]
        ok = (sub_c >= COVERAGE) & (sub_i > 0)
        d = np.where(ok, 100.0 - sub_i + 1e-3, 0.0)
        mst = minimum_spanning_tree(d).tocoo()
        for i, j, _ in zip(mst.row, mst.col, mst.data):
            edges.append((mem[i], mem[j], "sequence", sub_i[i, j]))
        # join pieces left apart by the coverage rule with their best pair
        parent = list(range(len(mem)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for i, j in zip(mst.row, mst.col):
            parent[find(i)] = find(j)
        while len({find(i) for i in range(len(mem))}) > 1:
            best = None
            for i in range(len(mem)):
                for j in range(len(mem)):
                    if find(i) != find(j) and (best is None or sub_i[i, j] > best[2]):
                        best = (i, j, sub_i[i, j])
            i, j, v = best
            edges.append((mem[i], mem[j], "weak", v))
            parent[find(i)] = find(j)

    # components inside a fold, joined on structure
    struct = {int(r["component"]): r["representative"] for r in
              read_tsv(os.path.join(TREE_WORK, "structures.tsv")) if r["representative"]}
    tm = {}
    for r in read_tsv(os.path.join(TREE_WORK, "homology", "pairs.tsv")):
        # no structure for this component (IS11): TM is NA, kept as -1 and
        # joined last; the edge is labelled "seq" since HHalign placed it
        v = float(r["tm"]) if r["tm"] not in ("NA", "") else -1.0
        tm[(int(r["a"]), int(r["b"]))] = tm[(int(r["b"]), int(r["a"]))] = v
    by_fold = collections.defaultdict(set)
    for t in ids:
        by_fold[folds.fold_key(tips[t])].add(int(tips[t]["component"].lstrip("c")))
    rep_of = {}
    for c, mem in comps.items():
        ci = int(c.lstrip("c"))
        rep_of[ci] = struct.get(ci) if struct.get(ci) in tips else mem[0]
    for k, cs in by_fold.items():
        cs = sorted(cs)
        if len(cs) < 2:
            continue
        d = np.array([[0.0 if a == b else (0.95 if tm.get((a, b), 0.0) < 0 else
                                           1.0 - tm.get((a, b), 0.0) + 1e-3) for b in cs] for a in cs])
        mst = minimum_spanning_tree(d).tocoo()
        for i, j in zip(mst.row, mst.col):
            edges.append((rep_of[cs[i]], rep_of[cs[j]], "structure", tm.get((cs[i], cs[j]), 0.0)))

    # --- layout -------------------------------------------------------------
    dot = ["graph G {", "  graph [overlap=false, pack=24, packmode=\"graph\", splines=false];",
           "  node [shape=point, width=0.05];"]
    for t in ids:
        dot.append(f'  n{t};')
    for a, b, kind, v in edges:
        ln = (100 - v) / 22 + 0.15 if kind != "structure" else 1.4 + 3.0 * (1 - max(v, 0.5))
        dot.append(f"  n{a} -- n{b} [len={ln:.3f}];")
    dot.append("}")
    res = subprocess.run(["neato", "-Tplain"], input="\n".join(dot), capture_output=True,
                         text=True, check=True)
    pos = {}
    for line in res.stdout.splitlines():
        p = line.split()
        if p and p[0] == "node":
            pos[p[1][1:]] = (float(p[2]), float(p[3]))

    # --- draw ---------------------------------------------------------------
    xs = np.array([pos[t] for t in ids])
    span = xs.max(axis=0) - xs.min(axis=0)
    W = 16.0
    H = W * span[1] / span[0]
    fig = plt.figure(figsize=(W + 3.4, H + 0.6))
    ax = fig.add_axes([0.01, 0.01, W / (W + 3.4), 0.98])
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor(SURFACE)
    pad = 0.35
    for k in by_fold:
        pts = np.array([pos[t] for t in ids if folds.fold_key(tips[t]) == k])
        circ = np.array([[np.cos(a), np.sin(a)] for a in np.linspace(0, 2 * np.pi, 16, endpoint=False)])
        cloud = (pts[:, None, :] + pad * circ[None]).reshape(-1, 2)
        hull = cloud[ConvexHull(cloud).vertices]
        ax.fill(hull[:, 0], hull[:, 1], color=folds.tint(k), zorder=0, ec=folds.band(k), lw=0.8)
        top = pts[np.argmax(pts[:, 1])]
        ax.text(top[0], top[1] + pad + 0.12, f"{folds.name(k)} ({len(pts)})", ha="center",
                va="bottom", fontsize=7.5, color=folds.strong(k), fontweight="bold", zorder=9)
    for a, b, kind, v in edges:
        (x0, y0), (x1, y1) = pos[a], pos[b]
        if kind == "sequence":
            g = 0.75 - 0.6 * min(max((v - 30) / 70, 0), 1)
            ax.plot([x0, x1], [y0, y1], color=(g, g, g), lw=0.4 + 2.0 * min(max((v - 30) / 70, 0), 1),
                    zorder=2, solid_capstyle="round")
        elif kind == "weak":
            ax.plot([x0, x1], [y0, y1], color=MUTED, lw=0.6, ls=(0, (1, 1.5)), zorder=2)
        else:
            k = folds.fold_key(tips[a])
            ax.plot([x0, x1], [y0, y1], color=folds.strong(k), lw=0.9, ls=(0, (4, 2.5)), zorder=1.5)
            ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{v:.2f}" if v >= 0 else "HHalign", fontsize=4.5, color=folds.strong(k),
                    ha="center", va="center", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.08", fc=folds.tint(k), ec="none"))
    by = collections.defaultdict(list)
    for t in ids:
        r = rows[tips[t]["records"].split(";")[0]]
        mk = "o" if (r["metagenomic"] or "").lower() == "yes" else SHAPE.get(r["superkingdom"], "o")
        nrec = len(tips[t]["records"].split(";"))
        by[(mk, SUBSTRATE_COLOR[tips[t]["substrate_group"]])].append((*pos[t], nrec))
    for (mk, col), p in by.items():
        p = np.array(p)
        ax.scatter(p[:, 0], p[:, 1], s=14 + 12 * (p[:, 2] - 1) ** 0.85, marker=mk, c=col,
                   edgecolors=NODE_STROKE, linewidths=0.35, zorder=5)
    seeds = {pid: nm for nm, pid in seed_nodes(rows).items()}
    for t in ids:
        nm = next((seeds[r] for r in tips[t]["records"].split(";") if r in seeds), None)
        if nm:
            x, y = pos[t]
            ax.scatter([x], [y], s=60, facecolors="none", edgecolors=INK, linewidths=0.9, zorder=6)
            ax.annotate(nm, (x, y), xytext=(5, 5), textcoords="offset points", fontsize=7,
                        color=INK, zorder=9, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.75))
    ax.margins(0.02)

    lg = fig.add_axes([W / (W + 3.4), 0.05, 3.3 / (W + 3.4), 0.9])
    lg.axis("off")
    h1 = [Line2D([], [], color=(0.15, 0.15, 0.15), lw=2.2, label="sequence, 100% identity"),
          Line2D([], [], color=(0.45, 0.45, 0.45), lw=1.2, label="sequence, 65%"),
          Line2D([], [], color=(0.75, 0.75, 0.75), lw=0.4, label="sequence, 30%"),
          Line2D([], [], color=MUTED, lw=0.6, ls=(0, (1, 1.5)), label="sequence, best pair under\nthe coverage rule"),
          Line2D([], [], color=folds.strong("3.40.50.1820"), lw=0.9, ls=(0, (4, 2.5)),
                 label="structure, between components\nof one fold (number: TM-score)")]
    h2 = [Line2D([], [], marker="s", ls="none", ms=7, mfc=SUBSTRATE_COLOR[s], mec="none", label=s)
          for s in ["PET", "PE", "PA", "PUR", "PLA", "PBAT", "PBS", "PHA", "PC", "NR", "Combi."]]
    h3 = [Line2D([], [], marker=m, ls="none", ms=6, mfc="none", mec=INK, label=l) for l, m in
          [("Bacterial", "D"), ("Archaeal", "^"), ("Eukaryotic", "h"), ("Metagenomic", "o")]]
    y = 1.0
    for hs, title, nc in [(h1, "edges (minimum spanning tree)", 1), (h2, "node colour (substrate)", 3),
                          (h3, "node shape (organism)", 1)]:
        L = lg.legend(handles=hs, title=title, loc="upper left", bbox_to_anchor=(0, y), frameon=False,
                      fontsize=7, title_fontsize=7.5, ncol=nc, alignment="left")
        lg.add_artist(L)
        y -= 0.06 + 0.045 * (len(hs) / nc + 1)
    lg.text(0, y - 0.02, "Node area: records with that sequence.\nBlob: fold, in the tree's colours.\n"
            "No edges between folds: no shared\nhomology was detected (tree step 4).\n"
            f"{n} sequences, {sum(1 for e in edges if e[2] != 'structure')} sequence edges,\n"
            f"{sum(1 for e in edges if e[2] == 'structure')} structure edges, {len(by_fold)} trees.",
            fontsize=7, color=MUTED, va="top", linespacing=1.4)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(os.path.join(OUT, f"figure_mst.{ext}"), dpi=300, facecolor=SURFACE,
                    **({"metadata": {"Date": None}} if ext == "svg" else {}))
    write_tsv(os.path.join(OUT, "mst_edges.tsv"), ["a", "b", "kind", "value"],
              [[a, b, k, f"{v:.3f}"] for a, b, k, v in edges])
    c = collections.Counter(e[2] for e in edges)
    print(f"figure_mst: {n} nodes, edges {dict(c)}, {len(by_fold)} folds")


if __name__ == "__main__":
    main()
