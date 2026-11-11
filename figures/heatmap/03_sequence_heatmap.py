#!/usr/bin/env python3
"""
Step 3: every sequence against every other, one cell per pair.

The same triangle as step 2 without the grouping: nothing is pooled into a
median and no sequence is left out. The default set is the 90 percent
centroids, 411 sequences and 84,255 pairs. PLASTICOME_TREE_SET=all draws the
492 distinct sequences instead.

Order
  Rows run fold by fold, largest fold first, and inside a fold in the order of
  its tree (../tree, order_in_tree). Each SSN component is one clade of that
  tree, so a component is a run of neighbouring rows and close relatives sit
  next to each other inside it.

Cells
  Identity of the pair from step 1, in the bins of step 2. A pair that does
  not align over half of the shorter sequence (the SSN rule) has no identity
  and is grey: darker when the two share a fold, lighter when they do not.

Writes figure_heatmap_sequences.{svg,pdf,png}, sequence_order.tsv and
sequence_matrix.tsv.

--centroids takes a cluster table from lib/02 clustering/
cluster_reference_seeded.py and draws its representatives instead, so the
same figure can be made at any clustering identity (step 6 does this for a
series). Those figures are ordered by the tree of all sequences, are scaled to
one matrix size, and go to --out-dir with --tag in the file names.
"""

import argparse
import collections
import importlib.util
import os
import sys

os.environ.setdefault("PLASTICOME_TREE_SET", "centroids")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(HERE, "..", "tree")
sys.path.insert(0, TREE)
from common import (INK, MUTED, OUT as TREE_OUT, OUT_ALL, SET, SUFFIX,  # noqa: E402
                    SURFACE, WORK_ALL, load_rows, read_tsv, seed_nodes, write_tsv)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# bins, colours and the pair table are owned by step 2, label spreading by
# the pie, so the three figures cannot drift apart
hm = load(os.path.join(HERE, "02_heatmap.py"), "heatmap_groups")
spread = load(os.path.join(TREE, "08_draw_pie.py"), "draw_pie").spread

OUT = os.path.join(HERE, "out" + SUFFIX)
DPI = 300
CELL = 6 / DPI          # a whole number of pixels, so cells stay even
SIDE = 8.0              # matrix side with --centroids, inches
STRIP = ("#8a8880", "#c9c7bf")     # neighbouring components alternate
FOLD_MIN = 3            # folds under this many sequences share one label


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-seed-labels", action="store_true")
    ap.add_argument("--centroids", help="cluster table; draw its representatives")
    ap.add_argument("--out-dir", help="where to write, with --centroids")
    ap.add_argument("--tag", default="", help="suffix for the output files")
    ap.add_argument("--title", default="", help="a line above the legend")
    a = ap.parse_args()
    tree_out, out = TREE_OUT, OUT
    if a.centroids:
        tree_out = OUT_ALL       # the tree of all sequences
        out = a.out_dir or os.path.join(HERE, "pid-series")

    rows = load_rows()
    full = {rows[t["tip"]]["aa_sequence"].strip().upper(): t["tip"]
            for t in read_tsv(os.path.join(WORK_ALL, "tips.tsv"))}
    tips = read_tsv(os.path.join(tree_out, "tree_tips.tsv"))
    tips.sort(key=lambda t: (t["fold_group"], int(t["order_in_tree"])))
    for t in tips:      # pairs.tsv is keyed by the tips of the full set
        t["key"] = full[rows[t["tip"]]["aa_sequence"].strip().upper()]
    if a.centroids:
        reps = {r["rep_aa_sequence"].strip().upper()
                for r in read_tsv(a.centroids)}
        have = {rows[t["tip"]]["aa_sequence"].strip().upper() for t in tips}
        if reps - have:
            print(f"  {len(reps - have)} representative(s) not in the tree "
                  "(fragments), left out")
        tips = [t for t in tips
                if rows[t["tip"]]["aa_sequence"].strip().upper() in reps]
    order = [t["key"] for t in tips]
    n = len(order)
    m, ix = hm.load_identity(order)
    assert [ix[k] for k in order] == list(range(n))

    seed = {}
    for name, pid in seed_nodes(rows).items():
        for i, t in enumerate(tips):
            if pid in t["records"].split(";"):
                seed[i] = name
    fold_name = {r["fold_group"]: r["fold"] for r in
                 read_tsv(os.path.join(tree_out, "fold_table.tsv"))}
    fold = np.array([t["fold_group"] for t in tips])
    comp = [t["component"] for t in tips]

    # ---------------------------------------------------------------- image
    lo = np.array([b[0] for b in hm.BINS])
    rgba = np.array([to_rgba(b[1]) for b in hm.BINS])
    which = np.clip(np.searchsorted(lo, m, side="right") - 1, 0, len(lo) - 1)
    img = rgba[which]
    same_fold = fold[:, None] == fold[None, :]
    img[(m == 0) & same_fold] = to_rgba(hm.NO_ALN_SAME_FOLD)
    img[(m == 0) & ~same_fold] = to_rgba(hm.NO_ALN_OTHER_FOLD)
    img[np.triu_indices(n)] = (0, 0, 0, 0)

    cell = max(1, int(SIDE * DPI / n)) / DPI if a.centroids else CELL
    side = n * cell
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    rend = fig.canvas.get_renderer()
    fs = 6.6

    def width(t):
        return t.get_window_extent(rend).width / fig.dpi

    ax.imshow(img, extent=(0, side, -side, 0), interpolation="nearest",
              zorder=1)

    def runs(keys):
        out, start = [], 0
        for i in range(1, n + 1):
            if i == n or keys[i] != keys[start]:
                out.append((keys[start], start, i))
                start = i
        return out

    # a hairline of surface colour where one component ends, a wider one
    # where a fold ends
    fold_edge = {e for _, _, e in runs(list(fold))}
    for _, _, e in runs(comp):
        if e < n:
            lw = 0.9 if e in fold_edge else 0.35
            ax.plot([0, e * cell, e * cell], [-e * cell, -e * cell, -side],
                    color=SURFACE, lw=lw, solid_joinstyle="miter", zorder=2)

    # component strip, left of the matrix
    x_strip = -0.10
    for k, (_, s, e) in enumerate(runs(comp)):
        ax.add_patch(Rectangle((x_strip, -e * cell), 0.05, (e - s) * cell,
                               facecolor=STRIP[k % 2], edgecolor="none"))

    # fold brackets and names; the small folds at the end share one label
    x_fold = x_strip - 0.07
    labels, small = [], []
    for f, s, e in runs(list(fold)):
        if e - s >= FOLD_MIN:
            labels.append((s, e, f"{fold_name[f]}  ({e - s})"))
        else:
            small.append((s, e))
    if small:
        k = sum(e - s for s, e in small)
        labels.append((small[0][0], small[-1][1],
                       f"{len(small)} other folds  ({k})"))
    row_h = fs * 1.35 / 72
    want = [-(s + e) / 2 * cell for s, e, _ in labels]
    x_min = x_fold
    for (s, e, txt), y in zip(labels, spread(want, [row_h] * len(labels))):
        y0, y1 = -s * cell - 0.008, -e * cell + 0.008
        ax.plot([x_fold, x_fold], [y0, y1], color=INK, lw=1.1,
                solid_capstyle="butt")
        ax.plot([x_fold, x_fold - 0.10, x_fold - 0.16],
                [(y0 + y1) / 2, y, y], color=MUTED, lw=0.5)
        t = ax.text(x_fold - 0.19, y, txt, ha="right", va="center_baseline",
                    fontsize=fs, color=INK)
        x_min = min(x_min, x_fold - 0.19 - width(t))

    # curated seeds, called out from the diagonal
    x_max = side
    if not a.no_seed_labels and seed:
        idx = sorted(seed)
        want = [-(i + 0.5) * cell for i in idx]
        for i, y in zip(idx, spread(want, [row_h] * len(idx))):
            d = (i + 0.5) * cell
            x = max(-y, d) + 0.22
            ax.plot([d + 0.02, x - 0.08, x - 0.03], [-d, y, y], color=MUTED,
                    lw=0.5, zorder=3)
            ax.plot([d + 0.012], [-d], marker="o", ms=1.8, color=INK, zorder=4)
            t = ax.text(x, y, seed[i], ha="left", va="center_baseline",
                        fontsize=fs, color=INK)
            x_max = max(x_max, x + width(t))

    # ---------------------------------------------------------------- legend
    iu = np.tril_indices(n, -1)
    vals = m[iu]
    aligned = vals[vals > 0]
    share = np.bincount(np.clip(np.searchsorted(lo, aligned, side="right") - 1,
                                0, len(lo) - 1), minlength=len(lo)) / len(aligned)

    L, GAP = 0.30, 0.035
    lx, ly = side * 0.50, -0.30
    if a.centroids:         # clear of the seed labels on a smaller matrix
        lx = max(lx, 4.7)
    if a.title:
        ax.text(lx, ly + 0.32, a.title, fontsize=10, color=INK, ha="left",
                va="bottom")
    ax.text(lx, ly, "Identity (aa %)", fontsize=8.5, color=INK, ha="left",
            va="bottom")
    ax.text(lx, ly - 0.17, f"share of the {len(aligned):,} aligned pairs",
            fontsize=fs, color=MUTED, ha="left", va="bottom")
    bar_h = 0.75
    y_bar = ly - 0.42 - bar_h          # foot of the bars
    y_sw = y_bar - 0.20 - L
    edges = (30, 50, 90)
    x, start, tiers = lx, lx, []
    for k, (b, col) in enumerate(hm.BINS):
        if b in edges:
            tiers.append((start, x))
            x += 0.20
            start = x
        ax.add_patch(Rectangle((x, y_sw), L - GAP, L - GAP, facecolor=col,
                               edgecolor="none"))
        ax.text(x, y_sw + L + 0.02, str(b) if k else "", fontsize=fs,
                color=INK, ha="center", va="bottom")
        h = bar_h * share[k] / share.max()
        ax.add_patch(Rectangle((x, y_bar), L - GAP, h, facecolor=col,
                               edgecolor="none"))
        pct = 100 * share[k]
        ax.text(x + (L - GAP) / 2, y_bar + h + 0.02,
                f"{pct:.0f}" if pct >= 1 else f"{pct:.1f}", fontsize=fs,
                color=INK, ha="center", va="bottom")
        x += L
    tiers.append((start, x))
    x_end = x
    ax.text(x_end - GAP, y_sw + L + 0.02, "100", fontsize=fs, color=INK,
            ha="center", va="bottom")
    ax.plot([lx, x_end - GAP], [y_bar, y_bar], color=MUTED, lw=0.5)
    for (x0, x1), lab in zip(tiers, ["~Fold", "~Family", "~Subfamily",
                                     "~Centroid"]):
        ax.plot([x0, x1 - GAP], [y_sw - 0.05, y_sw - 0.05], color=MUTED, lw=0.8)
        ax.text((x0 + x1 - GAP) / 2, y_sw - 0.20, lab, fontsize=fs, color=INK,
                ha="center", va="baseline")

    n_pairs = len(vals)
    not_same = int(((m == 0) & same_fold)[iu].sum())
    not_other = int(((m == 0) & ~same_fold)[iu].sum())
    y2 = y_sw - 0.85
    ax.text(lx, y2 + L + 0.08,
            f"Not aligned over half of the shorter sequence",
            fontsize=fs, color=INK, ha="left", va="bottom")
    for k, (col, lab) in enumerate(
            [(hm.NO_ALN_SAME_FOLD, f"same fold  ({100 * not_same / n_pairs:.0f}% of pairs)"),
             (hm.NO_ALN_OTHER_FOLD, f"different fold  ({100 * not_other / n_pairs:.0f}%)")]):
        xk = lx + k * 1.75
        ax.add_patch(Rectangle((xk, y2), L - GAP, L - GAP, facecolor=col,
                               edgecolor="none"))
        ax.text(xk + L + 0.03, y2 + (L - GAP) / 2, lab, fontsize=fs, color=INK,
                ha="left", va="center_baseline")
    y3 = y2 - 0.42
    for k, col in enumerate(STRIP):
        ax.add_patch(Rectangle((lx + k * 0.07, y3), 0.07, L - GAP,
                               facecolor=col, edgecolor="none"))
    ax.text(lx + L + 0.03, y3 + (L - GAP) / 2,
            "left strip: SSN components, neighbours alternate", fontsize=fs,
            color=INK, ha="left", va="center_baseline")

    pad = 0.2
    x0, x1 = x_min - pad, max(x_max, x_end + 0.2) + pad
    y0, y1 = -side - pad, pad + (0.25 if a.title else 0)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.set_size_inches(x1 - x0, y1 - y0)
    fig.patch.set_facecolor(SURFACE)
    os.makedirs(out, exist_ok=True)
    for ext in ("svg", "pdf", "png"):
        md = {"Date": None} if ext == "svg" else {}
        fig.savefig(os.path.join(out, f"figure_heatmap_sequences{a.tag}.{ext}"),
                    dpi=DPI, facecolor=SURFACE, metadata=md)
    plt.close(fig)

    # ---------------------------------------------------------------- tables
    write_tsv(os.path.join(out, f"sequence_order{a.tag}.tsv"),
              ["row", "tip", "enzyme_name", "accession", "fold_group", "fold",
               "component", "seed"],
              [[i + 1, t["tip"], t["enzyme_name"], t["accession"],
                t["fold_group"], fold_name[t["fold_group"]], t["component"],
                seed.get(i, "")] for i, t in enumerate(tips)])
    names = [t["tip"] for t in tips]
    write_tsv(os.path.join(out, f"sequence_matrix{a.tag}.tsv"), ["tip"] + names,
              [[names[i]] + ["" if i == j else f"{m[i, j]:.1f}" if m[i, j] else "NA"
                             for j in range(n)] for i in range(n)])

    print(f"set: {a.centroids or SET}. {n} sequences, {n_pairs:,} pairs")
    print(f"  aligned: {len(aligned):,} ({100 * len(aligned) / n_pairs:.1f}%), "
          f"median identity {np.median(aligned):.1f}")
    print(f"  not aligned, same fold: {not_same:,} ({100 * not_same / n_pairs:.1f}%)")
    print(f"  not aligned, different fold: {not_other:,} ({100 * not_other / n_pairs:.1f}%)")
    cm = np.array(comp)
    same_comp = (cm[:, None] == cm[None, :])[iu]
    print(f"  aligned pairs inside one component: "
          f"{int(((vals > 0) & same_comp).sum()):,}; across components: "
          f"{int(((vals > 0) & ~same_comp).sum()):,}, of them "
          f"{int(((vals > 0) & ~same_comp & ~same_fold[iu]).sum()):,} across folds")
    big = collections.Counter(comp).most_common(1)[0]
    inside = (vals > 0) & same_comp & (cm[iu[0]] == big[0])
    print(f"  component {big[0]} ({big[1]} sequences) holds "
          f"{100 * inside.sum() / len(aligned):.0f}% of the aligned pairs")
    for (b, _), s in zip(hm.BINS, share):
        print(f"    {b:>3}+  {100 * s:5.1f}%")


if __name__ == "__main__":
    main()
