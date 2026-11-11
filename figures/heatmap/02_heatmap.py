#!/usr/bin/env python3
"""
Step 2: the triangular matrix of median identity between plasticome groups,
after panel c of the Serratus Coronaviridae figure (coronaviridae-heatmap.webp).

Rows
  One row per SSN component (../ssn, 30 percent identity over 50 percent
  coverage) that has at least --min-seqs distinct sequences. Components named
  in --split are cut into subfamilies first, by average linkage at --sub-id
  percent identity, and each subfamily with at least --min-seqs sequences gets
  a row. By default that is component 0 only, which holds over half the set.
  Rows are grouped by fold, the fold groups of ../tree.

Cells
  Median amino acid identity over every pair of sequences, one from each
  group. The diagonal is the median inside a group. A pair counts as aligned
  when the alignment covers half of the shorter sequence, the SSN rule. Pairs
  that are not aligned rank below every aligned pair, so a cell where half or
  more of the pairs do not align has no median and is drawn in grey: darker
  grey when the two groups share a fold, lighter when they do not.

PLASTICOME_TREE_SET picks the sequences, as in ../tree: all (out/) or
centroids (out_centroids/).

Writes figure_heatmap.{svg,pdf,png}, heatmap_matrix.tsv, heatmap_groups.tsv.
"""

import argparse
import collections
import os
import statistics
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tree"))
from common import (INK, MUTED, OUT as TREE_OUT, SET, SUFFIX, SURFACE,  # noqa: E402
                    WORK_ALL, load_rows, read_tsv, seed_nodes, write_tsv)

WORK = os.path.join(HERE, "work")
OUT = os.path.join(HERE, "out" + SUFFIX)

matplotlib.rcParams["svg.hashsalt"] = "plasticome-v2"
matplotlib.rcParams["svg.fonttype"] = "none"   # labels stay text in Inkscape

COVERAGE = 0.5      # of the shorter sequence, the SSN rule

# Identity bins, lower edge -> colour. Lightness falls with identity across
# the whole scale so the order reads without hue. The hue changes at the
# three thresholds the project already uses: 30 (SSN component), --sub-id
# (subfamily) and 90 (centroid).
BINS = [(20, "#f7ecb0"), (30, "#b4dfc0"), (40, "#7fc7a4"),
        (50, "#74a4e6"), (60, "#4a84d6"), (70, "#2c64b4"), (80, "#1b458a"),
        (90, "#2c1a66")]
FLOOR = BINS[0][0]
NO_ALN_SAME_FOLD = "#c9c7bf"
NO_ALN_OTHER_FOLD = "#ecebe6"
LIGHT_TEXT_FROM = 60    # cells at or above this identity take white numbers


def bin_color(v):
    col = BINS[0][1]
    for lo, c in BINS:
        if v >= lo:
            col = c
    return col


def load_identity(tips):
    """Identity for every pair of `tips`, 0 where the pair is not aligned."""
    ix = {t: i for i, t in enumerate(tips)}
    m = np.zeros((len(tips), len(tips)))
    with open(os.path.join(WORK, "pairs.tsv")) as fh:
        next(fh)
        for line in fh:
            a, b, ident, _, _, cov, *_ = line.split("\t")
            if a in ix and b in ix and float(cov) >= COVERAGE:
                m[ix[a], ix[b]] = m[ix[b], ix[a]] = float(ident)
    return m, ix


def subfamilies(members, m, ix, sub_id):
    """Average linkage clusters of one component at `sub_id` percent."""
    idx = [ix[t] for t in members]
    d = 100.0 - m[np.ix_(idx, idx)]
    np.fill_diagonal(d, 0.0)
    lab = fcluster(linkage(squareform(d, checks=False), "average"),
                   100.0 - sub_id, "distance")
    out = collections.defaultdict(list)
    for t, l in zip(members, lab):
        out[l].append(t)
    return sorted(out.values(), key=lambda g: (-len(g), min(int(t) for t in g)))


def cell(ga, gb, m, ix):
    """Median identity between two groups, or inside one when ga is gb."""
    ia, ib = [ix[t] for t in ga], [ix[t] for t in gb]
    if ga is gb:
        vals = m[np.ix_(ia, ia)][np.triu_indices(len(ia), 1)]
    else:
        vals = m[np.ix_(ia, ib)].ravel()
    aligned = float((vals > 0).mean())
    # Unaligned pairs are 0 and rank lowest. The median is a real identity
    # only when more than half of the pairs align; at exactly half numpy
    # would average a 0 with an identity.
    med = float(np.median(vals)) if aligned > 0.5 else None
    return med, aligned, len(vals)


def name_group(members, tip_row, seed_of):
    """Curated seeds in the group, else its most central named enzyme."""
    seeds = [seed_of[t] for t in members if t in seed_of]
    if seeds:
        return " / ".join(seeds)
    for t in members:     # members arrive most central first
        nm = tip_row[t]["enzyme_name"].split(";")[0].strip()
        if not nm.isdigit():
            return nm
    return tip_row[members[0]]["accession"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-seqs", type=int, default=3)
    ap.add_argument("--split", default="0",
                    help="components to cut into subfamilies, comma separated")
    ap.add_argument("--sub-id", type=float, default=50.0)
    ap.add_argument("--no-numbers", action="store_true")
    a = ap.parse_args()
    split = {f"c{int(c):02d}" for c in a.split.split(",") if c.strip()}

    rows = load_rows()
    # pairs.tsv is keyed by the tips of the full set; a tip of this set is
    # the same protein under what may be another record id
    full = {rows[t["tip"]]["aa_sequence"].strip().upper(): t["tip"]
            for t in read_tsv(os.path.join(WORK_ALL, "tips.tsv"))}
    tip_row = {}
    for t in read_tsv(os.path.join(TREE_OUT, "tree_tips.tsv")):
        tip_row[full[rows[t["tip"]]["aa_sequence"].strip().upper()]] = t
    tips = sorted(tip_row, key=int)
    m, ix = load_identity(tips)

    seed_of = {}
    for name, pid in seed_nodes(rows).items():
        for t, r in tip_row.items():
            if pid in r["records"].split(";"):
                seed_of[t] = name
    fold_name = {r["fold_group"]: r["fold"] for r in
                 read_tsv(os.path.join(TREE_OUT, "fold_table.tsv"))}

    comp = collections.defaultdict(list)
    for t in tips:
        comp[tip_row[t]["component"]].append(t)

    def central_first(g):
        idx = [ix[t] for t in g]
        mean = m[np.ix_(idx, idx)].mean(axis=1)
        return [t for _, t in sorted(zip(-mean, g), key=lambda p: (p[0], int(p[1])))]

    groups, left_out = [], []
    for c in sorted(comp):
        parts = subfamilies(comp[c], m, ix, a.sub_id) if c in split else [comp[c]]
        for g in parts:
            if len(g) < a.min_seqs:
                left_out.append((c, g))
                continue
            g = central_first(g)
            groups.append({"component": c, "fold": tip_row[g[0]]["fold_group"],
                           "members": g, "split": c in split,
                           "name": name_group(g, tip_row, seed_of)})
    groups.sort(key=lambda g: (g["fold"], g["component"], -len(g["members"])))
    n = len(groups)

    os.makedirs(OUT, exist_ok=True)
    cells, table = {}, []
    for i in range(n):
        for j in range(i + 1):
            gi, gj = groups[i], groups[j]
            med, aligned, npairs = cell(gi["members"], gj["members"], m, ix)
            same_fold = gi["fold"] == gj["fold"]
            cells[i, j] = (med, same_fold)
            table.append([i + 1, j + 1, gi["name"], gj["name"],
                          gi["component"], gj["component"],
                          "" if med is None else f"{med:.1f}", npairs,
                          f"{aligned:.3f}", int(same_fold)])
    write_tsv(os.path.join(OUT, "heatmap_matrix.tsv"),
              ["row", "col", "row_name", "col_name", "row_component",
               "col_component", "median_identity", "pairs", "fraction_aligned",
               "same_fold"], table)
    write_tsv(os.path.join(OUT, "heatmap_groups.tsv"),
              ["row", "name", "fold_group", "fold", "component", "subfamily",
               "sequences", "median_identity", "members"],
              [[i + 1, g["name"], g["fold"], fold_name[g["fold"]], g["component"],
                int(g["split"]), len(g["members"]),
                "" if cells[i, i][0] is None else f"{cells[i, i][0]:.1f}",
                ";".join(g["members"])] for i, g in enumerate(groups)])

    # ---------------------------------------------------------------- draw
    S = 0.30                                  # one cell, inches
    GAP = 0.035
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    fs = 6.6
    rend = fig.canvas.get_renderer()

    def extent(t):
        """Width and height of a drawn text, inches."""
        bb = t.get_window_extent(rend)
        return bb.width / fig.dpi, bb.height / fig.dpi

    for (i, j), (med, same_fold) in cells.items():
        col = (bin_color(med) if med is not None else
               NO_ALN_SAME_FOLD if same_fold else NO_ALN_OTHER_FOLD)
        ax.add_patch(Rectangle((j * S + GAP / 2, -(i + 1) * S + GAP / 2),
                               S - GAP, S - GAP, facecolor=col, edgecolor="none",
                               gid=f"cell-{i + 1:02d}-{j + 1:02d}"))
        if med is not None and not a.no_numbers:
            ax.text((j + 0.5) * S, -(i + 0.5) * S, f"{med:.0f}", ha="center",
                    va="center_baseline", fontsize=5.6,
                    color="white" if med >= LIGHT_TEXT_FROM else INK)

    # row and column labels
    x_lab = -0.08
    lab_w = col_h = 0.0
    for i, g in enumerate(groups):
        txt = f"{g['name']}  ({len(g['members'])})"
        t = ax.text(x_lab, -(i + 0.5) * S, txt, ha="right", va="center_baseline",
                    fontsize=fs, color=INK)
        lab_w = max(lab_w, extent(t)[0])
        t = ax.text((i + 0.5) * S + 0.03, -n * S - 0.08, g["name"], ha="right",
                    va="top", rotation=55, rotation_mode="anchor", fontsize=fs,
                    color=INK)
        col_h = max(col_h, extent(t)[1])

    # component brackets (only where a component has more than one row) and
    # fold brackets, left of the row labels
    def span(key):
        out = collections.OrderedDict()
        for i, g in enumerate(groups):
            out.setdefault(g[key], []).append(i)
        return out

    x_comp = x_lab - lab_w - 0.12
    x_fold = x_comp
    for c, idx in span("component").items():
        if len(idx) > 1:
            y0, y1 = -idx[0] * S - 0.05, -(idx[-1] + 1) * S + 0.05
            ax.plot([x_comp, x_comp], [y0, y1], color=MUTED, lw=0.8,
                    solid_capstyle="butt")
            t = ax.text(x_comp - 0.07, (y0 + y1) / 2,
                        f"SSN component {int(c[1:])}\nsubfamilies at {a.sub_id:.0f}%",
                        ha="right", va="center", rotation=90, fontsize=fs,
                        color=MUTED, linespacing=1.25, multialignment="center")
            x_fold = min(x_fold, x_comp - 0.07 - extent(t)[0] - 0.12)
    x_min = x_fold
    for f, idx in span("fold").items():
        y0, y1 = -idx[0] * S - 0.05, -(idx[-1] + 1) * S + 0.05
        ax.plot([x_fold, x_fold], [y0, y1], color=INK, lw=1.1,
                solid_capstyle="butt")
        wrapped = textwrap.fill(fold_name[f], 22)
        t = ax.text(x_fold - 0.07, (y0 + y1) / 2, wrapped, ha="right",
                    va="center", fontsize=fs, color=INK, linespacing=1.15)
        x_min = min(x_min, x_fold - 0.07 - extent(t)[0])

    # legend, in the empty half of the square
    L = 0.26                                   # swatch
    lx = (n * 0.47) * S
    ly = -0.35
    ax.text(lx, ly, "Median identity (aa %)", fontsize=8.5, color=INK,
            ha="left", va="bottom")
    y_sw = ly - 0.30 - L

    def swatch(x, color, label=None):
        ax.add_patch(Rectangle((x, y_sw), L - GAP, L - GAP, facecolor=color,
                               edgecolor="none"))
        if label is not None:
            ax.text(x, y_sw + L + 0.02, label, fontsize=fs, color=INK,
                    ha="center", va="bottom")

    def tier(x0, x1, label):
        ax.plot([x0, x1 - GAP], [y_sw - 0.05, y_sw - 0.05], color=MUTED, lw=0.8)
        ax.text((x0 + x1 - GAP) / 2, y_sw - 0.20, label, fontsize=fs,
                color=INK, ha="center", va="baseline")

    edges = {"family": 30, "subfamily": a.sub_id, "centroid": 90}
    x = lx
    tiers, start = [], x
    prev = None
    for lo, col in BINS:
        t = ("below" if lo < edges["family"] else "family" if lo < edges["subfamily"]
             else "subfamily" if lo < edges["centroid"] else "centroid")
        if prev is not None and t != prev:
            tiers.append((prev, start, x))
            x += 0.20
            start = x
        swatch(x, col, str(lo))
        x += L
        prev = t
    tiers.append((prev, start, x))
    x_end = x
    ax.text(x_end - GAP, y_sw + L + 0.02, "100", fontsize=fs, color=INK,
            ha="center", va="bottom")
    tier_label = {"below": "~Fold", "family": "~Family",
                  "subfamily": "~Subfamily", "centroid": "~Centroid"}
    for t, x0, x1 in tiers:
        tier(x0, x1, tier_label[t])

    y2 = y_sw - 0.85
    ax.text(lx, y2 + L + 0.08, "No median: half or more of the pairs do not align",
            fontsize=fs, color=INK, ha="left", va="bottom")
    for k, (col, lab) in enumerate([(NO_ALN_SAME_FOLD, "same fold"),
                                    (NO_ALN_OTHER_FOLD, "different fold")]):
        xk = lx + k * 1.15
        ax.add_patch(Rectangle((xk, y2), L - GAP, L - GAP, facecolor=col,
                               edgecolor="none"))
        ax.text(xk + L + 0.03, y2 + (L - GAP) / 2, lab, fontsize=fs, color=INK,
                ha="left", va="center_baseline")

    pad = 0.2
    x0, x1 = x_min - pad, max(n * S, x_end + 0.2) + pad
    y0, y1 = -n * S - 0.08 - col_h - pad, pad
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.set_size_inches(x1 - x0, y1 - y0)
    fig.patch.set_facecolor(SURFACE)
    for ext in ("svg", "pdf", "png"):
        md = {"Date": None} if ext == "svg" else {}
        fig.savefig(os.path.join(OUT, f"figure_heatmap.{ext}"), dpi=300,
                    facecolor=SURFACE, metadata=md)
    plt.close(fig)

    # ---------------------------------------------------------------- report
    shown = sum(len(g["members"]) for g in groups)
    print(f"set: {SET}. {n} rows, {shown} of {len(tips)} sequences "
          f"({len(tips) - shown} in groups under {a.min_seqs} sequences)")
    by_c = collections.Counter()
    for c, g in left_out:
        by_c[c] += len(g)
    print("  left out: " + ", ".join(f"{c} {k}" for c, k in sorted(by_c.items())
                                     if c in split or k > 0)[:400])
    diag = [cells[i, i][0] for i in range(n) if cells[i, i][0] is not None]
    off = [v[0] for (i, j), v in cells.items() if i != j]
    num = [v for v in off if v is not None]
    print(f"  diagonal: median inside a group {min(diag):.0f} to {max(diag):.0f}, "
          f"middle {statistics.median(diag):.0f}")
    print(f"  off diagonal: {len(num)} of {len(off)} cells have a median"
          + (f", {min(num):.0f} to {max(num):.0f}" if num else ""))
    cross = [(groups[i]["name"], groups[j]["name"], v[0]) for (i, j), v in
             cells.items() if i != j and v[0] is not None
             and groups[i]["component"] != groups[j]["component"]]
    print(f"  cells with a median between two different components: {len(cross)}")
    for r in cross:
        print(f"    {r[0]} x {r[1]}: {r[2]:.1f}")


if __name__ == "__main__":
    main()
