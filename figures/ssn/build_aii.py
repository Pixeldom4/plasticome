#!/usr/bin/env python3
"""
Build panel a-ii: enzyme count per catalytic family, stacked by substrate.

Shares the palette, substrate mapping and family mapping with build_ssn.py so
the two panels cannot drift apart.

Usage:  python build_aii.py
"""

import collections
import os

import matplotlib
matplotlib.use("Agg")
# Fixed salt: matplotlib otherwise seeds SVG element ids (clip paths,
# marker defs) from a random uuid per figure, so two identical renders
# produce different files and the SVG cannot be diffed or version tracked.
matplotlib.rcParams["svg.hashsalt"] = "plasticome-v2"
# Emit real <text> elements in the SVG instead of glyph outlines. The default,
# svg.fonttype="path", traces every label as vector paths, so in Inkscape the
# wording cannot be edited or retyped. The font is referenced by name rather
# than embedded, so the editing machine needs DejaVu Sans (matplotlib ships it)
# or Inkscape substitutes and the metrics shift slightly.
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import build_ssn as S

OUT = S.OUT

# The published panel bins the long tail of families into "Other". Order is
# fixed (largest first, Other last) so colour never follows rank.
FAMILY_BINS = [
    ("A/B hydrolase", ["A/B hydrolase"]),
    ("Amidase", ["Amidase", "L-amino peptidase\nD-ALA esterase"]),
    ("Peptidase", ["Peptidase", "Trypsin-like"]),
    ("Oxygenase", ["Oxygenase-like", "Cytochrome C-like"]),
]


def bin_of(fam):
    for label, members in FAMILY_BINS:
        if fam in members:
            return label
    return "Other"


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = S.load()

    counts = collections.defaultdict(collections.Counter)
    for r in rows.values():
        counts[bin_of(S.family_of(r["domain"]))][S.substrate_of(r["substrate"])] += 1

    order = [lab for lab, _ in FAMILY_BINS] + ["Other"]
    order = [o for o in order if counts[o]]
    subs = [s for s in S.SUBSTRATE_ORDER
            if any(counts[o][s] for o in order)]

    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    fig.patch.set_facecolor(S.SURFACE)
    ax.set_facecolor(S.SURFACE)

    x = range(len(order))
    bottoms = [0.0] * len(order)
    for s in subs:
        vals = [counts[o][s] for o in order]
        ax.bar(x, vals, bottom=bottoms, width=0.66,
               color=S.SUBSTRATE_COLOR[s], edgecolor=S.SURFACE, linewidth=1.4,
               zorder=2)
        bottoms = [b + v for b, v in zip(bottoms, vals)]

    # total above each bar: the one number a reader always wants
    for xi, tot in zip(x, bottoms):
        ax.text(xi, tot + max(bottoms) * 0.015, f"{int(tot)}", ha="center",
                va="bottom", fontsize=9.5, color=S.INK, zorder=3)

    ax.set_ylabel("count", fontsize=10, color=S.INK)
    ax.set_xticks(list(x))
    ax.set_xticklabels(order, rotation=35, ha="right", fontsize=9.5,
                       color=S.INK)
    ax.set_ylim(0, max(bottoms) * 1.10)
    ax.tick_params(axis="y", labelsize=9, colors=S.MUTED, length=0)
    ax.tick_params(axis="x", length=0)
    ax.yaxis.grid(True, color="#e4e4e1", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    handles = [Line2D([], [], marker="s", linestyle="none", markersize=8,
                      markerfacecolor=S.SUBSTRATE_COLOR[s],
                      markeredgecolor=S.NODE_STROKE, markeredgewidth=0.4,
                      label=s) for s in subs]
    leg = ax.legend(handles=handles, loc="upper right", ncol=2, frameon=False,
                    fontsize=9, handletextpad=0.5, columnspacing=1.1,
                    labelspacing=0.5, title="substrate")
    leg.get_title().set_fontsize(9)
    leg.get_title().set_color(S.MUTED)

    fig.tight_layout()
    for ext in ("svg", "pdf", "png"):
        md = {"Date": None} if ext == "svg" else {}
        fig.savefig(os.path.join(OUT, f"figure_aii.{ext}"), dpi=400,
                    facecolor=S.SURFACE, metadata=md)
    plt.close(fig)

    # the table behind the bars
    import csv
    with open(os.path.join(OUT, "figure_aii_counts.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family"] + subs + ["total"])
        for o in order:
            w.writerow([o] + [counts[o][s] for s in subs] +
                       [sum(counts[o].values())])
        w.writerow(["TOTAL"] + [sum(counts[o][s] for o in order) for s in subs] +
                   [sum(sum(counts[o].values()) for o in order)])

    print("panel a-ii:")
    for o in order:
        print(f"   {o:<16} {sum(counts[o].values()):>4}")
    print(f"   {'TOTAL':<16} {sum(sum(counts[o].values()) for o in order):>4}")
    print(f"   wrote out/figure_aii.{{svg,pdf,png}}, figure_aii_counts.tsv")


if __name__ == "__main__":
    main()
