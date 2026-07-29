#!/usr/bin/env python3
"""Plot the cluster size distribution for the v260701-union reference-seeded clustering.
Reads results/cluster_membership.tsv, writes results/cluster_size_distribution.png."""
import csv, collections, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def p(*a): return os.path.join(ROOT, *a)

# one row per cluster -> its size
size = {}
for r in csv.DictReader(open(p("results", "cluster_membership.tsv")), delimiter="\t"):
    size[r["cluster_id"]] = int(r["size"])

dist = collections.Counter(size.values())          # size -> #clusters
sizes = sorted(dist)                                # 1,2,3,4,5,6,21
counts = [dist[s] for s in sizes]
records = [s * dist[s] for s in sizes]
total_clusters = len(size)

INK   = "#1a2027"    # primary text
MUTED = "#5b6670"    # secondary text
GRID  = "#e6e9ec"    # recessive grid
BAR   = "#2f6f9f"    # single-series hue (accessible on white, CVD-safe alone)

fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=150)
fig.patch.set_facecolor("white"); ax.set_facecolor("white")

x = range(len(sizes))
bars = ax.bar(x, counts, width=0.66, color=BAR, zorder=3)

ax.set_yscale("log")
ax.set_ylim(0.8, max(counts) * 2.2)

# direct value labels above each bar: #clusters (bold) + records held (muted)
for xi, s, c, rec in zip(x, sizes, counts, records):
    ax.annotate(f"{c}", (xi, c), textcoords="offset points", xytext=(0, 16),
                ha="center", va="bottom", fontsize=10.5, fontweight="bold",
                color=INK, zorder=4)
    ax.annotate(f"{rec} seqs", (xi, c), textcoords="offset points", xytext=(0, 4),
                ha="center", va="bottom", fontsize=7.5, color=MUTED, zorder=4)

ax.set_xticks(list(x))
ax.set_xticklabels([str(s) for s in sizes])
ax.set_xlabel("cluster size (sequences per cluster)", fontsize=11, color=INK)
ax.set_ylabel("number of clusters  (log scale)", fontsize=11, color=INK)
ax.set_title("Cluster size distribution — plasticome union v1.v260701 @ 90% id",
             fontsize=12.5, fontweight="bold", color=INK, pad=32)
ax.text(0.0, 1.045,
        f"611 sequences → {total_clusters} clusters (13 reference + 399 de novo)  ·  "
        f"{dist[1]} singletons ({100*dist[1]/total_clusters:.0f}%)  ·  largest = {max(sizes)}",
        transform=ax.transAxes, fontsize=9.5, color=MUTED)

# recessive grid & spines
ax.yaxis.grid(True, which="major", color=GRID, lw=1, zorder=0)
ax.yaxis.grid(True, which="minor", color=GRID, lw=0.5, alpha=0.6, zorder=0)
ax.yaxis.set_minor_locator(LogLocator(base=10, subs=tuple(range(2, 10)), numticks=12))
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"): ax.spines[sp].set_color(MUTED)
ax.tick_params(colors=MUTED)

fig.tight_layout()
out = p("results", "cluster_size_distribution.png")
fig.savefig(out, facecolor="white", bbox_inches="tight")
print("wrote", out)
print("size:count ->", dict(zip(sizes, counts)))
