#!/usr/bin/env python3
"""
Step 7: draw the forest.

The plasticome has no single tree, so the figure is a forest: one tree per
fold group. The largest fold group is drawn as a fan, the rest as small
rectangular trees beside it, and fold groups that are one sequence are listed.
Every tree shares one scale, so a branch length means the same everywhere.

Encoding (tips match panel a-i, ../ssn):
  tip colour    substrate
  tip shape     diamond bacterial, hexagon eukaryotic, circle metagenomic
  ringed tip    one of the curated seed enzymes
  dashed black  backbone between components (structure based, step 5), drawn
                at a fixed length and not to scale, with its UFBoot printed at
                the node. Unsupported backbone
                branches were already collapsed in step 6, so a fork that is
                drawn is a fork that is claimed
  dark grey     branch inside a component with ensemble support >= 70, or a
                terminal branch
  light grey    branch inside a component with ensemble support below 70
  outer arc     one SSN component (fan only), numbered when it has 4 or more
                tips

Writes out/figure_tree.{svg,pdf,png}.
"""

import collections
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "plasticome-v2-tree"
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge

from common import (INK, MUTED, NODE_STROKE, OUT, SHAPE_LABEL, SUBSTRATE_COLOR,
                    SUBSTRATE_ORDER, SURFACE, WORK, cdir, load_rows, read_tsv,
                    seed_nodes)
from build_ssn import BLOB, regroup_svg

sys.setrecursionlimit(20000)

SUPPORT_MIN = 70
STYLE = {                       # colour, linewidth, zorder
    "backbone": (INK, 1.0, 3),
    "solid": ("#555555", 0.6, 2),
    "weak": ("#c2c2be", 0.6, 1),
}
FAN_RADIUS = 4.4                # inches from the centre to the deepest tip of the fan
FAN_INNER = 0.75                # radius of the root. At zero, the stems of 36
                                # components meet in a point and their dashes
                                # merge into rings that are not in the tree
BIG_COMPONENT = 100             # tips. A component this large gets tighter
BIG_WEIGHT = 0.6                # tip spacing in the fan, so that one component
                                # does not take three quarters of the circle
ROW = 0.108                     # inches per tip in the rectangular trees
HEAD = 0.55                     # inches reserved for a tree heading
GAP = 0.28                      # inches between stacked trees
LABEL_W = 1.75                  # inches reserved for tip labels
COL_GAP = 0.35
LABEL_MIN_TIPS = 4              # fan components this size or larger get a number
BACKBONE_STEP = 0.42            # drawn length of a backbone branch, in inches


def walk(tree, multi, step, weight=None):
    """Flatten a tree into nodes with depth, tip-order position and the style
    of the branch above each node. `step` is the drawn length of a backbone
    branch in tree units. `weight` maps a tip to the share of a row it takes;
    y, lo and hi are in rows, y at the middle of the tip."""
    nodes, counter = [], [0.0]
    weight = weight or {}

    def rec(cl, parent, up, in_comp, comp):
        name = cl.name or "X"
        if name.startswith("S_"):
            comp = name[2:]
        node = {"clade": cl, "parent": parent, "comp": comp,
                "name": name, "tip": cl.is_terminal(), "boot": None}
        if not in_comp and not name.startswith("C_"):
            node["edge"] = "backbone"
            if name.startswith("B_"):
                node["boot"] = int(name.split("_")[1])
        elif name.startswith("C_"):
            node["edge"] = "solid" if int(name.split("_")[2]) >= SUPPORT_MIN else "weak"
        else:
            node["edge"] = "solid"
        # Backbone branches are drawn at a fixed length. Their estimates run
        # from 3 to 11 substitutions per site, which is saturation, not a
        # measurement, and to scale they would dwarf every component.
        if parent is None:
            depth = 0.0
        elif node["edge"] == "backbone":
            depth = up + step
        else:
            depth = up + (cl.branch_length or 0.0)
        node["depth"] = depth
        below = in_comp or name.startswith("S_")
        node["fork"] = "solid" if below else "backbone"
        nodes.append(node)
        if cl.is_terminal():
            w = weight.get(name, 1.0)
            node["lo"], node["hi"] = counter[0], counter[0] + w
            node["y"] = counter[0] + w / 2
            node["ntips"] = 1
            counter[0] += w
        else:
            kids = [rec(ch, node, depth, below, comp) for ch in cl.clades]
            node["kids"] = kids
            node["lo"] = min(k["lo"] for k in kids)
            node["hi"] = max(k["hi"] for k in kids)
            node["ntips"] = sum(k["ntips"] for k in kids)
            node["y"] = (kids[0]["y"] + kids[-1]["y"]) / 2
        return node

    rec(tree.root, None, 0.0, not multi, None)
    return nodes


class Canvas:
    """One axes in inches. Collects line segments by style and artists by tree,
    so each tree becomes one movable group in the SVG."""

    def __init__(self):
        self.fig, self.ax = plt.subplots()
        self.segs = collections.defaultdict(list)
        self.gids = collections.defaultdict(list)
        self.box = [0.0, 0.0, 0.0, 0.0]
        self.n = 0

    def grow(self, x, y):
        self.box = [min(self.box[0], x), min(self.box[1], y),
                    max(self.box[2], x), max(self.box[3], y)]

    def tag(self, artist, key):
        self.n += 1
        gid = f"plt-{key}-{self.n}"
        artist.set_gid(gid)
        self.gids[key].append(gid)
        return artist

    def line(self, key, style, pts):
        self.segs[(key, style)].append(pts)
        for x, y in pts:
            self.grow(x, y)

    def flush(self):
        for (key, style), segs in sorted(self.segs.items(), key=str):
            col, lw, z = STYLE[style] if style in STYLE else style
            dashed = style == "backbone"
            self.tag(self.ax.add_collection(LineCollection(
                segs, colors=col, linewidths=lw, zorder=z,
                linestyles=[(0, (2.2, 1.6))] if dashed else "solid",
                capstyle="butt" if dashed else "round", joinstyle="round")), key)

    def text(self, key, x, y, s, **kw):
        # matplotlib does not report text extent before drawing, so reserve
        # room from an estimate, or a long label runs off the canvas
        w = max(len(line) for line in s.split("\n")) * kw.get("fontsize", 10) * 0.0082
        ha = kw.get("ha", "left")
        self.grow(x - (w if ha == "right" else w / 2 if ha == "center" else 0), y)
        self.grow(x + (w if ha == "left" else w / 2 if ha == "center" else 0), y)
        return self.tag(self.ax.text(x, y, s, **kw), key)

    def markers(self, key, pts, z=5):
        """pts: (x, y, shape, colour, is_seed)."""
        by = collections.defaultdict(list)
        for x, y, mk, col, seed in pts:
            by[(mk, col)].append((x, y))
            self.grow(x, y)
        for (mk, col), xy in sorted(by.items()):
            a = np.array(xy)
            self.tag(self.ax.scatter(a[:, 0], a[:, 1], marker=mk, s=15, c=col,
                                     edgecolors=NODE_STROKE, linewidths=0.35,
                                     zorder=z), key)
        ring = collections.defaultdict(list)
        for x, y, mk, col, seed in pts:
            if seed:
                ring[mk].append((x, y))
        for mk, xy in sorted(ring.items()):
            a = np.array(xy)
            self.tag(self.ax.scatter(a[:, 0], a[:, 1], marker=mk, s=38,
                                     facecolors="none", edgecolors=INK,
                                     linewidths=0.9, zorder=z + 1), key)


def heading(cv, key, x, y, title, sub):
    cv.text(key, x, y, title, ha="left", va="top", fontsize=10.5, color=INK, zorder=8)
    cv.text(key, x, y - 0.21, sub, ha="left", va="top", fontsize=7.5, color=MUTED,
            zorder=8)


def draw_fan(cv, key, tree, tipinfo, seeds, cx, cy):
    # The backbone step is fixed in inches, the rest of the fan is to scale,
    # and the two together must fill FAN_RADIUS. Solve for the scale.
    size = collections.Counter(t["component"] for t in tipinfo.values())
    weight = {name: BIG_WEIGHT for name, t in tipinfo.items()
              if size[t["component"]] >= BIG_COMPONENT}
    # single sequence components sit at the end of a bare stem, close to the
    # centre where rows are narrow; two rows each keeps their markers apart
    weight.update({name: 2.0 for name, t in tipinfo.items()
                   if size[t["component"]] == 1})
    reach = FAN_RADIUS - FAN_INNER
    scale = reach / max(d["depth"] for d in walk(tree, True, 0.0))
    for _ in range(6):
        nodes = walk(tree, True, BACKBONE_STEP / scale, weight)
        scale = reach / max(d["depth"] for d in nodes)
    nodes = walk(tree, True, BACKBONE_STEP / scale, weight)
    n = max(d["hi"] for d in nodes)             # rows, after weighting
    gap = max(4.0, n * 0.02)                    # open wedge at the top, in rows

    def ang(y):
        return math.pi / 2 - 2 * math.pi * (y + gap / 2) / (n + gap)

    def xy(r, a):
        return (cx + r * math.cos(a), cy + r * math.sin(a))

    def arc(r, a0, a1):
        k = max(2, int(abs(a1 - a0) / math.radians(1.0)) + 2)
        return [xy(r, a) for a in np.linspace(a0, a1, k)]

    for d in nodes:
        d["r"], d["a"] = FAN_INNER + d["depth"] * scale, ang(d["y"])
    for d in nodes:
        if d["parent"] is not None:
            cv.line(key, d["edge"], [xy(d["parent"]["r"], d["a"]), xy(d["r"], d["a"])])
        if not d["tip"]:
            cv.line(key, d["fork"], arc(d["r"], d["kids"][0]["a"], d["kids"][-1]["a"]))
        if d["boot"] is not None:
            x, y = xy(d["r"], d["a"])
            cv.text(key, x, y, str(d["boot"]), fontsize=5.5, color=INK, zorder=9,
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.12", fc=SURFACE, ec="none"))

    # component shading and the outer arc that names it
    r_ring = FAN_RADIUS + 0.16
    pad = 0.2 * 2 * math.pi / (n + gap)
    stems = [d for d in nodes if d["name"].startswith("S_")]
    for i, d in enumerate(sorted(stems, key=lambda d: d["lo"])):
        a0, a1 = ang(d["lo"]) - pad, ang(d["hi"]) + pad
        cv.tag(cv.ax.add_patch(Wedge((cx, cy), r_ring - 0.05, math.degrees(a1),
                                     math.degrees(a0), width=r_ring - 0.05 - d["r"],
                                     facecolor=BLOB if i % 2 else "#f0f0ee",
                                     edgecolor="none", zorder=0)), key)
        cv.line(key, (MUTED, 1.4, 4), arc(r_ring, a0, a1))
        if d["ntips"] >= LABEL_MIN_TIPS:
            am = (a0 + a1) / 2
            x, y = xy(r_ring + 0.09, am)
            cv.text(key, x, y, d["comp"], fontsize=6.5, color=MUTED, zorder=8,
                    ha="left" if math.cos(am) > 0.3 else
                    "right" if math.cos(am) < -0.3 else "center",
                    va="bottom" if math.sin(am) > 0.3 else
                    "top" if math.sin(am) < -0.3 else "center")

    pts = []
    for d in nodes:
        if d["tip"]:
            t = tipinfo[d["name"]]
            x, y = xy(d["r"], d["a"])
            pts.append((x, y, t["shape"], SUBSTRATE_COLOR[t["substrate_group"]],
                        d["name"] in seeds))
    cv.markers(key, pts)

    # seed names in two columns outside the ring, spread so they do not collide
    r_lab = r_ring + 0.55
    called = [(d, seeds[d["name"]]) for d in nodes if d["tip"] and d["name"] in seeds]
    for side in (1, -1):
        col = sorted((c for c in called if (math.cos(c[0]["a"]) >= 0) == (side > 0)),
                     key=lambda c: -math.sin(c[0]["a"]))
        ys = [cy + r_lab * math.sin(d["a"]) for d, _ in col]
        for i in range(1, len(ys)):
            ys[i] = min(ys[i], ys[i - 1] - 0.17)
        for (d, name), ly in zip(col, ys):
            dy = max(-r_lab, min(r_lab, ly - cy))
            lx = cx + side * max(math.sqrt(max(r_lab ** 2 - dy ** 2, 0.0)), 0.6)
            # out along the tip's own radius first, so the leader does not
            # cut across neighbouring clades
            cv.line(key, (MUTED, 0.5, 6), [xy(d["r"] + 0.05, d["a"]),
                                           xy(r_ring + 0.12, d["a"]), (lx, ly)])
            cv.text(key, lx + side * 0.04, ly, name, fontsize=7.5, color=INK,
                    zorder=9, ha="left" if side > 0 else "right", va="center")
    return scale


def draw_rect(cv, key, tree, multi, tipinfo, seeds, x0, y_top, scale):
    nodes = walk(tree, multi, BACKBONE_STEP / scale)
    X = lambda d: x0 + d["depth"] * scale
    Y = lambda y: y_top - y * ROW
    pts = []
    for d in nodes:
        if d["parent"] is not None:
            cv.line(key, d["edge"], [(X(d["parent"]), Y(d["y"])), (X(d), Y(d["y"]))])
        if not d["tip"]:
            cv.line(key, d["fork"], [(X(d), Y(d["kids"][0]["y"])),
                                     (X(d), Y(d["kids"][-1]["y"]))])
        if d["boot"] is not None:
            cv.text(key, X(d) - 0.03, Y(d["y"]) + 0.02, str(d["boot"]), fontsize=5.5,
                    color=INK, zorder=9, ha="right", va="bottom")
        if d["tip"]:
            t = tipinfo[d["name"]]
            pts.append((X(d), Y(d["y"]), t["shape"],
                        SUBSTRATE_COLOR[t["substrate_group"]], d["name"] in seeds))
            name = t["enzyme_name"]
            cv.text(key, X(d) + 0.07, Y(d["y"]), name[:27] + ("..." if len(name) > 27 else ""),
                    fontsize=6, color=INK, ha="left", va="center", zorder=8)
    cv.markers(key, pts)
    n = sum(1 for d in nodes if d["tip"])
    return max(d["depth"] for d in nodes) * scale, n * ROW


def main():
    rows = load_rows()
    tips = read_tsv(os.path.join(OUT, "tree_tips.tsv"))
    tipinfo = {t["tip"]: t for t in tips}
    record_tip = {r: t["tip"] for t in tips for r in t["records"].split(";")}
    seeds = {record_tip[pid]: name for name, pid in seed_nodes(rows).items()
             if pid in record_tip}

    groups = collections.defaultdict(list)
    for t in tips:
        groups[t["fold_group"]].append(t)
    n_comp = {g: len({t["component"] for t in ts}) for g, ts in groups.items()}

    def title(g):
        doms = collections.Counter(
            (rows[r]["domain"] or "NA").split(";")[0].strip()
            for t in groups[g] for r in t["records"].split(";"))
        name = doms.most_common(1)[0][0]
        name = " ".join(name.replace("[topology]", "").split())
        key = groups[g][0]["fold_key"]
        cath = f"CATH {key}, " if not key.startswith("domain:") else ""
        nc, nt = n_comp[g], len(groups[g])
        return (name[0].upper() + name[1:],
                f"{cath}{nc} component{'s' if nc != 1 else ''}, "
                f"{nt} sequence{'s' if nt != 1 else ''}")

    trees = {g: Phylo.read(os.path.join(OUT, "trees", f"{g}.nwk"), "newick")
             for g in groups}
    order = sorted(groups, key=lambda g: (-len(groups[g]), g))
    fan, rest = order[0], order[1:]
    small = [g for g in rest if len(groups[g]) > 1]
    single = [g for g in rest if len(groups[g]) == 1]

    cv = Canvas()
    block = FAN_RADIUS + 1.45            # half width of the fan block
    cx, cy = block, -(block + HEAD)
    heading(cv, fan, 0.0, 0.0, *title(fan))
    scale = draw_fan(cv, fan, trees[fan], tipinfo, seeds, cx, cy)

    # --- small trees, stacked into columns as tall as the fan block ---------
    col_h = 2 * block + HEAD
    x, y = 2 * block + COL_GAP + 0.5, 0.0
    col_w = 0.0
    for g in small:
        need = HEAD + len(groups[g]) * ROW
        if y < 0 and -y + need > col_h:
            x += col_w + LABEL_W + COL_GAP
            y, col_w = 0.0, 0.0
        heading(cv, g, x, y, *title(g))
        w, h = draw_rect(cv, g, trees[g], n_comp[g] > 1, tipinfo, seeds,
                         x + 0.05, y - HEAD, scale)
        col_w = max(col_w, w, 1.6)
        y -= HEAD + h + GAP

    # --- fold groups that are a single sequence ------------------------------
    if single:
        need = HEAD + len(single) * ROW * 1.6
        if -y + need > col_h:
            x += col_w + LABEL_W + COL_GAP
            y, col_w = 0.0, 0.0
        key = "singles"
        heading(cv, key, x, y, "Single sequence fold groups",
                f"{len(single)} sequences with no homolog in the set")
        pts = []
        for i, g in enumerate(single):
            t = groups[g][0]
            yy = y - HEAD - (i + 0.5) * ROW * 1.6
            pts.append((x + 0.1, yy, t["shape"], SUBSTRATE_COLOR[t["substrate_group"]],
                        t["tip"] in seeds))
            note = ", homology unconfirmed" if t["homology"] == "unconfirmed" else ""
            cv.text(key, x + 0.2, yy, f"{t['enzyme_name']}  ({title(g)[0]}{note})",
                    fontsize=6, color=INK, ha="left", va="center", zorder=8)
        cv.markers(key, pts)
        col_w = max(col_w, 1.6)
    cv.grow(x + col_w + LABEL_W, 0.0)
    cv.flush()

    # --- legends in the free corners of the fan block ------------------------
    ax = cv.ax
    present = [s for s in SUBSTRATE_ORDER if any(t["substrate_group"] == s for t in tips)]
    h_sub = [Line2D([], [], marker="s", linestyle="none", markersize=6.5,
                    markerfacecolor=SUBSTRATE_COLOR[s], markeredgecolor=NODE_STROKE,
                    markeredgewidth=0.4, label=s) for s in present]
    h_shape = [Line2D([], [], marker=mk, linestyle="none", markersize=6.5,
                      markerfacecolor="none", markeredgecolor=INK,
                      markeredgewidth=0.8, label=lab) for lab, mk in SHAPE_LABEL]
    h_shape.append(Line2D([], [], marker="o", linestyle="none", markersize=8,
                          markerfacecolor="none", markeredgecolor=INK,
                          markeredgewidth=0.9, label="curated seed"))
    h_br = [Line2D([], [], color=STYLE["backbone"][0], linewidth=STYLE["backbone"][1],
                   linestyle=(0, (2.2, 1.6)),
                   label="backbone between components,\nnot to scale (number: UFBoot)"),
            Line2D([], [], color=STYLE["solid"][0], linewidth=1.0,
                   label=f"within component, support >= {SUPPORT_MIN}"),
            Line2D([], [], color=STYLE["weak"][0], linewidth=1.0,
                   label=f"within component, support < {SUPPORT_MIN}"),
            Line2D([], [], color=MUTED, linewidth=1.4, label="one SSN component")]
    y_bot = -(2 * block + HEAD)

    def legend(handles, ttl, x, y, loc, ncol=1):
        lg = ax.legend(handles=handles, title=ttl, loc=loc, bbox_to_anchor=(x, y),
                       bbox_transform=ax.transData, frameon=False, fontsize=7.5,
                       ncol=ncol, handletextpad=0.5, columnspacing=1.0,
                       labelspacing=0.5, borderaxespad=0.0, alignment="left")
        lg.get_title().set_fontsize(7.5)
        lg.get_title().set_color(MUTED)
        ax.add_artist(lg)
        cv.tag(lg, "legend")

    legend(h_sub, "substrate", 0.0, y_bot, "lower left", ncol=2)
    legend(h_shape, "origin", 2 * block, y_bot, "lower right")
    legend(h_br, "branches", 2 * block, -HEAD + 0.1, "upper right")

    # scale bar, valid for every solid branch in the figure
    bar = 0.5
    bx, by = 0.02, y_bot + 1.75
    cv.tag(ax.add_line(Line2D([bx, bx + bar * scale], [by, by], color=INK,
                              linewidth=1.0, zorder=8)), "legend")
    cv.text("legend", bx, by + 0.06, f"{bar} substitutions per site", fontsize=7.5,
            color=MUTED, ha="left", va="bottom")
    cv.grow(0.0, y_bot)

    pad = 0.2
    x0, y0, x1, y1 = cv.box
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    cv.fig.set_size_inches(x1 - x0 + 2 * pad, y1 - y0 + 2 * pad)
    cv.fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    cv.fig.patch.set_facecolor(SURFACE)

    for ext in ("svg", "pdf", "png"):
        md = {"Date": None} if ext == "svg" else {}
        cv.fig.savefig(os.path.join(OUT, f"figure_tree.{ext}"), dpi=300,
                       facecolor=SURFACE, metadata=md)
    plt.close(cv.fig)

    labels = {g: f"{g} {title(g)[0]} ({title(g)[1]})" for g in groups}
    labels.update(singles="single sequence fold groups", legend="legends and scale bar")
    moved = regroup_svg(os.path.join(OUT, "figure_tree.svg"),
                        [(f"tree-{k}", labels[k], v) for k, v in sorted(cv.gids.items())])
    print(f"fan: {title(fan)[0]}, {len(groups[fan])} tips; {len(small)} small trees; "
          f"{len(single)} single sequence fold groups")
    print(f"seeds called out: {len(seeds)}; scale {scale:.2f} in per substitution per site")
    print(f"canvas {x1 - x0 + 2 * pad:.1f} x {y1 - y0 + 2 * pad:.1f} in; "
          f"svg: {len(cv.gids)} groups, {moved} elements grouped")
    if os.path.exists(os.path.join(cdir(0), "PROVISIONAL")):
        print("WARNING: component 0 is using a provisional tree")


if __name__ == "__main__":
    main()
