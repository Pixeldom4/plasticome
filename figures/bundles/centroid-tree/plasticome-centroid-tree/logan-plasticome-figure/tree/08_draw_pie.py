#!/usr/bin/env python3
"""
Step 8: the forest as one pie, cut into slices by fold.

Same trees as step 7, different layout. Every fold group is one slice of a
single circle, the slices are separated by a gap, and the background of a
slice is coloured by its fold. Nothing joins two slices: each tree starts on
the inner circle on its own, because trees in different fold groups have no
common ancestor in this analysis.

Encoding
  slice           one fold group (one tree)
  slice colour    the fold. Named outside the slice as well, so colour is
                  never the only cue. The name is all the label carries; the
                  CATH accession and the component and sequence counts are in
                  the fold table (also written to out/fold_table.tsv). Fold groups that are a single sequence
                  share one neutral grey
  two tints       inside a slice with several components, neighbouring
                  components alternate between two tints of the fold colour
  tip, branches   as in step 7

Tip names are off by default; only the curated seeds are called out, in every
slice. --tip-names puts the names back along the rim (all slices except the
largest, where 377 names do not fit).

Branch lengths share one scale across all slices. Dashed backbone branches
are not to scale, as in step 7.

The fold table (name, CATH accession, components, sequences) is a figure of
its own, so it can be placed independently of the pie.

Writes out/figure_pie.{svg,pdf,png}, out/figure_pie_table.{svg,pdf,png} and
out/fold_table.tsv. Step 7 and out/figure_tree.* are left alone.

Usage: python 08_draw_pie.py [--tip-names]
"""

import argparse
import collections
import colorsys
import importlib
import math
import os

import numpy as np
from Bio import Phylo
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Wedge

from common import (INK, MUTED, NODE_STROKE, OUT, SHAPE_LABEL, SUBSTRATE_COLOR,
                    SUBSTRATE_ORDER, SURFACE, load_rows, read_tsv, seed_nodes,
                    write_tsv)
from build_ssn import regroup_svg

d7 = importlib.import_module("07_draw")     # walk(), Canvas and the branch styles
walk, Canvas, STYLE = d7.walk, d7.Canvas, d7.STYLE
SUPPORT_MIN, BACKBONE_STEP = d7.SUPPORT_MIN, d7.BACKBONE_STEP
plt = d7.plt

RIM = 5.0                   # inches, centre to the rim of the pie
INNER = 0.95                # inches, radius of the circle the roots sit on
GAP_ROWS = 3.5              # the cut between two slices, in tip rows
MIN_ROWS = 4.0              # narrowest slice, so a one sequence fold is visible
NAME_MAX_TIPS = 60          # slices up to this size are "small": wider rows,
NAME_ROWS = 1.35            # this many per tip, and tip names with --tip-names
NAME_CHARS = 24
BIG_COMPONENT, BIG_WEIGHT = d7.BIG_COMPONENT, d7.BIG_WEIGHT
LABEL_MIN_TIPS = d7.LABEL_MIN_TIPS

# Light tints, one per fold. They sit behind substrate coloured tips, so they
# are kept pale and neighbouring slices are given hues far apart. The largest
# slice is mostly PET blue tips and gets a warm tint for contrast.
FOLD_TINT = {
    "3.40.50.1820": "#f3ead6",                              # alpha/beta hydrolase
    "3.60.70.12": "#f5dde2",                                # L-amino peptidase
    "3.90.1300.10": "#e4dff3",                              # amidase signature
    "domain:oxygenase MpaB family protein": "#dcedd6",
    "3.40.50.200": "#d8e6f5",                               # peptidase S8/S53
    "2.40.10.10": "#fbe2cf",                                # trypsin like
    "3.40.710.10": "#d3ece8",                               # beta-lactamase like
    "1.10.760.10": "#f4efc2",                               # cytochrome c like
    "2.60.40.420": "#eddcee",                               # cupredoxin
    "3.40.50.1110": "#e2e8cc",                              # SGNH hydrolase
}
SINGLE_TINT = "#e9e9e6"


def shade(hex_, k):
    """Same hue, lightness scaled by k."""
    r, g, b = (int(hex_[i:i + 2], 16) / 255 for i in (1, 3, 5))
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    r, g, b = colorsys.hls_to_rgb(h, max(0.0, min(1.0, l * k)), s)
    return "#%02x%02x%02x" % (round(r * 255), round(g * 255), round(b * 255))


def spread(want, heights):
    """Move labels as little as possible so that none overlap.

    `want` is the y each label would like, top to bottom. Pushing only
    downward piles a crowded run below where it belongs and its leader lines
    cross. This is the least squares answer instead (pool adjacent violators),
    so a crowded run is centred on where its labels wanted to be."""
    if not want:
        return []
    need = [0.0]
    for i in range(1, len(want)):
        need.append(need[-1] + (heights[i - 1] + heights[i]) / 2)
    v = [-(y) - c for y, c in zip(want, need)]       # must end up non-decreasing
    blocks = []                                       # [sum, count]
    for x in v:
        blocks.append([x, 1])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            s_, n_ = blocks.pop()
            blocks[-1][0] += s_
            blocks[-1][1] += n_
    fit = [b[0] / b[1] for b in blocks for _ in range(b[1])]
    return [-(f + c) for f, c in zip(fit, need)]


def place(items, r, cx, cy):
    """Put labels in a left and a right column near radius r, each as close
    to its own angle as the others allow.

    items: (angle, height, payload). Returns (payload, x, y, side)."""
    out = []
    for side in (1, -1):
        col = sorted((i for i in items if (math.cos(i[0]) >= 0) == (side > 0)),
                     key=lambda i: -math.sin(i[0]))
        ys = spread([cy + r * math.sin(a) for a, _, _ in col], [h for _, h, _ in col])
        for (a, h, payload), y in zip(col, ys):
            dy = max(-r, min(r, y - cy))
            x = cx + side * max(math.sqrt(max(r * r - dy * dy, 0.0)), r * 0.3)
            out.append((payload, x, y, side))
    return out


def save(cv, stem, pad=0.25):
    """Fit the canvas to what was drawn and write svg, pdf and png."""
    x0, y0, x1, y1 = cv.box
    cv.ax.set_xlim(x0 - pad, x1 + pad)
    cv.ax.set_ylim(y0 - pad, y1 + pad)
    cv.ax.set_aspect("equal")
    cv.ax.axis("off")
    cv.fig.set_size_inches(x1 - x0 + 2 * pad, y1 - y0 + 2 * pad)
    cv.fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    cv.fig.patch.set_facecolor(SURFACE)
    for ext in ("svg", "pdf", "png"):
        md = {"Date": None} if ext == "svg" else {}
        cv.fig.savefig(os.path.join(OUT, f"{stem}.{ext}"), dpi=300,
                       facecolor=SURFACE, metadata=md)
    plt.close(cv.fig)
    return x1 - x0 + 2 * pad, y1 - y0 + 2 * pad


def draw_table(table, tint):
    """The fold table as its own figure. table: (group, name, cath,
    components, sequences) per fold, in slice order."""
    cv = Canvas()
    ax = cv.ax
    row_h, fs = 0.2, 7
    c_name, c_cath, c_comp, c_seq = 0.25, 3.3, 5.3, 6.0     # inches from the left

    def cell(key, x, y, txt, ha="left", color=INK):
        cv.text(key, x, y, str(txt), fontsize=fs, color=color, ha=ha, va="center",
                zorder=8)

    def rule(key, y):
        cv.tag(ax.add_line(Line2D([0, c_seq], [y, y], color=MUTED, linewidth=0.5,
                                  zorder=8)), key)

    cell("header", c_name, 0, "fold", color=MUTED)
    cell("header", c_cath, 0, "CATH", color=MUTED)
    cell("header", c_comp, 0, "components", "right", MUTED)
    cell("header", c_seq, 0, "sequences", "right", MUTED)
    rule("header", -row_h * 0.6)
    for i, (g, name, cath, nc, nt) in enumerate(table):
        y = -(i + 1.2) * row_h
        cv.tag(ax.add_patch(Rectangle((0, y - 0.065), 0.13, 0.13, facecolor=tint(g),
                                      edgecolor=MUTED, linewidth=0.4, zorder=8)),
               f"row-{g}")
        cell(f"row-{g}", c_name, y, name)
        cell(f"row-{g}", c_cath, y, cath)
        cell(f"row-{g}", c_comp, y, nc, "right")
        cell(f"row-{g}", c_seq, y, nt, "right")
    y = -(len(table) + 1.2) * row_h
    rule("total", y + row_h * 0.6)
    cell("total", c_name, y - row_h * 0.1, "total")
    cell("total", c_comp, y - row_h * 0.1, sum(r[3] for r in table), "right")
    cell("total", c_seq, y - row_h * 0.1, sum(r[4] for r in table), "right")
    cv.grow(0, y - row_h)
    size = save(cv, "figure_pie_table", pad=0.2)

    # one group per row, all under one parent, so a row can be moved or the
    # whole table can
    svg = os.path.join(OUT, "figure_pie_table.svg")
    names = {g: name for g, name, *_ in table}
    leaf = [(f"table-{k}", f"row: {names[k[4:]]}" if k.startswith("row-") else k, v)
            for k, v in cv.gids.items()]
    regroup_svg(svg, leaf)
    regroup_svg(svg, [("fold-table", "fold table", [i for i, _, _ in leaf])])
    return size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tip-names", action="store_true",
                    help="name every tip along the rim in the small slices")
    args = ap.parse_args()
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
    comp_size = collections.Counter(t["component"] for t in tips)

    def fold_name(g):
        doms = collections.Counter(
            (rows[r]["domain"] or "NA").split(";")[0].strip()
            for t in groups[g] for r in t["records"].split(";"))
        name = " ".join(doms.most_common(1)[0][0].replace("[topology]", "").split())
        # two slices carry the beta-lactamase name; the suffix tells them apart
        note = " (unconfirmed)" if groups[g][0]["homology"] == "unconfirmed" else ""
        return name[0].upper() + name[1:] + note

    def fold_cath(g):
        key = groups[g][0]["fold_key"]
        return "none" if key.startswith("domain:") else key.replace(";", "; ")

    def tint(g):
        if len(groups[g]) == 1:
            return SINGLE_TINT
        return FOLD_TINT.get(groups[g][0]["fold_key"], SINGLE_TINT)

    order = sorted(groups, key=lambda g: (-len(groups[g]), g))
    named = {g for g in order if len(groups[g]) <= NAME_MAX_TIPS}
    trees = {g: Phylo.read(os.path.join(OUT, "trees", f"{g}.nwk"), "newick")
             for g in order}

    weight = {}
    for t in tips:
        g = t["fold_group"]
        if g in named:
            weight[t["tip"]] = NAME_ROWS
        elif comp_size[t["component"]] >= BIG_COMPONENT:
            weight[t["tip"]] = BIG_WEIGHT
        elif comp_size[t["component"]] == 1:
            weight[t["tip"]] = 2.0

    # --- one scale for every slice: the largest that fits the deepest tree ---
    reach = RIM - INNER
    scale = None
    for g in order:
        multi = n_comp[g] > 1
        if max(d["depth"] for d in walk(trees[g], multi, 0.0)) == 0:
            continue
        s = reach / max(d["depth"] for d in walk(trees[g], multi, 0.0))
        for _ in range(8):
            s = reach / max(d["depth"] for d in walk(trees[g], multi, BACKBONE_STEP / s))
        scale = s if scale is None else min(scale, s)
    step = BACKBONE_STEP / scale

    # --- rows: slices one after another, clockwise from the top -------------
    slices, at = [], GAP_ROWS / 2
    for g in order:
        nodes = walk(trees[g], n_comp[g] > 1, step, weight)
        used = max(d["hi"] for d in nodes)
        width = max(used, MIN_ROWS)
        shift = at + (width - used) / 2
        for d in nodes:
            d["y"] += shift
            d["lo"] += shift
            d["hi"] += shift
        slices.append((g, nodes, at, at + width))
        at += width + GAP_ROWS
    total = at - GAP_ROWS / 2

    cx = cy = 0.0
    ang = lambda y: math.pi / 2 - 2 * math.pi * y / total
    xy = lambda r, a: (cx + r * math.cos(a), cy + r * math.sin(a))

    def arc(r, a0, a1):
        k = max(2, int(abs(a1 - a0) / math.radians(1.0)) + 2)
        return [xy(r, a) for a in np.linspace(a0, a1, k)]

    cv = Canvas()
    ax = cv.ax
    r_arc = RIM + 0.10          # component arcs
    r_name = RIM + 0.20         # where tip names start
    name_len = NAME_CHARS * 5 * 0.0082 + 0.1
    fold_items, seed_items = [], []

    for g, nodes, lo, hi in slices:
        base = tint(g)
        a_lo, a_hi = ang(lo), ang(hi)
        cv.tag(ax.add_patch(Wedge((cx, cy), RIM, math.degrees(a_hi), math.degrees(a_lo),
                                  width=RIM - INNER, facecolor=base,
                                  edgecolor="none", zorder=0)), g)
        for d in nodes:
            d["r"] = INNER + d["depth"] * scale
            d["a"] = ang(d["y"])
        single = len(nodes) == 1
        if single:
            nodes[0]["r"] = INNER + BACKBONE_STEP
            cv.line(g, "solid", [xy(INNER, nodes[0]["a"]), xy(nodes[0]["r"], nodes[0]["a"])])

        # components inside the slice: second tint on every other one, and an
        # arc outside the rim for each
        stems = sorted((d for d in nodes if d["name"].startswith("S_")),
                       key=lambda d: d["lo"])
        for i, d in enumerate(stems):
            a0, a1 = ang(d["lo"]), ang(d["hi"])
            if i % 2:
                cv.tag(ax.add_patch(Wedge((cx, cy), RIM, math.degrees(a1),
                                          math.degrees(a0), width=RIM - INNER,
                                          facecolor=shade(base, 0.955),
                                          edgecolor="none", zorder=0.5)), g)
            pad = 0.15 * 2 * math.pi / total
            cv.line(g, (MUTED, 1.2, 4), arc(r_arc, a0 - pad, a1 + pad))
            if d["ntips"] >= LABEL_MIN_TIPS and g not in named:
                am = (a0 + a1) / 2
                x, y = xy(r_arc + 0.09, am)
                cv.text(g, x, y, d["comp"], fontsize=6, color=MUTED, zorder=8,
                        ha="left" if math.cos(am) > 0.3 else
                        "right" if math.cos(am) < -0.3 else "center",
                        va="bottom" if math.sin(am) > 0.3 else
                        "top" if math.sin(am) < -0.3 else "center")

        pts = []
        for d in nodes:
            if d["parent"] is not None:
                cv.line(g, d["edge"], [xy(d["parent"]["r"], d["a"]), xy(d["r"], d["a"])])
            if not d["tip"]:
                cv.line(g, d["fork"], arc(d["r"], d["kids"][0]["a"], d["kids"][-1]["a"]))
            if d["boot"] is not None:
                x, y = xy(d["r"], d["a"])
                cv.text(g, x, y, str(d["boot"]), fontsize=5.5, color=INK, zorder=9,
                        ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.12", fc=base, ec="none"))
            if not d["tip"]:
                continue
            t = tipinfo[d["name"]]
            x, y = xy(d["r"], d["a"])
            pts.append((x, y, t["shape"], SUBSTRATE_COLOR[t["substrate_group"]],
                        d["name"] in seeds))
            if args.tip_names and g in named:
                # a hairline from the tip to its name, so a short branch is
                # still tied to the right label
                cv.line(g, ("#b9b9b4", 0.3, 1.5), [xy(d["r"] + 0.06, d["a"]),
                                                  xy(r_name - 0.04, d["a"])])
                name = t["enzyme_name"]
                name = name[:NAME_CHARS] + ("..." if len(name) > NAME_CHARS else "")
                left = math.cos(d["a"]) < 0
                x, y = xy(r_name, d["a"])
                cv.tag(ax.text(x, y, name, fontsize=5, color=INK, zorder=8,
                               rotation=math.degrees(d["a"]) + (180 if left else 0),
                               rotation_mode="anchor", va="center",
                               ha="right" if left else "left"), g)
                cv.grow(*xy(r_name + name_len, d["a"]))
            elif d["name"] in seeds:        # no name on the rim: call it out
                seed_items.append((d["a"], 0.17, (g, d, seeds[d["name"]])))
        cv.markers(g, pts)

        # the leader to the fold name starts where this slice's own tip names end
        longest = max(min(len(tipinfo[d["name"]]["enzyme_name"]), NAME_CHARS + 3)
                      for d in nodes if d["tip"])
        r_from = (r_name + longest * 5 * 0.0082 + 0.08
                  if args.tip_names and g in named else r_arc + 0.05)
        fold_items.append(((a_lo + a_hi) / 2, 0.32, (g, r_from)))

    # --- seed names, wherever the tips are not named ------------------------
    r_seed = r_arc + 0.5
    for (g, d, name), lx, ly, side in place(seed_items, r_seed, cx, cy):
        cv.line(g, (MUTED, 0.5, 6), [xy(d["r"] + 0.05, d["a"]),
                                     xy(r_arc + 0.1, d["a"]), (lx, ly)])
        cv.text(g, lx + side * 0.04, ly, name, fontsize=7.5, color=INK, zorder=9,
                ha="left" if side > 0 else "right", va="center")

    # --- fold names, one per slice, outside everything -----------------------
    # Name only. CATH accession and the counts are in the table. Each label
    # (leader, swatch, name) is its own SVG group, apart from its slice.
    r_fold = r_name + name_len + 0.45 if args.tip_names else r_arc + 1.3
    for (g, r_from), lx, ly, side in place(fold_items, r_fold, cx, cy):
        key = f"label-{g}"
        am = next((lo + hi) / 2 for gg, _, lo, hi in slices if gg == g)
        cv.line(key, (MUTED, 0.5, 6), [xy(r_from, ang(am)), (lx, ly)])
        cv.tag(ax.add_patch(Rectangle((lx + side * 0.05 - (0.15 if side < 0 else 0),
                                       ly - 0.075), 0.15, 0.15, facecolor=tint(g),
                                      edgecolor=MUTED, linewidth=0.4, zorder=7)), key)
        cv.text(key, lx + side * 0.27, ly, fold_name(g), fontsize=9, color=INK,
                zorder=9, ha="left" if side > 0 else "right", va="center")
    cv.flush()

    # --- legends in a strip to the right -------------------------------------
    # top right corner of the pie's bounding square, which the circle leaves free
    x_leg, y_top = r_fold + 0.2, cv.box[3]
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
            Line2D([], [], color=MUTED, linewidth=1.2, label="one SSN component"),
            Line2D([], [], color=FOLD_TINT["3.40.50.1820"], linewidth=7,
                   label="slice background: fold\n(two tints: neighbouring components)")]

    def legend(key, handles, ttl, y, ncol=1):
        lg = ax.legend(handles=handles, title=ttl, loc="upper left",
                       bbox_to_anchor=(x_leg, y), bbox_transform=ax.transData,
                       frameon=False, fontsize=7.5, ncol=ncol, handletextpad=0.5,
                       columnspacing=1.0, labelspacing=0.5, borderaxespad=0.0,
                       alignment="left")
        lg.get_title().set_fontsize(7.5)
        lg.get_title().set_color(MUTED)
        ax.add_artist(lg)
        cv.tag(lg, key)

    legend("legend-branches", h_br, "branches and slices", y_top)
    legend("legend-substrate", h_sub, "substrate", y_top - 1.75, ncol=2)
    legend("legend-origin", h_shape, "origin", y_top - 3.25)
    bar, by = 0.5, y_top - 4.55
    cv.tag(ax.add_line(Line2D([x_leg, x_leg + bar * scale], [by, by], color=INK,
                              linewidth=1.0, zorder=8)), "legend-scale")
    cv.text("legend-scale", x_leg, by + 0.06, f"{bar} substitutions per site",
            fontsize=7.5, color=MUTED, ha="left", va="bottom")
    cv.grow(x_leg + 2.6, y_top)

    w, h = save(cv, "figure_pie")

    # --- fold table: a tsv and a figure of its own ---------------------------
    table = [(g, fold_name(g), fold_cath(g), n_comp[g], len(groups[g])) for g in order]
    write_tsv(os.path.join(OUT, "fold_table.tsv"),
              ["fold_group", "fold", "cath", "components", "sequences", "records",
               "homology"],
              [[g, name, cath.replace("; ", ";"), nc, nt,
                sum(len(t["records"].split(";")) for t in groups[g]),
                groups[g][0]["homology"]] for g, name, cath, nc, nt in table])
    tw, th = draw_table(table, tint)

    # One named group per slice, per fold label and per legend.
    # A second pass nests the fold labels under one parent and the legends
    # under another, so either can be moved whole or piece by piece.
    svg = os.path.join(OUT, "figure_pie.svg")
    leaf = []
    for k, v in sorted(cv.gids.items()):
        if k.startswith("label-"):
            g = k[6:]
            leaf.append((f"fold-label-{g}", f"label: {fold_name(g)}", v))
        elif k.startswith("legend-"):
            leaf.append((k, f"legend: {k[7:]}", v))
        else:
            leaf.append((f"slice-{k}", f"{k} {fold_name(k)}", v))
    moved = regroup_svg(svg, leaf)
    regroup_svg(svg, [
        ("fold-labels", "fold labels", [i for i, _, _ in leaf if i.startswith("fold-label-")]),
        ("legend", "legend", [i for i, _, _ in leaf if i.startswith("legend-")])])
    share = {g: (hi - lo) / total for g, _, lo, hi in slices}
    print(f"{len(slices)} slices, {total:.0f} rows; largest slice "
          f"{100 * share[order[0]]:.0f}% of the circle; "
          f"tip names {'on' if args.tip_names else 'off'}, "
          f"{len(seed_items)} seeds called out")
    print(f"scale {scale:.2f} in per substitution per site; "
          f"canvas {w:.1f} x {h:.1f} in; "
          f"svg: {len(leaf)} groups, {moved} elements grouped")
    print(f"table: figure_pie_table.svg, {tw:.1f} x {th:.1f} in, {len(table)} rows")


if __name__ == "__main__":
    main()
