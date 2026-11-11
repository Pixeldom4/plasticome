#!/usr/bin/env python3
"""
Step 4: where the sequence diversity of the plasticome sits, fold by fold.

Four figures that make the same point in different forms: most distinct
sequences are alpha/beta hydrolases and the other folds are a long tail.

  figure_fold_cutoffs    heatmap, fold by identity cut-off, clusters per cell
  figure_fold_bars       ranked bars, centroids per fold, cut by SSN component
  figure_fold_treemap    treemap, area is centroids, fold then SSN component
  figure_fold_identity   identity of the pairs inside each fold, as shares

Clusters at a cut-off
  30   the SSN components of ../ssn, the family definition of the project
  40 to 80   single linkage inside each SSN component, on the pairs of step 1
       that reach the cut-off and cover half of the shorter sequence
  90   the sequences themselves, since the set is the 90 percent centroids
A cluster never crosses a component, so the counts can only grow with the
cut-off and a cluster always has one fold.

The default set is the centroids. Writes fold_cutoffs.tsv beside the figures.
"""

import collections
import importlib.util
import os
import sys
import textwrap

os.environ.setdefault("PLASTICOME_TREE_SET", "centroids")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(HERE, "..", "tree")
sys.path.insert(0, TREE)
from common import (INK, MUTED, OUT as TREE_OUT, SET, SUFFIX, SURFACE,  # noqa: E402
                    WORK_ALL, load_rows, read_tsv, write_tsv)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


hm = load(os.path.join(HERE, "02_heatmap.py"), "heatmap_groups")
pie = load(os.path.join(TREE, "08_draw_pie.py"), "draw_pie")

OUT = os.path.join(HERE, "out" + SUFFIX)
CUTS = [30, 40, 50, 60, 70, 80, 90]
FS = 6.6
BLUE = "#2a78d6"        # the fold the figure is about
GREY = "#b5b3ac"        # every other fold
RAMP = LinearSegmentedColormap.from_list(
    "share", ["#eef4fc", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
UNIT = "centroids" if SET == "centroids" else "sequences"


class Canvas:
    """A figure drawn in inches, fitted to its content when saved."""

    def __init__(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.rend = self.fig.canvas.get_renderer()
        self.box = [0.0, 0.0, 0.0, 0.0]

    def see(self, x0, y0, x1, y1):
        b = self.box
        self.box = [min(b[0], x0), min(b[1], y0), max(b[2], x1), max(b[3], y1)]

    def rect(self, x, y, w, h, color, **kw):
        self.ax.add_patch(Rectangle((x, y), w, h, facecolor=color,
                                    edgecolor="none", **kw))
        self.see(x, y, x + w, y + h)

    def text(self, x, y, s, ha="left", va="center_baseline", size=FS,
             color=INK, **kw):
        t = self.ax.text(x, y, s, ha=ha, va=va, fontsize=size, color=color, **kw)
        bb = t.get_window_extent(self.rend)
        w, h = bb.width / self.fig.dpi, bb.height / self.fig.dpi
        x0 = x - w if ha == "right" else x - w / 2 if ha == "center" else x
        self.see(x0, y - h, x0 + w, y + h)
        return w, h

    def width(self, s, size=FS):
        t = self.ax.text(0, 0, s, fontsize=size)
        w = t.get_window_extent(self.rend).width / self.fig.dpi
        t.remove()
        return w

    def line(self, xs, ys, color=MUTED, lw=0.6):
        self.ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="butt")
        self.see(min(xs), min(ys), max(xs), max(ys))

    def save(self, stem, pad=0.2):
        x0, y0, x1, y1 = self.box
        self.ax.set_xlim(x0 - pad, x1 + pad)
        self.ax.set_ylim(y0 - pad, y1 + pad)
        self.ax.set_aspect("equal")
        self.ax.axis("off")
        self.fig.set_size_inches(x1 - x0 + 2 * pad, y1 - y0 + 2 * pad)
        self.fig.patch.set_facecolor(SURFACE)
        for ext in ("svg", "pdf", "png"):
            md = {"Date": None} if ext == "svg" else {}
            self.fig.savefig(os.path.join(OUT, f"{stem}.{ext}"), dpi=300,
                             facecolor=SURFACE, metadata=md)
        plt.close(self.fig)


def squarify(areas, x, y, w, h):
    """Squarified treemap (Bruls et al.). `areas` largest first, summing to
    w * h. Returns one (x, y, w, h) per area, y measured downward."""
    out, areas = [], list(areas)
    while areas:
        side = min(w, h)

        def worst(row):
            s = sum(row)
            return max(max(row) * side * side / (s * s),
                       s * s / (min(row) * side * side))

        k = 1
        while k < len(areas) and worst(areas[:k + 1]) <= worst(areas[:k]):
            k += 1
        row, areas = areas[:k], areas[k:]
        s = sum(row)
        if w >= h:
            cw, yy = s / h, y
            for a in row:
                out.append((x, yy, cw, a / cw))
                yy += a / cw
            x, w = x + cw, w - cw
        else:
            rh, xx = s / w, x
            for a in row:
                out.append((xx, y, a / rh, rh))
                xx += a / rh
            y, h = y + rh, h - rh
    return out


def pct(x):
    return "<1%" if x < 0.005 else f"{100 * x:.0f}%"


def short(name):
    return (name.replace(" superfamily", "").replace(" - blue copper proteins", "")
            .replace(" family protein", "").replace(", domain 2", "")
            .replace(" domain", ""))


# ------------------------------------------------------------------ figures
def fig_cutoffs(folds, counts):
    cv = Canvas()
    cw, ch, gap = 0.50, 0.27, 0.025
    tot = {c: sum(counts[f["id"]][c] for f in folds) for c in CUTS}
    top = max(counts[f["id"]][c] / tot[c] for f in folds for c in CUTS)
    first = counts[folds[0]["id"]]
    for i, f in enumerate(folds):
        y = -(i + 1) * ch
        cv.text(-0.10, y + ch / 2, f["name"], ha="right")
        for j, c in enumerate(CUTS):
            k = counts[f["id"]][c]
            share = k / tot[c]
            cv.rect(j * cw + gap / 2, y + gap / 2, cw - gap, ch - gap,
                    RAMP(share / top))
            cv.text((j + 0.5) * cw, y + ch / 2, str(k), ha="center",
                    color="white" if share / top > 0.55 else INK)
    n, w = len(folds), len(CUTS) * cw
    for j, c in enumerate(CUTS):
        cv.text((j + 0.5) * cw, 0.06, str(c), ha="center", va="bottom")
    cv.text(0.5 * cw, 0.22, "SSN\nfamilies", ha="center", va="bottom",
            color=MUTED, linespacing=1.15)
    cv.text((len(CUTS) - 0.5) * cw, 0.22, UNIT, ha="center", va="bottom",
            color=MUTED)
    cv.text(w / 2, 0.52, "Clusters at identity cut-off (aa %)", ha="center",
            va="bottom", size=8)

    y = -n * ch - 0.10
    cv.line([0, w], [y, y])
    for lab, fn in [("all folds", lambda c: str(tot[c])),
                    (f"share in {folds[0]['name']}",
                     lambda c: f"{100 * first[c] / tot[c]:.0f}%")]:
        y -= ch * 0.85
        cv.text(-0.10, y + ch / 2, lab, ha="right")
        for j, c in enumerate(CUTS):
            cv.text((j + 0.5) * cw, y + ch / 2, fn(c), ha="center")

    # the last column again as length, so the tail can be seen
    bx, bw = w + 0.25, 2.0
    cv.text(bx, 0.06, UNIT.capitalize(), va="bottom")
    for i, f in enumerate(folds):
        cv.rect(bx, -(i + 1) * ch + 0.06, bw * f["n"] / folds[0]["n"],
                ch - 0.12, BLUE if i == 0 else GREY)

    # colour key
    kx, ky = bx + 0.55, -n * ch * 0.55
    cv.text(kx, ky + 0.34, "share of the clusters\nat that cut-off", va="bottom",
            linespacing=1.2)
    steps = 40
    for s in range(steps):
        cv.rect(kx + s * 1.3 / steps, ky + 0.14, 1.3 / steps + 0.002, 0.14,
                RAMP(s / (steps - 1)))
    cv.text(kx, ky + 0.02, "0%", ha="left")
    cv.text(kx + 1.3, ky + 0.02, f"{100 * top:.0f}%", ha="right")
    cv.save("figure_fold_cutoffs")


def fig_bars(folds, comps_of, n_seq):
    cv = Canvas()
    rh, bw = 0.27, 5.2
    big = max(f["n"] for f in folds)
    cv.text(0, 0.10, f"{UNIT.capitalize()} per fold. Each block is one SSN family.",
            va="bottom", size=8)
    for i, f in enumerate(folds):
        y = -(i + 1) * rh
        cv.text(-0.10, y + rh / 2, f["name"], ha="right")
        x = 0.0
        for k in comps_of[f["id"]]:
            w = bw * k / big
            cv.rect(x, y + 0.05, max(w - 0.012, 0.012), rh - 0.10,
                    BLUE if i == 0 else GREY)
            x += w
        nf = len(comps_of[f["id"]])
        cv.text(x + 0.06, y + rh / 2,
                f"{f['n']}  ({pct(f['n'] / n_seq)}),  "
                f"{nf} {'family' if nf == 1 else 'families'}")
    cv.save("figure_fold_bars")


def fig_treemap(folds, comps_of, n_seq):
    cv = Canvas()
    W, H = 7.2, 4.6
    scale = W * H / n_seq
    key = []
    for f, (x, y, w, h) in zip(folds, squarify([f["n"] * scale for f in folds],
                                               0, 0, W, H)):
        tint = pie.FOLD_TINT.get(f["key"], pie.SINGLE_TINT)
        g = 0.035                       # gap between folds
        ix, iy, iw, ih = x + g / 2, y + g / 2, w - g, h - g
        sub = comps_of[f["id"]]
        s2 = iw * ih / sum(sub)
        for k, (cx, cy, c_w, c_h) in zip(sub, squarify([k * s2 for k in sub],
                                                       ix, iy, iw, ih)):
            cg = 0.012 if len(sub) > 1 else 0.0
            cv.rect(cx + cg / 2, -(cy + c_h) + cg / 2, max(c_w - cg, 0.004),
                    max(c_h - cg, 0.004), pie.shade(tint, 0.90))
        lab = None
        for chars in range(60, 7, -4):
            lines = textwrap.wrap(short(f["name"]), chars,
                                  break_long_words=False) + [str(f["n"])]
            if (max(cv.width(l) for l in lines) < iw - 0.12
                    and len(lines) * FS * 1.25 / 72 < ih - 0.10):
                lab = "\n".join(lines)
                break
        if lab:
            cv.text(ix + 0.06, -iy - 0.06, lab, va="top", linespacing=1.2)
        else:
            key.append(f)
            cv.text(ix + iw / 2, -iy - ih / 2, str(len(key)), ha="center")
    cv.text(0, 0.10, f"Area is the number of {UNIT} ({n_seq}). "
            "Blocks inside a fold are its SSN families.", va="bottom", size=8)
    for k, f in enumerate(key):
        cv.text(W + 0.18, -0.12 - k * 0.20,
                f"{k + 1}   {short(f['name'])}  ({f['n']})")
    cv.save("figure_fold_treemap")


def fig_identity(folds, m, fold_of):
    cv = Canvas()
    rh, bw = 0.30, 5.0
    lo = np.array([b[0] for b in hm.BINS])
    cols = [hm.NO_ALN_SAME_FOLD] + [b[1] for b in hm.BINS]
    cv.text(0, 0.10, "Identity of the pairs inside each fold, as a share of its pairs",
            va="bottom", size=8)
    row = 0
    singles = []
    for f in folds:
        idx = np.where(fold_of == f["id"])[0]
        if len(idx) < 2:
            singles.append(f["name"])
            continue
        v = m[np.ix_(idx, idx)][np.triu_indices(len(idx), 1)]
        b = np.where(v > 0, np.clip(np.searchsorted(lo, v, side="right") - 1,
                                    0, len(lo) - 1) + 1, 0)
        share = np.bincount(b, minlength=len(cols)) / len(v)
        y = -(row + 1) * rh
        cv.text(-0.10, y + rh / 2, f["name"], ha="right")
        x = 0.0
        for s, col in zip(share, cols):
            if s > 0:
                cv.rect(x, y + 0.05, max(bw * s - 0.012, 0.006), rh - 0.10, col)
            x += bw * s
        cv.text(bw + 0.08, y + rh / 2,
                f"{len(v):,} {'pair' if len(v) == 1 else 'pairs'}", color=MUTED)
        row += 1
    y = -row * rh - 0.08
    for t in (0, 0.25, 0.5, 0.75, 1.0):
        cv.line([bw * t, bw * t], [y, y - 0.04])
        cv.text(bw * t, y - 0.07, f"{100 * t:.0f}%", ha="center", va="top")
    # legend
    L = 0.24
    lx, ly = 0.0, y - 0.62
    cv.rect(lx, ly, L - 0.03, L - 0.03, cols[0])
    cv.text(lx + L + 0.03, ly + L / 2 - 0.015, "not aligned")
    x = lx + 1.25
    cv.text(x - 0.08, ly + L / 2 - 0.015, "identity (aa %)", ha="left")
    x += cv.width("identity (aa %)") + 0.15
    for k, (b, col) in enumerate(hm.BINS):
        cv.rect(x, ly, L - 0.03, L - 0.03, col)
        if k:
            cv.text(x, ly - 0.04, str(b), ha="center", va="top")
        x += L
    cv.text(x - 0.03, ly - 0.04, "100", ha="center", va="top")
    if singles:
        cv.text(0, ly - 0.42, f"Not shown, one sequence each: "
                + ", ".join(short(s) for s in singles) + ".", color=MUTED)
    cv.save("figure_fold_identity")


def load_data():
    """Everything the fold figures are drawn from. Step 5 reads it too."""
    rows = load_rows()
    full = {rows[t["tip"]]["aa_sequence"].strip().upper(): t["tip"]
            for t in read_tsv(os.path.join(WORK_ALL, "tips.tsv"))}
    tips = read_tsv(os.path.join(TREE_OUT, "tree_tips.tsv"))
    order = [full[rows[t["tip"]]["aa_sequence"].strip().upper()] for t in tips]
    m, _ = hm.load_identity(order)
    n_seq = len(tips)
    fold_of = np.array([t["fold_group"] for t in tips])
    comp_of = np.array([t["component"] for t in tips])
    name = {r["fold_group"]: r["fold"] for r in
            read_tsv(os.path.join(TREE_OUT, "fold_table.tsv"))}
    fkey = {t["fold_group"]: t["fold_key"] for t in tips}

    size = collections.Counter(fold_of)
    folds = [{"id": g, "name": name[g], "n": size[g], "key": fkey[g]}
             for g in sorted(size, key=lambda g: (-size[g], g))]
    comps_of = {}
    for g in size:
        c = collections.Counter(comp_of[fold_of == g])
        comps_of[g] = sorted(c.values(), reverse=True)

    same_comp = comp_of[:, None] == comp_of[None, :]
    by_id = {g: {} for g in size}
    for c in CUTS:
        if c == CUTS[0]:
            lab = comp_of
        elif c == CUTS[-1] and SET == "centroids":
            lab = np.arange(n_seq)
        else:
            _, lab = connected_components(csr_matrix((m >= c) & same_comp),
                                          directed=False)
        for g in size:
            by_id[g][c] = len(set(lab[fold_of == g]))
    return dict(tips=tips, m=m, n_seq=n_seq, fold_of=fold_of, comp_of=comp_of,
                folds=folds, comps_of=comps_of, clusters=by_id)


def main():
    d = load_data()
    folds, comps_of, by_id = d["folds"], d["comps_of"], d["clusters"]
    n_seq, m, fold_of, size = d["n_seq"], d["m"], d["fold_of"], d["clusters"]

    os.makedirs(OUT, exist_ok=True)
    write_tsv(os.path.join(OUT, "fold_cutoffs.tsv"),
              ["fold_group", "fold", UNIT, "ssn_components"]
              + [f"clusters_{c}" for c in CUTS],
              [[f["id"], f["name"], f["n"], len(comps_of[f["id"]])]
               + [by_id[f["id"]][c] for c in CUTS] for f in folds])

    fig_cutoffs(folds, by_id)
    fig_bars(folds, comps_of, n_seq)
    fig_treemap(folds, comps_of, n_seq)
    fig_identity(folds, m, fold_of)

    print(f"set: {SET}. {n_seq} {UNIT}, {len(folds)} folds")
    print("  cut-off  clusters  in " + folds[0]["name"])
    for c in CUTS:
        tot = sum(by_id[g][c] for g in size)
        k = by_id[folds[0]["id"]][c]
        print(f"  {c:>7}  {tot:>8}  {k:>4}  {100 * k / tot:.0f}%")
    print(f"  folds with 5 {UNIT} or fewer: {sum(f['n'] <= 5 for f in folds)}")


if __name__ == "__main__":
    main()
