#!/usr/bin/env python3
"""
Step 5: the fold figures of step 4, each one as a heatmap.

  figure_heat_cutoffs    fold by identity cut-off, clusters per cell
                         (the step 4 heatmap, copied under this name)
  figure_heat_families   fold by SSN family rank, centroids per cell
                         (the ranked bars as a heatmap)
  figure_heat_tiles      one tile per centroid, blocked by fold, coloured by
                         its identity to the nearest other centroid of its fold
                         (the treemap as a heatmap)
  figure_heat_identity   fold by identity bin, share of the pairs inside the
                         fold (the identity distributions as a heatmap)
  figure_heat_substrate  fold by substrate, centroids per cell
                         (what the pie shows with tip colour, as a heatmap)

Rows are the 16 folds in the same order in every figure, largest first.
Data and cluster definitions come from step 4.
"""

import collections
import importlib.util
import os
import shutil
import sys

os.environ.setdefault("PLASTICOME_TREE_SET", "centroids")

import numpy as np
from matplotlib.colors import to_hex

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "fold_diversity", os.path.join(HERE, "04_fold_diversity.py"))
fd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fd)

from common import INK, MUTED, SUBSTRATE_ORDER  # noqa: E402  (path set by step 4)

hm, Canvas, RAMP, OUT, UNIT, FS = fd.hm, fd.Canvas, fd.RAMP, fd.OUT, fd.UNIT, fd.FS
EMPTY = "#f1f0ec"       # a cell with nothing in it
CW, CH, GAP = 0.50, 0.27, 0.025


def row_labels(cv, folds):
    for i, f in enumerate(folds):
        cv.text(-0.10, -(i + 0.5) * CH, f["name"], ha="right")


def cell(cv, i, j, color, label=None, light=False, w=CW):
    cv.rect(j * w + GAP / 2, -(i + 1) * CH + GAP / 2, w - GAP, CH - GAP, color)
    if label is not None:
        cv.text((j + 0.5) * w, -(i + 0.5) * CH, label, ha="center",
                color="white" if light else INK)


def ramp_key(cv, x, y, title, lo, hi, width=1.3):
    cv.text(x, y + 0.34, title, va="bottom", linespacing=1.2)
    steps = 40
    for s in range(steps):
        cv.rect(x + s * width / steps, y + 0.14, width / steps + 0.002, 0.14,
                RAMP(s / (steps - 1)))
    cv.text(x, y + 0.02, lo)
    cv.text(x + width, y + 0.02, hi, ha="right")


def log_shade(k, top):
    """Position on the ramp for a count, log scaled so 1 is still visible."""
    return 0.12 + 0.88 * np.log(k) / np.log(top) if top > 1 else 1.0


# ------------------------------------------------------------------ figures
def heat_families(folds, comps_of):
    cv = Canvas()
    w = 0.27                                     # square cells
    ncol = max(len(comps_of[f["id"]]) for f in folds)
    top = max(max(comps_of[f["id"]]) for f in folds)
    row_labels(cv, folds)
    for i, f in enumerate(folds):
        sizes = comps_of[f["id"]]
        for j in range(ncol):
            if j < len(sizes):
                t = log_shade(sizes[j], top)
                cell(cv, i, j, RAMP(t), str(sizes[j]), light=t > 0.55, w=w)
            else:
                cell(cv, i, j, EMPTY, w=w)
        cv.text(ncol * w + 0.12, -(i + 0.5) * CH,
                f"{len(sizes)} {'family' if len(sizes) == 1 else 'families'}, "
                f"{f['n']} {UNIT if f['n'] > 1 else UNIT[:-1]}", color=MUTED)
    for j in range(ncol):
        if j == 0 or (j + 1) % 5 == 0:
            cv.text((j + 0.5) * w, 0.06, str(j + 1), ha="center", va="bottom")
    cv.text(ncol * w / 2, 0.26,
            f"SSN families of each fold, largest first. Cell: {UNIT} in the family",
            ha="center", va="bottom", size=8)
    ramp_key(cv, ncol * w + 2.3, -len(folds) * CH * 0.55,
             f"{UNIT} in the family\n(log scale)", "1", str(top))
    cv.save("figure_heat_families")


def heat_tiles(d):
    folds, m = d["folds"], d["m"]
    cv = Canvas()
    T, rows_n, g = 0.17, 12, 0.02
    # nearest aligned neighbour inside the fold; an alignment across folds
    # is chance, not a relative
    same = d["fold_of"][:, None] == d["fold_of"][None, :]
    nn = np.where(same, m, 0).max(axis=1)
    lo = np.array([b[0] for b in hm.BINS])
    comp_size = collections.Counter(d["comp_of"])
    x = depth = 0.0
    left = 0.0
    cv.text(0, 0.12, f"One tile per centroid ({d['n_seq']}), blocked by fold. "
            "Colour: identity to the nearest other centroid of its fold",
            va="bottom", size=8)
    for f in folds:
        idx = [i for i in range(d["n_seq"]) if d["fold_of"][i] == f["id"]]
        # families together, largest first, then tree order
        idx.sort(key=lambda i: (-comp_size[d["comp_of"][i]], d["comp_of"][i],
                                int(d["tips"][i]["order_in_tree"])))
        for k, i in enumerate(idx):
            c, r = divmod(k, rows_n)
            col = (hm.NO_ALN_SAME_FOLD if nn[i] == 0 else
                   hm.BINS[max(0, np.searchsorted(lo, nn[i], side="right") - 1)][1])
            cv.rect(x + c * T + g / 2, -(r + 1) * T + g / 2, T - g, T - g, col)
        ncol = -(-len(idx) // rows_n)
        lab = f"{fd.short(f['name'])}  ({f['n']})"
        ax_ = x + ncol * T / 2 + 0.03
        cv.ax.text(ax_, -rows_n * T - 0.08, lab, ha="right", va="top",
                   rotation=50, rotation_mode="anchor", fontsize=FS, color=INK)
        # a rotated label reaches down and to the left of its anchor
        wl = cv.width(lab)
        depth = max(depth, wl * np.sin(np.radians(50)) + 0.12)
        left = min(left, ax_ - wl * np.cos(np.radians(50)) - 0.05)
        x += ncol * T + 0.10
    # the rotated labels hang below; reserve their depth
    cv.see(left, -rows_n * T - 0.08 - depth, x, 0)
    L = 0.24
    lx, ly = x + 0.25, -0.9
    cv.text(lx, ly + 0.42, "identity to the nearest other\ncentroid of the fold (aa %)",
            va="bottom", linespacing=1.2)
    for k, (b, col) in enumerate(hm.BINS):
        cv.rect(lx + k * L, ly, L - 0.03, L - 0.03, col)
        if k:
            cv.text(lx + k * L, ly - 0.04, str(b), ha="center", va="top")
    cv.text(lx + len(hm.BINS) * L - 0.03, ly - 0.04, "100", ha="center", va="top")
    cv.rect(lx, ly - 0.62, L - 0.03, L - 0.03, hm.NO_ALN_SAME_FOLD)
    cv.text(lx + L + 0.03, ly - 0.62 + L / 2 - 0.015, "aligns with no other centroid of its fold")
    cv.save("figure_heat_tiles")
    return nn


def heat_identity(d):
    folds, m = d["folds"], d["m"]
    cv = Canvas()
    lo = np.array([b[0] for b in hm.BINS])
    heads = ["not\naligned", "<30"] + [f"{b}-{b + 10}" for b, _ in hm.BINS[1:]]
    row_labels(cv, folds)
    for i, f in enumerate(folds):
        idx = np.where(d["fold_of"] == f["id"])[0]
        if len(idx) < 2:
            for j in range(len(heads)):
                cell(cv, i, j, EMPTY)
            cv.text(len(heads) * CW + 0.12, -(i + 0.5) * CH, "one sequence",
                    color=MUTED)
            continue
        v = m[np.ix_(idx, idx)][np.triu_indices(len(idx), 1)]
        b = np.where(v > 0, np.clip(np.searchsorted(lo, v, side="right") - 1,
                                    0, len(lo) - 1) + 1, 0)
        share = np.bincount(b, minlength=len(heads)) / len(v)
        for j, s in enumerate(share):
            if s == 0:
                cell(cv, i, j, EMPTY)
            else:
                cell(cv, i, j, RAMP(s), "<1" if s < 0.005 else f"{100 * s:.0f}",
                     light=s > 0.55)
        cv.text(len(heads) * CW + 0.12, -(i + 0.5) * CH,
                f"{len(v):,} {'pair' if len(v) == 1 else 'pairs'}", color=MUTED)
    for j, h in enumerate(heads):
        cv.text((j + 0.5) * CW, 0.06, h, ha="center", va="bottom", linespacing=1.1)
    cv.text(len(heads) * CW / 2, 0.42,
            "Identity of the pairs inside each fold (aa %). Cell: % of the fold's pairs",
            ha="center", va="bottom", size=8)
    ramp_key(cv, len(heads) * CW + 1.3, -len(folds) * CH * 0.55,
             "share of the\nfold's pairs", "0%", "100%")
    cv.save("figure_heat_identity")


def heat_substrate(d):
    folds = d["folds"]
    cv = Canvas()
    sub = np.array([t["substrate_group"] for t in d["tips"]])
    total = collections.Counter(sub)
    cols = sorted((s for s in SUBSTRATE_ORDER if total[s]),
                  key=lambda s: (s in ("NR", "Combi."), -total[s]))
    top = max(collections.Counter(zip(d["fold_of"], sub)).values())
    row_labels(cv, folds)
    for i, f in enumerate(folds):
        here = collections.Counter(sub[d["fold_of"] == f["id"]])
        for j, s in enumerate(cols):
            k = here[s]
            if k:
                t = log_shade(k, top)
                cell(cv, i, j, RAMP(t), str(k), light=t > 0.55)
            else:
                cell(cv, i, j, EMPTY)
        cv.text(len(cols) * CW + 0.12, -(i + 0.5) * CH,
                f"{sum(1 for s in cols if here[s])} of {len(cols)}", color=MUTED)
    for j, s in enumerate(cols):
        cv.text((j + 0.5) * CW, 0.06, s, ha="center", va="bottom")
    cv.text(len(cols) * CW + 0.12, 0.06, "substrates", va="bottom", color=MUTED)
    cv.text(len(cols) * CW / 2, 0.26, f"Substrate. Cell: {UNIT} of the fold",
            ha="center", va="bottom", size=8)
    y = -len(folds) * CH - 0.10
    cv.line([0, len(cols) * CW], [y, y])
    first = collections.Counter(sub[d["fold_of"] == folds[0]["id"]])
    for lab, fn in [("all folds", lambda s: str(total[s])),
                    (f"share in {folds[0]['name']}",
                     lambda s: f"{100 * first[s] / total[s]:.0f}%")]:
        y -= CH * 0.85
        cv.text(-0.10, y + CH / 2, lab, ha="right")
        for j, s in enumerate(cols):
            cv.text((j + 0.5) * CW, y + CH / 2, fn(s), ha="center")
    ramp_key(cv, len(cols) * CW + 1.0, -len(folds) * CH * 0.55,
             f"{UNIT}\n(log scale)", "1", str(top))
    cv.save("figure_heat_substrate")
    return cols, total, first


def main():
    d = fd.load_data()
    os.makedirs(OUT, exist_ok=True)
    for ext in ("svg", "pdf", "png"):
        src = os.path.join(OUT, f"figure_fold_cutoffs.{ext}")
        if not os.path.exists(src):
            raise SystemExit("run 04_fold_diversity.py first")
        shutil.copyfile(src, os.path.join(OUT, f"figure_heat_cutoffs.{ext}"))
    heat_families(d["folds"], d["comps_of"])
    nn = heat_tiles(d)
    heat_identity(d)
    cols, total, first = heat_substrate(d)

    hyd = d["fold_of"] == d["folds"][0]["id"]
    print(f"set: {fd.SET}. {d['n_seq']} {UNIT}")
    print(f"  no aligned neighbour: {int((nn == 0).sum())} "
          f"({int(((nn == 0) & hyd).sum())} in {d['folds'][0]['name']})")
    for name, mask in (("hydrolase", hyd), ("other folds", ~hyd)):
        v = nn[mask & (nn > 0)]
        print(f"  nearest neighbour identity, {name}: median {np.median(v):.0f}")
    print("  substrates: " + ", ".join(
        f"{s} {first[s]}/{total[s]}" for s in cols))


if __name__ == "__main__":
    main()
