#!/usr/bin/env python3
"""
Step 7: the heatmap set after the reviewer feedback of 2026-10-27.

Main figure, out_v2/figure_heatmap_v2.*
  The 411 centroids against each other, as step 3, with these changes:
  - one colour for "not alignable" (the alignment covers under half of the
    shorter sequence), whatever the folds. Step 3 used two greys, same fold
    and different fold, which read as two kinds of result. The fold is shown
    by the margins and outlines instead
  - each fold's block on the diagonal is outlined in the fold's colour, and
    the margins carry a fold strip in the same colours as the tree
    (../tree/folds.py), so the clustering of the folds is visible
  - rows follow ../tree's ring figure (ring_order.tsv), the same order as
    the circular tree
  - more labels: every fold with its counts, every SSN component of 5 or
    more, the curated seeds, and boxes for the three zooms
  - a worked legend, "how a cell gets its colour", from pair to colour

Zooms, out_v2/figure_zoom_<name>.*
  All distinct sequences (not only centroids) of one group, every row named:
  ispetase   the 50 percent subfamily of component 0 that holds IsPETase
  studied    the 50 percent subfamily with the most records, the group the
             literature has returned to most often
  dense      the 70 percent subfamily with the highest median identity among
             those of 20 or more sequences

Cluster sizes, out_v2/figure_cluster_sizes.*
  Every SSN component ranked by size, coloured by fold, and the largest 90
  percent clusters, each named.

Length distribution, out_v2/figure_lengths.*
  Plasticome lengths against PETadex-nr family representatives and the Logan
  sequences sent to the halo assay, see LENGTH_SOURCES.

Run ../tree/18_ring_tree.py for both sets first.
"""

import collections
import importlib.util
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["svg.hashsalt"] = "plasticome-heatmap-v2"
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(HERE, "..", "tree")
sys.path.insert(0, TREE)
os.environ.setdefault("PLASTICOME_TREE_SET", "centroids")
from common import INK, MUTED, SURFACE, RUNS, load_rows, read_tsv, seed_nodes, write_tsv  # noqa: E402
import folds  # noqa: E402


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


hm = load(os.path.join(HERE, "02_heatmap.py"), "heatmap_groups")

OUT = os.path.join(HERE, "out_v2")
NO_ALN = "#d6d4ce"
DPI = 300
TIERS = [("~Fold", 20, 30), ("~Family", 30, 50), ("~Subfamily", 50, 90), ("~Centroid", 90, 100)]

PROJECTS = os.path.join(HERE, "..", "..", "..")
LENGTH_SOURCES = {
    # full ORFs of the PETadex atlas family representatives (nr origin). The
    # representative sequence itself is a trimmed catalytic domain, so its
    # length is drawn too, to show what trimming does to the comparison
    "PETadex-nr families": os.path.join(
        PROJECTS, "petadex-adjacency-landscape", "data", "family_rep_orfs.csv.gz"),
    # the only Logan protein sequences on this machine: the 199 genes sent to
    # the halo assay. A chosen sample, not a draw from Logan
    "Logan, halo assay set": os.path.join(
        PROJECTS, "petai", "inputs", "halo", "logan.haloAssay.v241123.msta.fa"),
}


def ring_order(run):
    return read_tsv(os.path.join(RUNS, run, "results", "ring_order.tsv"))


def tips_of(run):
    return {t["tip"]: t for t in read_tsv(os.path.join(RUNS, run, "results", "tree_tips.tsv"))}


def keyed(rows, tips_all, order):
    """pairs.tsv is keyed by the tips of the all-sequences run; map by sequence."""
    by_seq = {rows[t]["aa_sequence"].strip().upper(): t for t in tips_all}
    return [by_seq[rows[o["tip"]]["aa_sequence"].strip().upper()] for o in order]


def colour_matrix(m):
    lo = np.array([b[0] for b in hm.BINS])
    rgba = np.array([to_rgba(b[1]) for b in hm.BINS])
    which = np.clip(np.searchsorted(lo, m, side="right") - 1, 0, len(lo) - 1)
    img = rgba[which]
    img[m == 0] = to_rgba(NO_ALN)
    d = np.arange(len(m))
    img[d, d] = to_rgba("#ffffff")        # a sequence against itself
    return img


def runs_of(keys):
    out, s = [], 0
    for i in range(1, len(keys) + 1):
        if i == len(keys) or keys[i] != keys[s]:
            out.append((keys[s], s, i))
            s = i
    return out


def save(fig, stem):
    for ext in ("svg", "pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"), dpi=DPI, facecolor=SURFACE,
                    **({"metadata": {"Date": None}} if ext == "svg" else {}))
    plt.close(fig)


def identity_legend(ax, x, y, w=3.4, bars=None, n_pairs=None):
    """Swatches with tier brackets, the not-alignable swatch, and the bars."""
    sw = w / 9
    ax.text(x, y + 0.42, "Pairwise amino acid identity (%)", fontsize=8.5, color=INK,
            fontweight="bold", va="bottom")
    if bars is not None:
        top = max(bars)
        for i, (b, v) in enumerate(zip(hm.BINS, bars)):
            h = 0.55 * v / top
            ax.add_patch(Rectangle((x + i * sw + 0.02, y + 0.95), sw - 0.04, h, fc=b[1], ec="none"))
            ax.text(x + i * sw + sw / 2, y + 0.97 + h, f"{100 * v:.0f}" if v >= 0.005 else "<1",
                    fontsize=5.6, ha="center", va="bottom", color=MUTED)
        ax.text(x, y + 1.75, f"share of the {n_pairs:,} alignable pairs, per bin (%)",
                fontsize=6, color=MUTED, va="bottom")
    for i, (lo, c) in enumerate(hm.BINS):
        ax.add_patch(Rectangle((x + i * sw + 0.02, y), sw - 0.04, 0.3, fc=c, ec="none"))
        ax.text(x + i * sw, y - 0.05, "<30" if lo == 20 else str(lo), fontsize=6, ha="center",
                va="top", color=INK)
    ax.text(x + 8 * sw, y - 0.05, "100", fontsize=6, ha="center", va="top", color=INK)
    for lab, a, b in TIERS:
        i0 = [lo for lo, _ in hm.BINS].index(a)
        i1 = 8 if b == 100 else [lo for lo, _ in hm.BINS].index(b)
        x0, x1 = x + i0 * sw + 0.03, x + i1 * sw - 0.03
        ax.plot([x0, x1], [y - 0.28, y - 0.28], color=MUTED, lw=0.7)
        ax.text((x0 + x1) / 2, y - 0.33, lab, fontsize=6.2, ha="center", va="top", color=MUTED)
    ax.add_patch(Rectangle((x + 0.02, y - 0.95), sw - 0.04, 0.3, fc=NO_ALN, ec="none"))
    ax.text(x + sw + 0.05, y - 0.8, "not alignable: the local alignment covers\nless than "
            "half of the shorter sequence", fontsize=6.2, va="center", color=INK)


def how_legend(ax, x, y):
    """Pair -> alignment -> coverage check -> identity -> colour."""
    ax.text(x, y, "How a cell gets its colour", fontsize=8.5, color=INK, fontweight="bold",
            va="bottom")
    steps = [("1  one pair", "row sequence i,\ncolumn sequence j"),
             ("2  align", "Smith-Waterman, BLOSUM62,\ngap 11/1 (BLAST defaults)"),
             ("3  coverage", "aligned columns /\nlength of the shorter"),
             ("4  identity", "identical columns /\naligned columns")]
    bw, bh, gap = 1.45, 0.62, 0.22
    for i, (t, s) in enumerate(steps):
        bx = x + i * (bw + gap)
        ax.add_patch(Rectangle((bx, y - 0.1 - bh), bw, bh, fc="white", ec=MUTED, lw=0.6))
        ax.text(bx + 0.07, y - 0.17, t, fontsize=6.6, va="top", color=INK, fontweight="bold")
        ax.text(bx + 0.07, y - 0.36, s, fontsize=5.6, va="top", color=MUTED, linespacing=1.25)
        if i:
            ax.add_patch(FancyArrowPatch((bx - gap + 0.02, y - 0.1 - bh / 2), (bx - 0.02, y - 0.1 - bh / 2),
                                         arrowstyle="-|>", mutation_scale=6, color=MUTED, lw=0.6))
    yy = y - 0.1 - bh - 0.18
    x3 = x + 2 * (bw + gap) + bw / 2
    ax.add_patch(FancyArrowPatch((x3, y - 0.1 - bh), (x3, yy - 0.12), arrowstyle="-|>",
                                 mutation_scale=6, color=MUTED, lw=0.6))
    ax.add_patch(Rectangle((x3 - 0.75, yy - 0.42), 0.16, 0.16, fc=NO_ALN, ec="none"))
    ax.text(x3 - 0.55, yy - 0.34, "under 50%: grey, no identity", fontsize=5.8, va="center",
            color=INK)
    x4 = x + 3 * (bw + gap) + bw / 2
    ax.add_patch(FancyArrowPatch((x4, y - 0.1 - bh), (x4, yy - 0.12), arrowstyle="-|>",
                                 mutation_scale=6, color=MUTED, lw=0.6))
    for i, (lo, c) in enumerate(hm.BINS[1:6]):
        ax.add_patch(Rectangle((x4 - 0.4 + i * 0.16, yy - 0.42), 0.16, 0.16, fc=c, ec="none"))
    ax.text(x4 + 0.45, yy - 0.34, "50% or more: binned\nidentity, scale above", fontsize=5.8,
            va="center", color=INK)


# --- main heatmap --------------------------------------------------------------

def main_heatmap(rows, zooms):
    order = ring_order("centroids")
    tips_all = tips_of("all-sequences")
    keys = keyed(rows, tips_all, order)
    n = len(keys)
    m, ix = hm.load_identity(keys)
    img = colour_matrix(m)
    img[np.triu_indices(n, 1)] = (0, 0, 0, 0)

    side = 10.0
    c = side / n
    W, H = side + 9.5, side + 2.4
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-4.6, W - 4.6)
    ax.set_ylim(-side - 1.4, 1.0)
    ax.axis("off")
    ax.set_aspect("equal")
    fig.patch.set_facecolor(SURFACE)
    ax.imshow(img, extent=(0, side, -side, 0), interpolation="nearest", zorder=1)

    fk = [o["fold_key"] for o in order]
    comp = [o["component"] for o in order]
    hue = [o["hue_hex"] for o in order]
    # margins: fold, component, on the left and along the bottom
    strips = [("fold", 0.16, lambda i: folds.band(fk[i])),
              ("component", 0.10, None)]
    xo = -0.06
    cr = runs_of(comp)
    alt = {cc: ("#8a8880", "#c9c7bf")[j % 2] for j, (cc, _, _) in enumerate(cr)}
    for name, w, colfn in strips:
        xo -= w
        for i in range(n):
            col = colfn(i) if colfn else alt[comp[i]]
            ax.add_patch(Rectangle((xo, -(i + 1) * c), w, c, fc=col, ec="none", lw=0, zorder=2))
            ax.add_patch(Rectangle((i * c, -side - 0.06 - (-xo - 0.06)), c, w, fc=col, ec="none",
                                   lw=0, zorder=2))
        ax.text(xo + w / 2, 0.08, name, rotation=90, fontsize=5.5, ha="center", va="bottom",
                color=MUTED)
        xo -= 0.03
    x_lab = xo - 0.08

    # fold outlines on the diagonal, labels in the fold colour
    labels = []
    for k, s, e in runs_of(fk):
        x0, x1, y0, y1 = s * c, e * c, -s * c, -e * c
        ax.add_patch(Polygon([(x0, y0), (x0, y1), (x1, y1)], closed=True, fill=False,
                             ec=folds.strong(k), lw=1.1, zorder=4))
        nseq = e - s
        nrec = sum(int(order[i]["size"]) for i in range(s, e))
        labels.append(((y0 + y1) / 2, f"{folds.name(k)}  ({nseq} seq, {nrec} rec)", folds.strong(k)))
    ys = hm_spread([y for y, _, _ in labels], 0.17)
    for (y, txt, col), ly in zip(labels, ys):
        ax.plot([x_lab + 0.02, x_lab - 0.12, x_lab - 0.22], [y, y, ly], color=col, lw=0.5)
        ax.text(x_lab - 0.26, ly, txt, fontsize=6.6, ha="right", va="center", color=col,
                fontweight="bold")

    # components of 5 or more, named just above the diagonal
    for cc, s, e in cr:
        if e - s >= 5:
            ax.text((s + e) / 2 * c + 0.06, -(s + e) / 2 * c + 0.06, cc, fontsize=5.4,
                    color=MUTED, rotation=-45, ha="left", va="bottom", rotation_mode="anchor")
    # seeds and the zoom boxes
    seed = {}
    for nm, pid in seed_nodes(rows).items():
        for i, o in enumerate(order):
            if pid in tips_of("centroids")[o["tip"]]["records"].split(";"):
                seed[i] = nm
    ss = sorted(seed)
    ly = hm_spread([-(i + 0.5) * c for i in ss], 0.15)
    for i, y in zip(ss, ly):
        x = (i + 0.5) * c
        ax.plot([x, x + 0.5, x + 0.62], [-(i + 0.5) * c, y, y], color=INK, lw=0.45, zorder=6)
        ax.scatter([x], [-(i + 0.5) * c], s=6, c=INK, zorder=7)
        ax.text(x + 0.66, y, seed[i], fontsize=6.4, va="center", color=INK, zorder=7)
    pos = {o["tip"]: i for i, o in enumerate(order)}
    for tag, (title, members) in zooms.items():
        idx = sorted(pos[t] for t in members if t in pos)
        if not idx:
            continue
        s, e = idx[0], idx[-1] + 1
        ax.add_patch(Rectangle((s * c, -e * c), (e - s) * c, (e - s) * c, fill=False, ec=INK,
                               lw=0.9, ls=(0, (3, 1.5)), zorder=8))
        ax.text(s * c - 0.05, -e * c + 0.02, f"zoom: {title}", fontsize=6.4, ha="right",
                va="bottom", color=INK, fontweight="bold", zorder=8)

    # legends
    iu = np.tril_indices(n, -1)
    vals = m[iu]
    al = vals[vals > 0]
    lo = np.array([b[0] for b in hm.BINS])
    which = np.clip(np.searchsorted(lo, al, side="right") - 1, 0, len(lo) - 1)
    bars = np.bincount(which, minlength=len(lo)) / len(al)
    lx = side * 0.52
    identity_legend(ax, lx, -1.6, bars=bars, n_pairs=len(al))
    how_legend(ax, lx, -3.4)
    ax.text(lx + 1.9, -5.15, f"{n} centroids, {len(vals):,} pairs; {100 * len(al) / len(vals):.0f}% "
            f"alignable, {100 - 100 * len(al) / len(vals):.0f}% not.\nOutlined triangles: one fold "
            "each, in the fold colour used by the tree.\nRows run in the order of the circular tree.", fontsize=6.4, color=MUTED,
            va="top", linespacing=1.4)
    save(fig, "figure_heatmap_v2")
    same = np.array(fk)[:, None] == np.array(fk)[None, :]
    lower = np.tril(np.ones((n, n), bool), -1)
    return {"pairs": int(lower.sum()), "alignable": int((m[lower] > 0).sum()),
            "not alignable, same fold": int(((m == 0) & same & lower).sum()),
            "not alignable, different fold": int(((m == 0) & ~same & lower).sum())}


def hm_spread(want, h):
    """Non overlapping label positions, top to bottom, as close as possible."""
    if not want:
        return []
    o = sorted(range(len(want)), key=lambda i: -want[i])
    w = [want[i] for i in o]
    v = [-(y) - k * h for k, y in enumerate(w)]
    blocks = []
    for x in v:
        blocks.append([x, 1])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            s_, n_ = blocks.pop()
            blocks[-1][0] += s_
            blocks[-1][1] += n_
    fit = [b[0] / b[1] for b in blocks for _ in range(b[1])]
    res = [0.0] * len(want)
    for k, i in enumerate(o):
        res[i] = -(fit[k] + k * h)
    return res


# --- zooms ---------------------------------------------------------------------

def find_zooms(rows):
    """The three groups, from all distinct sequences."""
    order = ring_order("all-sequences")
    tips = tips_of("all-sequences")
    keys = [o["tip"] for o in order]
    m, ix = hm.load_identity(keys)
    bycomp = collections.defaultdict(list)
    for o in order:
        bycomp[o["component"]].append(o["tip"])
    nrec = {t: len(tips[t]["records"].split(";")) for t in tips}
    sub50 = [g for c, mem in bycomp.items() if len(mem) >= 6
             for g in hm.subfamilies(mem, m, ix, 50)]
    sub70 = [g for c, mem in bycomp.items() if len(mem) >= 6
             for g in hm.subfamilies(mem, m, ix, 70)]
    isp = next(t for t in tips if tips[t]["enzyme_name"] == "ISPETase")
    g_isp = next(g for g in sub50 if isp in g)
    g_studied = max(sub50, key=lambda g: (sum(nrec[t] for t in g), len(g)))

    def med(g):
        idx = [ix[t] for t in g]
        return float(np.median(m[np.ix_(idx, idx)][np.triu_indices(len(g), 1)]))
    g_dense = max((g for g in sub70 if len(g) >= 20), key=med)
    return {"ispetase": ("IsPETase subfamily", g_isp),
            "studied": ("most studied subfamily", g_studied),
            "dense": ("densest subfamily", g_dense)}, m, ix, order, tips, nrec, med


def zoom_heatmap(rows, tag, title, members, m, ix, order, tips, nrec, med):
    pos = {o["tip"]: i for i, o in enumerate(order)}
    members = sorted(members, key=lambda t: pos[t])
    idx = [ix[t] for t in members]
    sub = m[np.ix_(idx, idx)]
    n = len(members)
    img = colour_matrix(sub)
    img[np.triu_indices(n, 1)] = (0, 0, 0, 0)
    c = min(0.16, 7.5 / n)
    side = n * c
    seeds = set(seed_nodes(rows).values())
    fig = plt.figure(figsize=(side + 6.5, side + 2.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-2.7, side + 3.8)
    ax.set_ylim(-side - 1.2, 1.0)
    ax.axis("off")
    ax.set_aspect("equal")
    fig.patch.set_facecolor(SURFACE)
    ax.imshow(img, extent=(0, side, -side, 0), interpolation="nearest")
    hue = {o["tip"]: o["hue_hex"] for o in order}
    fs = max(4.0, min(7.0, c * 46))
    numbers = n <= 40
    for i, t in enumerate(members):
        nm = tips[t]["enzyme_name"]
        is_seed = any(r in seeds for r in tips[t]["records"].split(";"))
        lab = f"{nm} ({nrec[t]})" if nrec[t] > 1 else nm
        ax.text(-0.2, -(i + 0.5) * c, lab, fontsize=fs, ha="right", va="center",
                color=INK, fontweight="bold" if is_seed or nm == "ISPETase" else "normal")
        if numbers:
            for j in range(i):
                v = sub[i, j]
                if v > 0:
                    ax.text((j + 0.5) * c, -(i + 0.5) * c, f"{v:.0f}", fontsize=fs * 0.62,
                            ha="center", va="center",
                            color="white" if v >= hm.LIGHT_TEXT_FROM else INK)
    k = order[pos[members[0]]]["fold_key"]
    ax.text(0, 0.62, f"Zoom: {title}", fontsize=10, fontweight="bold", color=INK)
    comps = sorted({order[pos[t]]["component"] for t in members})
    ax.text(0, 0.3, f"{n} distinct sequences, {sum(nrec[t] for t in members)} records, "
            f"component {', '.join(comps)}, {folds.name(k)}; median identity {med(members):.1f}%. "
            "Bold: curated seed. (n): records for that sequence.", fontsize=6.6, color=MUTED)
    identity_legend(ax, side + 0.3, -side + 2.6, w=3.2)
    save(fig, f"figure_zoom_{tag}")
    return n


# --- cluster sizes -------------------------------------------------------------

def cluster_sizes(rows):
    tips = tips_of("all-sequences")
    comp_n = collections.Counter()
    comp_rec = collections.Counter()
    comp_key = {}
    comp_names = collections.defaultdict(list)
    seeds = {pid: nm for nm, pid in seed_nodes(rows).items()}
    for t in tips.values():
        recs = t["records"].split(";")
        comp_n[t["component"]] += 1
        comp_rec[t["component"]] += len(recs)
        comp_key[t["component"]] = folds.fold_key(t)
        for r in recs:
            if r in seeds:
                comp_names[t["component"]].append(seeds[r])
    comps = sorted(comp_n, key=lambda c: (-comp_rec[c], c))
    members = collections.Counter()
    cname = {}
    for pid, r in rows.items():
        cid = (r["centroid_id"] or "").strip()
        if cid and cid.lower() not in ("nan", "na"):
            members[str(int(float(cid)))] += 1
    ct = tips_of("centroids")
    csize = {t: 1 + members.get(t, 0) for t in ct}
    top = sorted(csize, key=lambda t: (-csize[t], int(t)))[:25]

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(13, 7.6), gridspec_kw={"hspace": 0.55})
    fig.patch.set_facecolor(SURFACE)
    x = np.arange(len(comps))
    a1.bar(x, [comp_rec[c] for c in comps], color=[folds.band(comp_key[c]) for c in comps],
           edgecolor=[folds.strong(comp_key[c]) for c in comps], linewidth=0.5, width=0.82)
    a1.bar(x, [comp_n[c] for c in comps], color="none", edgecolor=INK, linewidth=0.0,
           width=0.82, hatch="")
    for i, c in enumerate(comps):
        a1.text(i, comp_rec[c] + 4, f"{comp_rec[c]}", ha="center", fontsize=5.6, color=MUTED)
        if i < 14:
            nm = ", ".join(comp_names[c][:3]) if comp_names[c] else \
                next(t["enzyme_name"] for t in tips.values() if t["component"] == c)
            a1.text(i, comp_rec[c] + 22 if i == 0 else comp_rec[c] + 14, nm, rotation=50,
                    fontsize=6, ha="left", va="bottom", color=INK)
    a1.set_xticks(x)
    a1.set_xticklabels([c.replace("c", "") for c in comps], fontsize=5.6, rotation=90)
    a1.set_ylabel("records", fontsize=8)
    a1.set_xlabel("SSN component (30% identity family), ranked by records", fontsize=8)
    a1.set_title("Family sizes: all 56 SSN components, coloured by fold", fontsize=10,
                 loc="left", color=INK)
    a1.set_ylim(0, max(comp_rec.values()) * 1.6)
    present = []
    for c in comps:
        if comp_key[c] not in present:
            present.append(comp_key[c])
    h = [Rectangle((0, 0), 1, 1, fc=folds.band(k), ec=folds.strong(k), lw=0.5) for k in present]
    a1.legend(h, [folds.name(k) for k in present], fontsize=6, ncol=2, frameon=False,
              loc="upper right", title="fold", title_fontsize=7)
    x2 = np.arange(len(top))
    a2.bar(x2, [csize[t] for t in top], color=[folds.band(folds.fold_key(ct[t])) for t in top],
           edgecolor=[folds.strong(folds.fold_key(ct[t])) for t in top], linewidth=0.5, width=0.8)
    a2.set_xticks(x2)
    a2.set_xticklabels([ct[t]["enzyme_name"] for t in top], rotation=55, ha="right", fontsize=6.5)
    for i, t in enumerate(top):
        a2.text(i, csize[t] + 0.3, str(csize[t]), ha="center", fontsize=6, color=MUTED)
    a2.set_ylabel("records", fontsize=8)
    a2.set_title("The 25 largest 90% clusters, named by their centroid", fontsize=10, loc="left",
                 color=INK)
    dist = collections.Counter(csize.values())
    a2.text(0.99, 0.95, "cluster size: " + ", ".join(f"{k}: {dist[k]}" for k in sorted(dist)),
            transform=a2.transAxes, ha="right", va="top", fontsize=6.5, color=MUTED)
    for a in (a1, a2):
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(labelsize=6.5)
        a.set_facecolor(SURFACE)
    save(fig, "figure_cluster_sizes")
    write_tsv(os.path.join(OUT, "component_sizes.tsv"),
              ["component", "fold", "distinct_sequences", "records", "seeds"],
              [[c, folds.name(comp_key[c]), comp_n[c], comp_rec[c], ";".join(comp_names[c])]
               for c in comps])
    return comps[:3], [(ct[t]["enzyme_name"], csize[t]) for t in top[:5]]


# --- lengths -------------------------------------------------------------------

def read_lengths():
    out = {}
    p = LENGTH_SOURCES["PETadex-nr families"]
    if os.path.exists(p):
        import pandas as pd
        d = pd.read_csv(p, usecols=["origin", "length", "dom_start", "dom_end"])
        d = d[d.origin == "NR"]
        out["PETadex-nr families, full ORF"] = d.length.to_numpy()
        out["PETadex-nr families, catalytic domain only"] = \
            (d.dom_end - d.dom_start + 1).dropna().astype(int).to_numpy()
    p = LENGTH_SOURCES["Logan, halo assay set"]
    if os.path.exists(p):
        seqs, cur = [], []
        for line in open(p):
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
        out["Logan, halo assay set"] = np.array([len(s.replace("-", "").replace("*", ""))
                                                 for s in seqs])
    return out


def lengths(rows):
    tips = tips_of("all-sequences")
    plast = np.array([len(rows[t]["aa_sequence"].strip()) for t in tips])
    by_fold = collections.defaultdict(list)
    for t in tips.values():
        by_fold[folds.fold_key(t)].append(len(rows[t["tip"]]["aa_sequence"].strip()))
    ext = read_lengths()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
    fig.patch.set_facecolor(SURFACE)
    bins = np.arange(0, 1001, 20)
    series = [("Plasticome (492 distinct)", plast, "#275ca0")] + \
        [(k, v, c) for (k, v), c in zip(ext.items(), ["#6f6f6f", "#bdbdbd", "#a3842a"])]
    stats = []
    for lab, v, col in series:
        a1.hist(np.clip(v, 0, 1000), bins=bins, density=True, histtype="step", lw=1.4,
                color=col, label=f"{lab}: n={len(v):,}, median {np.median(v):.0f} aa")
        a1.axvline(np.median(v), color=col, lw=0.8, ls=(0, (3, 2)))
        stats.append([lab, len(v), int(np.median(v)), int(np.percentile(v, 5)),
                      int(np.percentile(v, 95)), int(v.min()), int(v.max())])
    a1.set_xlabel("length (aa), values over 1000 piled at 1000", fontsize=8)
    a1.set_ylabel("density", fontsize=8)
    a1.legend(fontsize=6.6, frameon=False)
    a1.set_title("Sequence length: plasticome against PETadex-nr and Logan", fontsize=9.5,
                 loc="left")
    ks = sorted(by_fold, key=lambda k: -len(by_fold[k]))
    data = [by_fold[k] for k in ks]
    bp = a2.boxplot(data, orientation="horizontal", patch_artist=True, widths=0.6, showfliers=True,
                    flierprops=dict(markersize=2, markeredgewidth=0.3))
    for p, k in zip(bp["boxes"], ks):
        p.set_facecolor(folds.band(k))
        p.set_edgecolor(folds.strong(k))
    for med_ in bp["medians"]:
        med_.set_color(INK)
    a2.set_yticks(range(1, len(ks) + 1))
    a2.set_yticklabels([f"{folds.name(k)} ({len(by_fold[k])})" for k in ks], fontsize=6.3)
    a2.invert_yaxis()
    a2.set_xlabel("length (aa)", fontsize=8)
    a2.set_title("Plasticome length by fold", fontsize=9.5, loc="left")
    for a in (a1, a2):
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(labelsize=6.5)
        a.set_facecolor(SURFACE)
    fig.tight_layout()
    save(fig, "figure_lengths")
    write_tsv(os.path.join(OUT, "length_summary.tsv"),
              ["set", "n", "median", "p5", "p95", "min", "max"], stats)
    return stats


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = load_rows()
    zooms, m, ix, order, tips, nrec, med = find_zooms(rows)
    # the zoom groups are named by all-sequences tips; the centroid figure
    # marks the centroids among them
    cent_by_seq = {rows[t]["aa_sequence"].strip().upper(): t for t in tips_of("centroids")}
    zc = {k: (title, [cent_by_seq.get(rows[t]["aa_sequence"].strip().upper()) for t in g])
          for k, (title, g) in zooms.items()}
    counts = main_heatmap(rows, zc)
    print("main:", counts)
    for tag, (title, g) in zooms.items():
        n = zoom_heatmap(rows, tag, title, g, m, ix, order, tips, nrec, med)
        print(f"zoom {tag}: {n} sequences, {sum(nrec[t] for t in g)} records, median {med(g):.1f}")
    write_tsv(os.path.join(OUT, "zoom_groups.tsv"), ["zoom", "tip", "enzyme_name", "records"],
              [[k, t, tips[t]["enzyme_name"], tips[t]["records"]] for k, (_, g) in zooms.items()
               for t in g])
    print("largest:", cluster_sizes(rows))
    for s in lengths(rows):
        print("length:", s)


if __name__ == "__main__":
    main()
