#!/usr/bin/env python3
"""
Step 18: the forest as one tree with a ring of plastics and the catalytic
mechanism of every fold beside it (Figure 2a).

Reviewer feedback of 2026-10-28 (Figure 2a), applied on top of the list below:
folds are grouped and named by type (folds.TYPES), the fold, accessory domain
and catalytic residue rings are gone (catalytic_conservation.tsv is still
written), every tip is a plain circle (no organism shapes), the structures are
larger and sit in a grid matched to the slices by numbered badges, and
--vertical draws the same tree rectangular, root on the left. The sections
below on the shapes, the mechanism layer and the accessory domains describe the
version of 2026-10-27.

Built for the reviewer feedback of 2026-10-27 on the pie (step 8 and the hand
edited copy). What changed against the pie, item by item:

  every fold shown       the six single sequence folds are slices of their
                         own, each with its own colour (folds.py), not one
                         merged "Other folds" slice
  full dataset           default is the all-sequences run: 492 tips standing
                         for all 609 records. A tip's area is the number of
                         records it stands for. With PLASTICOME_TREE_SET=
                         centroids the tips are the 411 centroids and the area
                         is the size of each 90 percent cluster instead
  ordering               folds clockwise from 12 o'clock, largest first; inside
                         a tree every node lists its larger child first, so
                         the components of a fold run largest to smallest.
                         The order is written to ring_order.tsv, which the
                         heatmap reads, so the two figures agree
  substrate rings        one ring per plastic, filled where any record of the
                         tip was reported on that plastic. Ring names sit in
                         the gap at 12 o'clock
  organism               tip shape: diamond bacterial, triangle archaeal,
                         hexagon eukaryotic, circle metagenomic
  bootstrap values       printed at the nodes. Inside a component the number is
                         the ensemble bootstrap of step 2 (the share of 100
                         FastTree trees, each on columns drawn from all 16
                         Muscle5 alignments, that contain the branch), shown at
                         70 or more. On the dashed backbone it is the UFBoot of
                         step 5 (1000 ultrafast bootstraps, IQ-TREE on the
                         FoldMason amino acid + 3Di alignment)
  line types             four kinds, each in the legend: dashed backbone (not
                         to scale), solid branch (to scale), pale branch
                         (support under 70), dotted leader to a structure, and
                         the faint coloured guides from each tip to the rings
  faint rainbow          every tip has a hue from its position around the
                         circle. Its guide to the rings and the outermost
                         "position" ring carry that hue, and the heatmap axes
                         use the same hues, so a sequence can be followed from
                         tip to rings to heatmap
  structures             the step 17 renders, larger, with the catalytic
                         residues, the modelled plastic and the cofactors
  mechanism layer        a ring for whether the tip keeps the catalytic
                         residues of its fold's annotated structure, from a
                         pairwise alignment to that structure's sequence. Where
                         the pairwise alignment is too weak to say, the cell is
                         marked as such; that gap is what the HMM layer is for
  accessory domains      a ring for the accessory CATH domain beside the
                         catalytic one

Writes runs/<set>/results/figures/<design>/tree-circle.{svg,pdf,png,html} (tree-vertical.*
with --vertical), ring_order.tsv and catalytic_conservation.tsv. The designs are plain,
fold-groups (--type-bars) and leaders (--leaders); the .html is the same figure as an editor
(ring_scene.py). Drag the structures and legends, double click a label to
reword it, switch layers off, then Export layout; --layout reads that file
back so the edits carry into the PDF.

The structures go in a grid left of the ring, in clockwise order, and a
numbered badge in the fold's colour on the slice and on the picture joins them
(tree-circle.*). --leaders gives the earlier circle with dotted leaders
(figures/leaders/); --vertical the rectangular tree (tree-vertical.*).
Each has its own layout file.

Usage: python 18_ring_tree.py [--rainbow] [--no-structures] [--vertical] [--leaders]
                              [--layout ring_layout.json]
"""

import argparse
import collections
import colorsys
import importlib
import json
import math
import os
import re
import textwrap

import numpy as np
from Bio import Align, Phylo
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib.patches import PathPatch, Rectangle, Wedge
from matplotlib.path import Path as MPath
from matplotlib.textpath import TextPath

import catalytic_struct
import folds
import ring_scene
from common import (INK, MUTED, NODE_STROKE, OUT, SET, SUBSTRATE_COLOR, SURFACE, WORK,
                    load_rows, read_tsv, seed_nodes, write_tsv)

d7 = importlib.import_module("07_draw")
d8 = importlib.import_module("08_draw_pie")
cards = importlib.import_module("17_mechanism_cards")
grey = importlib.import_module("17b_mechanism_colour_renders")
walk, STYLE, BACKBONE_STEP = d7.walk, d7.STYLE, d7.BACKBONE_STEP
plt = d7.plt

HERE = os.path.dirname(os.path.abspath(__file__))

INNER, RIM = 1.0, 5.0           # inches: root circle and deepest tip
GAP_ROWS = 3.0                  # between two folds
TOP_GAP_ROWS = 22.0             # the opening at 12 o'clock for the ring names
R_RING = RIM + 0.22             # first ring
RING_W, RING_GAP = 0.105, 0.025
STRUCT_R, STRUCT_IN = 9.9, 1.8    # structure centres and their size, inches
STAGGER = 2.55                    # inches: the outer row of structures sits this much further out
CAPTION_H = 0.36                  # inches below a picture: name and two lines of caption
PANEL_COLS, BADGE = 4, 9          # structure grid columns, and badge diameter in points
PANEL_IN = 2.8                    # inches: a structure in the grid (1.8 with --leaders)
ROW_H = 0.03                      # --vertical: inches per record (a row of the tree)
VPANEL_COLS = 3                   # --vertical: grid columns
# the bars that name the fold types: hydrolases slate, oxidoreductases plum
TYPE_COLORS = {"activity": ["#2b3648", "#4f6178", "#6d7f96", "#5a3a6e", "#85629b"],
               "mechanism": ["#2b3648", "#4f6178", "#6d7f96", "#8796a8", "#4a6a62", "#5a3a6e"]}

PLASTICS = ["PET", "PBAT", "PLA", "PUR", "PBS", "PHA", "PEF", "PC", "PA", "PE", "NR"]
RING_COLOR = dict(SUBSTRATE_COLOR, PEF="#7a1f2b")
MONO_PLASTIC = "#555555"        # --mono-plastics: one grey for every plastic column
ACCESSORY = [("Trefoil", "trefoil", "#e69f00"),
             ("Ig-like", "immunoglobulin", "#56b4e9"),
             ("Galactose-binding", "galactose", "#009e73"),
             ("Serralysin C-term", "serralysin", "#cc79a7"),
             ("other accessory", "", "#7f7f7f")]
CONSERVED = [("all catalytic residues kept", "#2b2b2b"),
             ("some kept", "#8f8f8f"),
             ("none kept", "#d8d8d4"),
             ("not assessed", "none"),
             ("no annotated structure", "none")]
SHAPES = [("Bacterial", "D"), ("Archaeal", "^"), ("Eukaryotic", "h"), ("Metagenomic", "o")]
LINES = {"backbone": (INK, 0.9), "solid": ("#555555", 0.55), "weak": ("#c4c4c0", 0.55)}


def shape(rows, recs):
    r = rows[recs[0]]
    if (r["metagenomic"] or "").strip().lower() == "yes":
        return "o"
    return {"Bacteria": "D", "Archaea": "^", "Eukaryota": "h"}.get(
        (r["superkingdom"] or "").strip(), "o")


def plastics_of(rows, recs):
    out = set()
    for p in recs:
        for s in (rows[p]["substrate"] or "").split(";"):
            s = s.strip()
            if s:
                out.add(s)
    return out


def accessory_of(rows, recs):
    """The accessory domain category of a tip, or None."""
    doms = set()
    for p in recs:
        parts = [x.strip() for x in (rows[p]["domain"] or "").split(";")[1:] if x.strip()]
        doms |= set(parts)
    if not doms:
        return None
    for lab, key, _ in ACCESSORY:
        if key and any(key in d.lower() for d in doms):
            return lab
    return "other accessory"


# --- catalytic residues kept? --------------------------------------------------

ALIGNER = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1,
                                end_gap_score=0)
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")


def reference(a):
    """Sequence of the annotated model (all protein chains joined) and the
    index of every catalytic residue in it."""
    model = next(iter(PDBParser(QUIET=True).get_structure("r", a["pdb"])))
    seq, index = [], {}
    for ch in model:
        for r in ch:
            if r.id[0] == " " and "CA" in r:
                index[(ch.id, r.id[1])] = len(seq)
                seq.append(seq1(r.get_resname()))
    want = [(index[(ch, n)], aa) for _, res in a["primary"] for ch, aa, n, _ in res
            if (ch, n) in index]
    return "".join(seq), want


TOLERANCE = 2      # residues either side, the same slip allowed on the structure route
RULE = f"tol{TOLERANCE}-idshort"      # the cached table is recomputed when this changes


def _kept(seq, col, i, aa):
    """Whether the sequence has residue `aa` where the model's column i falls,
    or within TOLERANCE of it. A column the sequence has a gap at falls before
    the sequence's next aligned residue."""
    j = col.get(i, -1)
    if j < 0:
        j = next((col[c] for c in sorted(col) if c > i and col[c] >= 0), len(seq))
        near = range(j - TOLERANCE, j + TOLERANCE)
    else:
        near = range(j - TOLERANCE, j + TOLERANCE + 1)
    return any(0 <= x < len(seq) and seq[x] == aa for x in near)


def conservation(tips, rows, keyof):
    """tip -> (category, identity, kept, total); cached in OUT.

    By sequence first: a pairwise alignment to the fold's annotated model,
    a residue counting as kept when it is at its aligned position or within
    TOLERANCE of it, as on the structure route. A
    sequence under 25 percent identity to the model (identical residues over
    the shorter of the two) is then placed by
    structure (catalytic_struct.py: its component's representative structure
    superposed on the model). What neither can place is "not assessed"."""
    path = os.path.join(OUT, "catalytic_conservation.tsv")
    if os.path.exists(path):
        cached = read_tsv(path)
        if cached and cached[0].get("rule") == RULE:
            return {r["tip"]: (r["category"], float(r["identity"]), int(r["kept"]),
                               int(r["total"])) for r in cached}
    A = folds.annotations()
    refs = {k: reference(a) for k, a in A.items()}
    out, method, extra, weak = {}, {}, {}, []
    for t in tips:
        k = keyof(t)
        seq = rows[t["records"].split(";")[0]]["aa_sequence"].strip().upper()
        if k not in refs:
            out[t["tip"]] = ("no annotated structure", 0.0, 0, 0)
            method[t["tip"]] = "-"
            continue
        ref, want = refs[k]
        aln = ALIGNER.align(ref, seq)[0]
        idx = aln.indices
        col = {int(i): j for i, j in zip(idx[0], idx[1]) if i >= 0}
        pairs = [(i, j) for i, j in zip(idx[0], idx[1]) if i >= 0 and j >= 0]
        # identical residues over the shorter sequence: over the aligned pairs
        # alone, a sequence that overlaps the model by a few residues (free
        # end gaps) would score as identical
        ident = sum(ref[i] == seq[j] for i, j in pairs) / min(len(ref), len(seq))
        kept = sum(1 for i, aa in want if _kept(seq, col, i, aa))
        if ident < 0.25:
            rep_pdb = os.path.join(WORK, "structures", f"{t['component']}.pdb")
            sites = [(ch, n, aa) for _, res in A[k]["primary"] for ch, aa, n, _ in res]
            weak.append((t["tip"], seq, rep_pdb if os.path.exists(rep_pdb) else None,
                         A[k]["pdb"], sites))
            out[t["tip"]] = (None, ident, 0, len(want))
            continue
        cat = ("all catalytic residues kept" if kept == len(want) else
               "some kept" if kept else "none kept")
        out[t["tip"]] = (cat, ident, kept, len(want))
        method[t["tip"]] = "sequence"
    if weak:
        got = catalytic_struct.assess(weak)
        for tip, (cat, tm, rid, kept, total) in got.items():
            out[tip] = (cat, out[tip][1], kept, total)
            method[tip] = "structure" if cat != "not assessed" else "-"
            extra[tip] = (tm, rid)
    table = []
    for t in tips:
        c, i, kk, n = out[t["tip"]]
        tm, rid = extra.get(t["tip"], ("", ""))
        table.append([t["tip"], t["enzyme_name"], keyof(t), c, method[t["tip"]], f"{i:.3f}",
                      f"{tm:.3f}" if tm != "" else "", f"{rid:.3f}" if rid != "" else "",
                      kk, n, RULE])
    write_tsv(path, ["tip", "enzyme_name", "fold_key", "category", "method", "identity",
                     "tm_rep_to_model", "identity_to_rep", "kept", "total", "rule"], table)
    return out


# --- layout helpers -------------------------------------------------------------

def spread_angles(want, sep):
    """Angles as close to `want` as possible, at least `sep` apart, in order.
    Least squares by pool adjacent violators, as d8.spread does for heights."""
    order = sorted(range(len(want)), key=lambda i: want[i])
    w = [want[i] for i in order]
    v = [x - k * sep for k, x in enumerate(w)]
    blocks = []
    for x in v:
        blocks.append([x, 1])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            s, n = blocks.pop()
            blocks[-1][0] += s
            blocks[-1][1] += n
    fit = [b[0] / b[1] for b in blocks for _ in range(b[1])]
    res = [0.0] * len(want)
    for k, i in enumerate(order):
        res[i] = fit[k] + k * sep
    return res


LANE_IN, LANE_GAP = 0.32, 0.10      # inches: first leader lane beyond the rings, and the step


def place_structures(want, sep, lo=0.0, hi=2 * math.pi):
    """Clockwise angle of each structure: as near its slice as possible, at
    least `sep` apart, in slice order (so leaders can be routed without
    crossing), and inside [lo, hi] so none drifts into the opening at 12
    o'clock where the ring names are. A panel layout would replace this."""
    # no structure beyond the first or the last slice: past them a leader
    # would run back toward the opening and along the alpha/beta hydrolase
    lo, hi = max(lo, min(want)), min(hi, max(want))
    got = spread_angles(want, sep)
    order = sorted(range(len(got)), key=lambda i: got[i])
    x = [got[i] for i in order]
    for _ in range(2):                       # clamp each end, push the rest along
        for n in range(len(x) - 1, -1, -1):
            x[n] = min(x[n], hi if n == len(x) - 1 else x[n + 1] - sep)
        for n in range(len(x)):
            x[n] = max(x[n], lo if n == 0 else x[n - 1] + sep)
    out = [0.0] * len(got)
    for n, i in enumerate(order):
        out[i] = x[n]
    return out


def assign_lanes(start, end):
    """A lane per leader, so none cross. A leader runs out from its slice to
    its lane, along the lane to the structure's angle, then out to the
    structure. Two leaders whose lane runs overlap travel the same way (the
    order is kept), and the one further along gets the inner lane: then
    neither crosses the other. Lane 0 is the innermost."""
    n = len(start)
    lane = [0] * n
    order = sorted(range(n), key=lambda i: start[i])
    span = lambda i: (min(start[i], end[i]), max(start[i], end[i]))
    overlap = lambda i, j: span(i)[0] < span(j)[1] and span(j)[0] < span(i)[1]
    # clockwise movers: furthest along first; anticlockwise: the reverse
    for direction, seq in ((1, order[::-1]), (-1, order)):
        done = []
        for i in seq:
            if (end[i] - start[i]) * direction <= 1e-9:
                continue
            ahead = [j for j in done if overlap(i, j)]
            lane[i] = 1 + max((lane[j] for j in ahead), default=-1)
            done.append(i)
    return lane


def leader(a_mid, r0, r_lane, cx, cy, box, arc):
    """The leader path: radial from the ring to the lane, along the lane, then
    straight towards the structure, stopping at the edge of its box (picture
    and caption, (x0, y0, x1, y1)) so it never runs through its own text."""
    a_s = math.atan2(cy, cx)
    a_s = a_mid + (a_s - a_mid + math.pi) % (2 * math.pi) - math.pi     # short way round
    pts = [(r0 * math.cos(a_mid), r0 * math.sin(a_mid))] + arc(r_lane, a_mid, a_s)
    ex, ey = pts[-1]
    dx, dy = cx - ex, cy - ey
    # the first t in (0, 1] where the segment enters the box, then a small gap
    t = 1.0
    for lo_, hi_, p0, dp in ((box[0], box[2], ex, dx), (box[1], box[3], ey, dy)):
        if abs(dp) > 1e-9:
            t_in = min((lo_ - p0) / dp, (hi_ - p0) / dp)
            if 0 < t_in < t:
                other = (ey + dy * t_in) if p0 == ex else (ex + dx * t_in)
                lo2, hi2 = ((box[1], box[3]) if p0 == ex else (box[0], box[2]))
                if lo2 - 1e-6 <= other <= hi2 + 1e-6:
                    t = t_in
    t = max(t - 0.06 / (math.hypot(dx, dy) or 1), 0)
    if t > 0:
        pts.append((ex + dx * t, ey + dy * t))
    return pts


def struct_box(cx, cy, size_in, pictured):
    """The area a structure and its caption take, in inches."""
    if pictured:
        return (cx - size_in / 2, cy - size_in / 2 - CAPTION_H, cx + size_in / 2, cy + size_in / 2)
    return (cx - 0.9, cy - 0.12, cx + 0.9, cy + 0.12)


def _cross(p, q, r, s):
    d = lambda a, b, c: (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (d(p, q, r) * d(p, q, s) < 0) and (d(r, s, p) * d(r, s, q) < 0)


def _box(t):
    """A text's corners, in inches: width from the character count (an
    estimate; good enough to catch collisions, not to typeset)."""
    lines = str(t["s"]).split("\n")
    w = max(len(l) for l in lines) * t["size"] * 0.55 / 72
    h = len(lines) * t["size"] * 1.25 / 72
    x0 = {"left": 0, "right": -w, "center": -w / 2}[t["ha"]]
    y0 = {"top": -h, "bottom": 0, "baseline": 0, "center": -h / 2}[t["va"]]
    a = math.radians(t.get("rotation") or 0)
    c, s_ = math.cos(a), math.sin(a)
    return [(t["x"] + c * x - s_ * y, t["y"] + s_ * x + c * y)
            for x, y in ((x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h))]


def _overlap(A, B):
    """Separating axis test for two convex quadrilaterals."""
    for poly in (A, B):
        for i in range(4):
            (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % 4]
            nx, ny = y1 - y2, x2 - x1
            pa = [nx * x + ny * y for x, y in A]
            pb = [nx * x + ny * y for x, y in B]
            if max(pa) <= min(pb) or max(pb) <= min(pa):
                return False
    return True


def check_layout(sc):
    """Counts of what makes the figure look sloppy: leaders crossing leaders,
    text on text, text on a structure. The aim is zero of each."""
    out = collections.Counter()
    leads = [st["leader"] for st in sc.structures]
    for i in range(len(leads)):
        for j in range(i + 1, len(leads)):
            if any(_cross(p, q, r, s) for p, q in zip(leads[i], leads[i][1:])
                   for r, s in zip(leads[j], leads[j][1:])):
                out["leader crossings"] += 1
    sq_of = lambda st, f=0.45: [(st["x"] - st["size"] * f, st["y"] - st["size"] * f),
                                (st["x"] + st["size"] * f, st["y"] - st["size"] * f),
                                (st["x"] + st["size"] * f, st["y"] + st["size"] * f),
                                (st["x"] - st["size"] * f, st["y"] + st["size"] * f)]
    for i, st in enumerate(sc.structures):
        for j, other in enumerate(sc.structures):
            if i == j or not other["image"]:
                continue
            q = sq_of(other)
            edges = list(zip(q, q[1:] + q[:1]))
            if any(_cross(p, r, a, b) for p, r in zip(st["leader"], st["leader"][1:])
                   for a, b in edges):
                out["leaders through a structure"] += 1
    texts = [it for it in sc.items if it["kind"] == "text"]
    for st in sc.structures:
        ty = st["y"] - st["size"] * 0.5 - 0.02 if st["image"] else st["y"] + 0.05
        texts.append({"s": st["name"], "size": 6.6, "x": st["x"], "y": ty, "ha": "center",
                      "va": "top"})
        texts.append({"s": st["lines"], "size": 4.8, "x": st["x"], "y": ty - 0.14,
                      "ha": "center", "va": "top"})
    boxes = [_box(t) for t in texts]
    for st in sc.structures:
        for b in boxes:
            edges = list(zip(b, b[1:] + b[:1]))
            if any(_cross(p, r, a, c) for p, r in zip(st["leader"], st["leader"][1:])
                   for a, c in edges):
                out["leaders through text"] += 1
                if os.environ.get("RING_DEBUG"):
                    print("  leader of", st["name"], "crosses text", texts[boxes.index(b)]["s"][:40])
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if _overlap(boxes[i], boxes[j]):
                out["text overlaps"] += 1
    for st in sc.structures:
        if not st["image"]:
            continue
        h = st["size"] * 0.45
        sq = [(st["x"] - h, st["y"] - h), (st["x"] + h, st["y"] - h),
              (st["x"] + h, st["y"] + h), (st["x"] - h, st["y"] + h)]
        out["text on a structure"] += sum(_overlap(sq, b) for b in boxes)
    return {k: v for k, v in out.items() if v}


def rainbow(frac):
    r, g, b = colorsys.hsv_to_rgb(frac % 1.0, 0.55, 0.92)
    return (r, g, b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rainbow", action="store_true",
                    help="tracking hue: guide lines from tips to rings and a position ring "
                         "(off by default)")
    ap.add_argument("--no-structures", action="store_true")
    ap.add_argument("--vertical", action="store_true",
                    help="a rectangular tree, root on the left, tips stacked top to bottom "
                         "(writes tree-vertical.*)")
    ap.add_argument("--type-bars", action="store_true",
                    help="draw a bar and a name for each fold type (folds.TYPES); off since "
                         "2026-10-28, the folds are still ordered by type")
    ap.add_argument("--grouping", choices=sorted(folds.GROUPINGS), default="activity",
                    help="how --type-bars and the slice order group the folds: enzyme "
                         "activity (EC-style) or catalytic mechanism")
    ap.add_argument("--leaders", action="store_true",
                    help="the earlier circle: structures around the ring on dotted leaders "
                         "(writes figures/leaders/tree-circle.*)")
    ap.add_argument("--layout", default=None,
                    help="edits exported from the editor (default ring_layout.<set>.json, "
                         "or ring_layout.<set>.panel.json with --panel; used when it exists)")
    ap.add_argument("--pictures", default=None,
                    help="folder with full/<slug>.png to use for the structures instead of step 17's")
    ap.add_argument("--fill-cells", action="store_true",
                    help="vertical mechanism figure: each structure as wide as its grid cell "
                         "instead of at most 5.5 in")
    ap.add_argument("--fit-content", action="store_true",
                    help="vertical mechanism figure: each row runs from the column headers to "
                         "the scale bar, and each protein is scaled until its drawn extent fills "
                         "its column or its row; captions align along the row")
    ap.add_argument("--text-pt", type=float, default=None,
                    help="vertical mechanism figure: every text at this size (default 37.5)")
    ap.add_argument("--line-scale", type=float, default=1.0,
                    help="thicken the tree's branches and the scale bar by this factor (print)")
    ap.add_argument("--no-group-names", action="store_true",
                    help="vertical mechanism figure: only the numbered circles beside the tree, "
                         "the names stay under the structures")
    ap.add_argument("--min-inset-mag", type=float, default=0.0,
                    help="with --insets --fit-content: leave out a close up that magnifies less")
    ap.add_argument("--mono-plastics", action="store_true",
                    help="vertical mechanism figure: every plastic column (and header) in one grey, "
                         "so colour means only the mechanism group")
    ap.add_argument("--no-tips", action="store_true", help="leave out the tip dots")
    ap.add_argument("--label-residues", action="store_true",
                    help="vertical mechanism figure: name the catalytic residues on each protein, "
                         "with a thin line from each name to its residue")
    ap.add_argument("--residues", action="store_true",
                    help="vertical mechanism figure: the catalytic residues under each title, in "
                         "the group colour (from the structure renders' scenes.json)")
    ap.add_argument("--order-by-size", action="store_true",
                    help="groups from the most tips to the fewest, not in folds.TYPES order "
                         "(the numbers 1 to 6 follow)")
    ap.add_argument("--row-h", type=float, default=None,
                    help="--vertical: inches per tree row (default 0.03); smaller makes the figure shorter")
    ap.add_argument("--grid-share", type=float, default=0.6,
                    help="--fit-content: the structures take at most this share of the width")
    ap.add_argument("--no-ids", action="store_true",
                    help="vertical mechanism figure: no structure ID under each title (in the caption)")
    ap.add_argument("--font", default=None, help="font family for every text, e.g. Arial")
    ap.add_argument("--band", action="store_true",
                    help="vertical mechanism figure: draw the lavender band to the next step "
                         "(removed by hand 2026-11-14, off by default)")
    ap.add_argument("--insets", action="store_true",
                    help="vertical mechanism figure: a round close up of the active site in a "
                         "corner of each structure, joined to a ring around the site")
    ap.add_argument("--stem", default=None,
                    help="write the figure to this path stem instead of results/figures/...")
    args = ap.parse_args()
    V = args.vertical
    args.panel = not args.leaders or V         # a grid of structures with numbered badges
    # one folder per design, the circle and vertical trees of a design side by side
    folds.use_grouping(args.grouping)
    few = args.grouping == "mechanism"         # six pictures: two columns, larger
    if args.font:                              # e.g. Arial for submission; TrueType kept in the PDF
        plt.rcParams["font.family"] = args.font
        plt.rcParams["pdf.fonttype"] = 42
    if args.mono_plastics:
        for p_ in PLASTICS:
            RING_COLOR[p_] = MONO_PLASTIC
    if V and few and args.type_bars:
        # vertical mechanism figure: two branch styles only, backbone (dotted) and
        # the rest; the support split is left out (2026-11-14)
        LINES["weak"] = LINES["solid"]
    TS = 3.6 if few else 1.0             # text size, for a figure shrunk into a four panel composite
    VFS = (args.text_pt or 37.5) if V and few else None    # vertical mechanism figure: every text at 37.5 pt (30 until 2026-11-14)
    F = lambda size: VFS or size
    global ROW_H
    if args.row_h:
        ROW_H = args.row_h
    if few:                                    # a wide tree and wide rings, so the figure is square
        global INNER, RIM, R_RING, RING_W, BADGE
        BADGE *= TS
        RIM, RING_W = 7.0, 0.085 * TS         # a ring as wide as its name is tall
        R_RING = RIM + 0.22
        if VFS:                                # 30 pt text: a column fits its name, a badge its number
            RING_W = max(RING_W, 1.1 * VFS / 72)
            BADGE = max(BADGE, 1.45 * VFS)

    def text_w(s_, size, bold=False):
        """Width in inches of the widest line of s_ at size points."""
        fp_ = FontProperties(family=plt.rcParams["font.family"], weight="bold" if bold else "normal")
        return max((TextPath((0, 0), l_, size=size, prop=fp_).get_extents().width / 72
                    if l_.strip() else 0.0) for l_ in s_.split("\n"))
    TYPE_COLOR = TYPE_COLORS[args.grouping]
    if args.grouping == "mechanism":
        # one pictured fold per group, rendered in the group's colour by
        # PLASTICOME_GROUPING=mechanism 17_mechanism_cards.py
        cards.OUTDIR = os.path.join(HERE, "runs", "mechanisms-mechanism")
        if V:                                  # grey protein, mechanism in group colour (17b, 2026-11-14)
            cards.OUTDIR = grey.OUTDIR
    if args.pictures:                          # other renders of the same scenes (a colour demo)
        cards.OUTDIR = os.path.abspath(args.pictures)
    DESIGN = ("leaders" if args.leaders else "plain" if not args.type_bars
              else "fold-groups" if args.grouping == "activity" else "mechanism-groups")
    name_ = "tree-vertical" if V else "tree-circle"
    # mechanism groups: each layout in its own folder (2026-11-14)
    STEM = os.path.join("figures", DESIGN, *([name_] if DESIGN == "mechanism-groups" else []), name_)
    if args.stem:
        STEM = os.path.abspath(args.stem)
    os.makedirs(os.path.join(OUT, os.path.dirname(STEM)), exist_ok=True)
    args.layout = args.layout or os.path.join(
        HERE, f"ring_layout.{SET}{'.vertical' if V else '.leaders' if args.leaders else ''}.json")
    layout = {}
    if os.path.exists(args.layout):
        with open(args.layout) as fh:
            layout = json.load(fh)
        print(f"layout: {os.path.relpath(args.layout, HERE)}")
    retext = lambda tid, s: layout.get("text", {}).get(tid, s)
    sc = ring_scene.Scene()
    # matplotlib's three branch styles and the dash of the backbone, for the scene
    def line_style(st):
        col, lw = LINES[st] if isinstance(st, str) else st[:2]
        return col, lw, ([2.2 * lw, 1.6 * lw] if st == "backbone" else None)

    rows = load_rows()
    tips = read_tsv(os.path.join(OUT, "tree_tips.tsv"))
    tipinfo = {t["tip"]: t for t in tips}
    keyof = folds.fold_key
    record_tip = {r: t["tip"] for t in tips for r in t["records"].split(";")}
    seeds = {record_tip[pid]: name for name, pid in seed_nodes(rows).items()
             if pid in record_tip}

    # what the area of a tip stands for
    if SET == "all":
        size = {t["tip"]: len(t["records"].split(";")) for t in tips}
        size_label = "records per distinct sequence"
    else:
        members = collections.Counter()
        for pid, r in rows.items():
            cid = (r["centroid_id"] or "").strip()
            if cid and cid.lower() not in ("nan", "na"):
                members[str(int(float(cid)))] += 1
        size = {t["tip"]: 1 + members.get(t["tip"], 0) for t in tips}
        size_label = "records in the 90% cluster"

    groups = collections.defaultdict(list)
    for t in tips:
        groups[t["fold_group"]].append(t)
    gkey = {g: keyof(groups[g][0]) for g in groups}
    # by fold type, then the larger fold first inside a type. --order-by-size:
    # the types themselves from the most tips to the fewest (2026-11-14)
    type_tips = collections.Counter()
    for g_, ts_ in groups.items():
        type_tips[folds.type_rank(gkey[g_])] += len(ts_)
    type_pos = ({r_: i for i, r_ in enumerate(sorted(type_tips, key=lambda r_: (-type_tips[r_], r_)))}
                if args.order_by_size else {r_: r_ for r_ in type_tips})
    order = sorted(groups, key=lambda g: (type_pos[folds.type_rank(gkey[g])], -len(groups[g]), g))
    n_comp = {g: len({t["component"] for t in ts}) for g, ts in groups.items()}
    comp_size = collections.Counter(t["component"] for t in tips)

    trees = {}
    for g in order:
        tr = Phylo.read(os.path.join(OUT, "trees", f"{g}.nwk"), "newick")
        tr.ladderize(reverse=True)        # larger child first, all the way down
        trees[g] = tr

    # a tip's share of the circle is the share of the dataset it stands for:
    # its records (all sequences) or its 90 percent cluster (centroids), so
    # every fold's slice is proportional to its records
    weight = {t["tip"]: float(size[t["tip"]]) for t in tips}

    reach = RIM - INNER
    scale = None
    for g in order:
        multi = n_comp[g] > 1
        if len(groups[g]) == 1 or max(d["depth"] for d in walk(trees[g], multi, 0.0)) == 0:
            continue
        s = reach / max(d["depth"] for d in walk(trees[g], multi, 0.0))
        for _ in range(8):
            s = reach / max(d["depth"] for d in walk(trees[g], multi, BACKBONE_STEP / s, weight))
        scale = s if scale is None else min(scale, s)
    step = BACKBONE_STEP / scale

    def stub(tip, y):
        w = weight[tip]
        return {"name": tip, "tip": True, "parent": None, "depth": 0.0, "boot": None,
                "comp": None, "edge": "solid", "fork": "solid", "ntips": 1,
                "stub": True, "y": y + w / 2, "lo": y, "hi": y + w}

    gap_rows = 6.0 if V else GAP_ROWS
    top_gap = 14.0 if V else TOP_GAP_ROWS * (1.6 if TS > 1 else 1.0)   # large ring names need a wider opening
    slices, at = [], top_gap / 2
    named = lambda g: TS > 1 and len(groups[g]) > 5    # vertical, large text: folds of more than 5 tips are named
    for g in order:
        nodes = ([stub(groups[g][0]["tip"], 0.0)] if len(groups[g]) == 1
                 else walk(trees[g], n_comp[g] > 1, step, weight))
        used = max(d["hi"] for d in nodes)
        width = used
        shift = at + (width - used) / 2
        for d in nodes:
            for f in ("y", "lo", "hi"):
                d[f] += shift
        slices.append((g, nodes, at, at + width))
        at += width + gap_rows
    total = at - gap_rows + top_gap / 2

    # Two geometries share all of the drawing code below. Polar: `r` is the
    # radius and `a` the angle. Vertical: `r` is x and `a` is y (inches, the
    # top of the tree at 0, going down), so a tip's row is its y and a
    # "fork arc" is a vertical line.
    if V:
        ang = lambda y: -y * ROW_H
        xy = lambda r, a: (r, a)
        arc = lambda r, a0, a1: [(r, a0), (r, a1)]
    else:
        ang = lambda y: math.pi / 2 - 2 * math.pi * y / total
        xy = lambda r, a: (r * math.cos(a), r * math.sin(a))

        def arc(r, a0, a1):
            k = max(2, int(abs(a1 - a0) / math.radians(0.8)) + 2)
            return [xy(r, a) for a in np.linspace(a0, a1, k)]

    def sector(r0, r1, a0, a1, **kw):
        """Annular sector (polar) or box (vertical) between r0 and r1."""
        if V:
            return ax.add_patch(Rectangle((r0, min(a0, a1)), r1 - r0, abs(a1 - a0),
                                          linewidth=0, **kw))
        return ax.add_patch(Wedge((0, 0), r1, math.degrees(min(a0, a1)),
                                  math.degrees(max(a0, a1)), width=r1 - r0, **kw))

    def ring_cell(r0, w, a0, a1, color, z=1, ec="none", **kw):
        return sector(r0, r0 + w, a0, a1, facecolor=color, edgecolor=ec, zorder=z, **kw)

    def scwedge(layer, group, r0, w, a0, a1, fill, title=None):
        if V:
            sc.rect(layer, group, r0, min(a0, a1), r0 + w, max(a0, a1), fill, title=title)
        else:
            sc.wedge(layer, group, r0, w, a0, a1, fill, title=title)

    def rim(a, r):
        """Where text set outward from r at a goes: x, y, rotation, ha."""
        if V:
            return r, a, 0.0, "left"
        x, y = xy(r, a)
        left = math.cos(a) < 0
        return x, y, math.degrees(a) + (180 if left else 0), "right" if left else "left"

    cons = conservation(tips, rows, keyof)

    cv = d8.PieCanvas()
    ax = cv.ax
    fold_mid, ordered_tips = {}, []
    fold_reach = {}

    # rings: one per plastic (the fold, accessory domain and catalytic residue
    # rings were removed on the reviewer feedback of 2026-10-28)
    rings = [(p, RING_W) for p in PLASTICS]
    if args.rainbow:
        rings.append(("position", 0.07))
    r_at, ring_r = R_RING, {}
    for name, w in rings:
        if name == "position":
            r_at += RING_GAP * 2
        ring_r[name] = (r_at, w)
        r_at += w + RING_GAP
    r_out = r_at
    if V:
        cv.grow(INNER - 0.3, 0.9)
        cv.grow(r_out + 0.9, -total * ROW_H - 0.4)
    else:
        cv.grow(-r_out - 0.9, -r_out - 0.9)       # wedges do not register their extent
        cv.grow(r_out + 0.9, r_out + 0.9)

    has_structure = set(folds.annotations())
    fold_span, seed_angles = {}, collections.defaultdict(list)
    for g, nodes, lo, hi in slices:
        k = gkey[g]
        a_lo, a_hi = ang(lo), ang(hi)
        fold_mid[g] = (a_lo + a_hi) / 2
        fold_span[g] = (min(a_lo, a_hi), max(a_lo, a_hi))
        cv.tag(sector(INNER, RIM + 0.12, a_lo, a_hi, facecolor=folds.tint(k),
                      edgecolor="none", zorder=0), f"bg|{g}", "slice")
        scwedge("slices", k, INNER, RIM + 0.12 - INNER, a_lo, a_hi, folds.tint(k),
                title=folds.title(k))
        # empty ring tracks across the slice
        for name, (r0, w) in ring_r.items():
            if name == "position":
                continue
            cv.tag(ring_cell(r0, w, a_lo, a_hi, "#f1f1ee", z=0.5), f"rings|{g}",
                   f"{name} track")
            scwedge("ring tracks", k, r0, w, a_lo, a_hi, "#f1f1ee")

        for d in nodes:
            d["a"] = ang(d["y"])
            d["r"] = INNER + (BACKBONE_STEP if d.get("stub") else d["depth"] * scale)
            if d.get("stub"):
                cv.line(f"br|{g}", ("#555555", 0.55 * args.line_scale, 2), [xy(INNER, d["a"]), xy(d["r"], d["a"])])
                sc.line("branches", k, [xy(INNER, d["a"]), xy(d["r"], d["a"])], "#555555", 0.55)
        stems = sorted((d for d in nodes if d["name"].startswith("S_")), key=lambda d: d["lo"])
        for i, d in enumerate(stems):
            if i % 2:
                a0, a1 = ang(d["lo"]), ang(d["hi"])
                cv.tag(sector(INNER, RIM + 0.12, a0, a1,
                              facecolor=d8.shade(folds.tint(k), 0.955),
                              edgecolor="none", zorder=0.2),
                       f"bg|{g}", f"component {d['comp']}, second tint")
                scwedge("slices", k, INNER, RIM + 0.12 - INNER, a0, a1,
                        d8.shade(folds.tint(k), 0.955), title=f"component {d['comp']}")

        pts = []
        for d in nodes:
            if d["parent"] is not None:
                pts_ = [xy(d["parent"]["r"], d["a"]), xy(d["r"], d["a"])]
                cv.line(f"br|{g}", d["edge"], pts_)
                sc.line("branches", k, pts_, *line_style(d["edge"]))
            if not d["tip"]:
                pts_ = arc(d["r"], d["kids"][0]["a"], d["kids"][-1]["a"])
                cv.line(f"br|{g}", d["fork"], pts_)
                sc.line("branches", k, pts_, *line_style(d["fork"]))
                # support is shown by line style only (pale under 70); no
                # numbers, they crowded the tree
                continue
            t = tipinfo[d["name"]]
            recs = t["records"].split(";")
            frac = (d["y"] - top_gap / 2) / (total - top_gap)
            ordered_tips.append((t, frac, d["a"]))
            x, y = xy(d["r"], d["a"])
            n = size[d["name"]]
            pts.append((x, y, "o",
                        SUBSTRATE_COLOR[t["substrate_group"]], n, d["name"] in seeds))
            sc.marker("tips", k, x, y, "o", SUBSTRATE_COLOR[t["substrate_group"]],
                      6 + 7 * (n - 1) ** 0.85 if n > 1 else 6,
                      title=f"{t['enzyme_name'] or d['name']} | {t['substrate_group']} | "
                            f"component {t['component']} | {n} {size_label}")
            half = (ROW_H / 2 if V else math.pi / total) * (d["hi"] - d["lo"]) * 0.98   # half the tip's share
            # guide from tip to the rings, in the tip's position hue
            if args.rainbow:
                cv.line(f"guide|{g}", (rainbow(frac) + (0.55,), 0.35, 0.8),
                        [xy(d["r"] + 0.06, d["a"]), xy(R_RING - 0.02, d["a"])])
                sc.line("guides", k, [xy(d["r"] + 0.06, d["a"]), xy(R_RING - 0.02, d["a"])],
                        rainbow(frac) + (0.55,), 0.35)
                r0, w = ring_r["position"]
                cv.tag(ring_cell(r0, w, d["a"] - half, d["a"] + half, rainbow(frac), z=1),
                       f"rings|{g}", "position")
                scwedge("rings", k, r0, w, d["a"] - half, d["a"] + half, rainbow(frac),
                        title=f"position {frac:.3f}")
            have = plastics_of(rows, recs)
            for p in PLASTICS:
                if p in have:
                    r0, w = ring_r[p]
                    cv.tag(ring_cell(r0, w, d["a"] - half, d["a"] + half, RING_COLOR[p], z=1),
                           f"rings|{g}", p)
                    scwedge("rings", k, r0, w, d["a"] - half, d["a"] + half, RING_COLOR[p],
                            title=f"{p}: {t['enzyme_name']}")
        # tips: area by size (the named enzymes are not marked or labelled)
        by = collections.defaultdict(list)
        for x, y, mk, col, n, seed in pts:
            by[(mk, col)].append((x, y, n, seed))
            cv.grow(x, y)
        for (mk, col), p in sorted(by.items()):
            if args.no_tips:
                break
            a = np.array([(x, y) for x, y, _, _ in p])
            s = [6 + 7 * (n - 1) ** 0.85 if n > 1 else 6 for _, _, n, _ in p]
            cv.tag(ax.scatter(a[:, 0], a[:, 1], marker=mk, s=s, c=col, edgecolors=NODE_STROKE,
                              linewidths=0.3, zorder=5), f"tip|{g}", col)

    # --- ring names at 12 o'clock -------------------------------------------
    for name, (r0, w) in ring_r.items():
        lab = retext(f"ring|{name}", {"fold": "fold", "position": "position (tracking hue)"}
                     .get(name, name))
        col = RING_COLOR.get(name, INK) if name in PLASTICS else INK
        # polar: in the opening at 12 o'clock; vertical: above the top, reading up
        nx, ny, nrot, nha = (r0 + w / 2, 0.08, 90, "left") if V else (0, r0 + w / 2, 0, "center")
        cv.text("ringnames", nx, ny, lab, fontsize=F((5.2 if name in PLASTICS else 4.8) * TS),
                color=col, ha=nha, va="center", zorder=9, rotation=nrot, rotation_mode="anchor",
                fontweight="bold" if name in PLASTICS else "normal")
        if VFS:
            cv.grow(nx, ny + text_w(lab, VFS, name in PLASTICS) + 0.1)
        sc.text("ring names", "", nx, ny, lab, F((5.2 if name in PLASTICS else 4.8) * TS), col, ha=nha,
                rotation=nrot, weight="bold" if name in PLASTICS else "normal",
                tid=f"ring|{name}")
    tx, ty, trot, tha = (INNER + 0.08, 0.08, 90, "left") if V else (0, INNER + 0.25, 0, "center")
    if TS == 1:
        cv.text("ringnames", tx, ty, "tree", fontsize=5 * TS, color=MUTED, ha=tha, rotation=trot,
                rotation_mode="anchor")
        sc.text("ring names", "", tx, ty, "tree", 5 * TS, MUTED, ha=tha, rotation=trot,
                va="baseline" if not V else "center")

    # vertical, large text: each fold's name just above the top left of its own
    # tree, in the gap between folds (gap_rows leaves room for one line)
    left_ext = r_out + 0.12                    # polar: how far left anything outside the rings reaches
    if TS > 1 and not V:                       # polar: the name runs outward from the slice's middle
        fs_ = 5.5 * TS
        for g, nodes, lo, hi in slices:
            if not named(g):
                continue
            k = gkey[g]
            name_ = retext(f"foldname|{k}", re.sub(r" *[(][^)]*[)]", "", folds.title(k)))
            if len(name_) > 22:                # a long name on two lines, so the group bars stay close
                cut_ = name_.rfind(" ", 0, len(name_) // 2 + 4)
                name_ = name_[:cut_] + "\n" + name_[cut_ + 1:]
            xs_, ys_, rot_, ha_ = rim(fold_mid[g], r_out + 0.12)
            cv.text(f"foldname|{g}", xs_, ys_, name_, fontsize=fs_, color=folds.strong(k), ha=ha_,
                    va="center", fontweight="bold", zorder=9, rotation=rot_, rotation_mode="anchor")
            sc.text("fold labels", k, xs_, ys_, name_, fs_, folds.strong(k), ha=ha_, rotation=rot_,
                    weight="bold", tid=f"foldname|{k}")
            reach_ = r_out + 0.12 + max(len(l_) for l_ in name_.split("\n")) * fs_ * 0.62 / 72
            cv.grow(*xy(reach_, fold_mid[g]))
            left_ext = max(left_ext, -xy(reach_, fold_mid[g])[0])
            fold_reach[g] = reach_
    # vertical: no fold names above the trees (removed 2026-11-05)

    # --- fold types: a bar and a name for each group of folds ---------------
    # vertical: right of the rings, so the grid of structures moves over by
    # type_w. Vertical mechanism groups (left_side): the bars, numbers and names
    # sit left of the tree and the structures left of them, so only the legend
    # is right of the rings. The order there reads name, number, bar, tree.
    type_w = 0.0
    name_y = {}
    badge_on_bar = few and args.type_bars    # the number sits on the group's bar, before its name
    left_side = V and badge_on_bar
    bar_x = INNER - 0.32 if left_side else r_out + 0.12   # vertical: the group bars' left edge
    name_a = {}                                    # polar: the angle each group's name and number sit at
    rb = r_out + 0.62
    if args.type_bars:
        type_span = collections.defaultdict(list)
        for g, nodes, lo, hi in slices:
            type_span[folds.type_rank(gkey[g])] += [ang(lo), ang(hi)]
        name_y = {r_: (min(v_) + max(v_)) / 2 for r_, v_ in type_span.items()}
        rb = max([r_out + 0.62] + [r_ + 0.25 for r_ in fold_reach.values()])   # polar: bars clear of the fold names
        if not V:
            cv.grow(-rb - 0.9, -rb - 0.9)      # the bars are wedges, which do not register their extent
            cv.grow(rb + 0.9, rb + 0.9)
        bar_mid = dict(name_y)       # each group's bar middle, before the names are pushed apart
        badge_col = {}               # --no-group-names: 1 for a circle in the outer column
        if V and args.no_group_names:
            # numbers alone: each circle level with the middle of its bar; a circle
            # that would overlap one already placed moves to the next column out
            # (2026-11-14), joined to its bar by a horizontal line under the circles
            placed_ = []
            for r_ in sorted(name_y, key=lambda r_: -name_y[r_]):
                col_ = 0
                while any(c_ == col_ and abs(y_ - name_y[r_]) < BADGE * 1.08 / 72 for c_, y_ in placed_):
                    col_ += 1
                badge_col[r_] = col_
                placed_.append((col_, name_y[r_]))
        elif V or TS > 1:   # push neighbouring names apart so two short bars do not collide
            sep_ = F(8 * TS) * (1.3 if V else 1.9) / 72 / (1.0 if V else rb + 0.4)
            ys_ = sorted(name_y, key=lambda r_: -name_y[r_])
            for prev_, r_ in zip(ys_, ys_[1:]):
                name_y[r_] = min(name_y[r_], name_y[prev_] - sep_)
        for rank_, (tname, _) in enumerate(folds.TYPES):
            if rank_ not in type_span:
                continue
            a0, a1 = min(type_span[rank_]), max(type_span[rank_])
            col = (folds.group_dark(rank_) if args.grouping == "mechanism"
                   else TYPE_COLOR[rank_ % len(TYPE_COLOR)])
            if V:
                bx0 = bar_x
                cv.tag(sector(bx0, bx0 + 0.07, a0, a1, facecolor=col, edgecolor="none", zorder=2),
                       f"type|{rank_}", tname)
                scwedge("groups", "", bx0, 0.07, a0, a1, col, title=tname)
                gap_ = 0.1 + 0.1 * TS + (BADGE / 72 + 0.1 if badge_on_bar else 0.0)   # bar to name
                if left_side:                  # mirrored: the name ends where the number begins
                    lx, ly, lrot, lha = bx0 - gap_, name_y[rank_], 0.0, "right"
                else:
                    lx, ly, lrot, lha = bx0 + 0.07 + gap_, name_y[rank_], 0.0, "left"
                if TS > 1 and len(tname) > 18 and abs(a1 - a0) > 1.2:   # a tall bar: two lines
                    cut = tname.rfind(" ", 0, len(tname) // 2 + 4)
                    tname = tname[:cut] + "\n" + tname[cut + 1:]
                wide = max(len(t_) for t_ in tname.split("\n"))
                if left_side and args.no_group_names:     # the number alone; the name is under the structure
                    more_ = max(badge_col.values(), default=0)      # columns beyond the first
                    type_w = max(type_w, BADGE / 72 * (1 + more_ * 1.08) + 0.25)
                    cv.grow(bx0 - 0.12 - BADGE / 72 * (1 + more_ * 1.08), ly)
                    continue
                type_w = max(type_w, 0.36 + 0.1 * TS + (BADGE / 72 + 0.1 if badge_on_bar else 0.0) + wide * F(8 * TS) * 0.6 / 72)
                cv.grow(lx + (-1 if left_side else 1) * wide * F(8 * TS) * 0.6 / 72, ly)   # room for the name
            else:
                cv.tag(sector(rb, rb + 0.07, a0, a1, facecolor=col, edgecolor="none", zorder=2),
                       f"type|{rank_}", tname)
                scwedge("groups", "", rb, 0.07, a0, a1, col, title=tname)
                am = (a0 + a1) / 2
                bshift = BADGE / 72 + 0.1 if badge_on_bar else 0.0   # the number sits between the bar and the name
                name_a[rank_] = (a0 + a1) / 2 if rank_ not in name_y else name_y[rank_]
                am = name_a[rank_]
                fs_in = 8 * TS / 72
                span_in = abs(a1 - a0) * (rb + 0.2)
                if span_in * 0.92 > len(tname) * 8 * TS * 0.58 / 72:    # fits along the bar
                    lx, ly = xy(rb + 0.17 + 0.1 * TS + bshift, am)
                    lrot = math.degrees(am) - 90
                    if math.sin(am) < 0:               # lower half: keep it upright
                        lrot += 180
                    lha = "center"
                    reach_ = rb + 0.17 + 0.1 * TS + bshift + fs_in
                    left_ext = max(left_ext, -xy(reach_, am)[0])
                else:                                               # radial, outside the bar
                    lx, ly, lrot, lha = rim(am, rb + 0.14 + bshift)
                    reach_ = rb + 0.2 + bshift + len(tname) * 8 * TS * 0.58 / 72
                    cv.grow(*xy(reach_, am))   # room for the name
                    left_ext = max(left_ext, -xy(reach_, am)[0])
            cv.tag(ax.text(lx, ly, tname, fontsize=F(8 * TS), color=col, fontweight="bold", zorder=9,
                           rotation=lrot, rotation_mode="anchor", va="center", ha=lha),
                   f"type|{rank_}", tname)
            sc.text("groups", "", lx, ly, tname, F(8 * TS), col, ha=lha, rotation=lrot, weight="bold",
                    tid=f"type|{rank_}")

    # --- structures and fold labels, outside -------------------------------
    A = folds.annotations()
    with open(os.path.join(cards.OUTDIR, "scenes.json")) as fh:
        inside = {sc_["key"]: sc_["contents"] for sc_ in json.load(fh)}
    # every fold gets a slot around the ring; a fold with no annotated
    # structure gets its name there and no picture
    def anchor_angle(g):
        """Where the leader leaves the slice: the middle, unless a seed name
        sits there, then the point of the slice furthest from every seed."""
        lo_, hi_ = fold_span[g]
        if not seed_angles[g]:
            return fold_mid[g]
        cand = np.linspace(lo_ + (hi_ - lo_) * 0.1, hi_ - (hi_ - lo_) * 0.1, 41)
        gap = [min(abs(c - a) for a in seed_angles[g]) - 0.05 * abs(c - fold_mid[g])
               for c in cand]
        return float(cand[int(np.argmax(gap))])

    # folds with a structure get a slot around the ring; one without gets its
    # name in small grey text just outside the rings, out of the way
    items = []
    group_span = collections.defaultdict(list)       # mechanism groups: one badge per group,
    for g, nodes, lo, hi in slices:                  # at the middle of the whole group
        group_span[folds.type_rank(gkey[g])] += [ang(lo), ang(hi)]
    for g, nodes, lo, hi in slices:
        k = gkey[g]
        if k in A and not args.no_structures and (
                args.grouping != "mechanism" or k in folds.GROUP_PICTURED):
            if args.grouping == "mechanism":
                sp_ = group_span[folds.type_rank(k)]
                items.append((g, (min(sp_) + max(sp_)) / 2))
            else:
                items.append((g, anchor_angle(g)))
            continue
        if V or few:                      # the group names sit outside the rings instead
            continue
        # on the slice's edge facing the opening at 12 o'clock when it borders
        # it (no leader runs there), otherwise at the middle of the slice
        lo_a, hi_a = fold_span[g]
        if not V and abs(math.cos(lo_a)) < math.sin(math.pi * top_gap / total) + 0.05:
            a0 = lo_a - 0.012
        else:
            a0 = anchor_angle(g)          # clear of the seed names
        xs, ys, rot, ha = rim(a0, r_out + 0.06)
        lab = retext(f"fold|{k}|rim", folds.title(k))
        cv.tag(ax.text(xs, ys, lab, fontsize=4.2, color=MUTED, zorder=9, rotation=rot,
                       rotation_mode="anchor", va="center", ha=ha,
                       fontstyle="italic"), f"label|{g}", lab)
        sc.text("fold labels", k, xs, ys, lab, 4.2, MUTED, ha=ha,
                rotation=rot, tid=f"fold|{k}|rim")
    # where each structure goes: in a grid beside the figure, in order, each
    # matched to its fold by the same number on the slice and on the picture
    # (default and --vertical); or around the ring on leaders (--leaders)
    sz = (3.6 if few else PANEL_IN) if args.panel else STRUCT_IN
    if V:
        want = [-a for _, a in items]                    # top to bottom
    else:
        want = [-(a - math.pi / 2) % (2 * math.pi) for _, a in items]   # clockwise from 12
    if args.panel:
        got, lanes, rank, r_struct = list(want), [0] * len(want), {}, 0.0
    else:
        opening = math.pi * top_gap / total          # half the gap at 12 o'clock
        # two staggered rows: neighbours alternate inner and outer, so a structure
        # needs only half the arc it would in one row and can sit by its slice
        # 0.6 rather than 0.5 of a one-row spacing, so a leader running out to
        # the outer row passes clear of the captions of its inner neighbours
        got = place_structures(want, (STRUCT_IN + 0.6) / STRUCT_R * 0.6, opening,
                               2 * math.pi - opening)
        rank = {i: n for n, i in enumerate(sorted(range(len(got)), key=lambda i: got[i]))}
        lanes = assign_lanes(want, got)
        # the inner row clears the outermost lane, caption included (a caption
        # below a structure at 12 o'clock faces the rings)
        r_struct = max(STRUCT_R, r_out + LANE_IN + max(lanes) * LANE_GAP + STRUCT_IN / 2 + 0.8)
    cols = 3 if V and few else 2 if few else VPANEL_COLS if V else PANEL_COLS   # vertical mechanism: 1 2 3 over 4 5 6 (2026-11-14)
    n_rows = math.ceil(len(items) / cols)
    def catalytic_of(k):
        """The catalytic residues drawn in k's structure, as the renders label
        them (scenes.json), without chain prefixes, in the sheet's order."""
        if not hasattr(catalytic_of, "sites"):
            with open(os.path.join(cards.OUTDIR, "scenes.json")) as fh:
                catalytic_of.sites = {s_["key"]: s_["sites"] for s_ in json.load(fh)}
        names_ = [x_["label"].split(":")[-1] for x_ in catalytic_of.sites.get(k, [])
                  if x_["kind"] == "primary" and x_["ok"]]
        return list(dict.fromkeys(names_))

    def captions(k, pictured):
        """The title and the lines under a structure; wrapped at 30 pt text."""
        if args.grouping == "mechanism":       # the group named, its pictured fold beneath
            title = retext(f"fold|{k}|name", folds.type_name(k))
            sub = retext(f"fold|{k}|sub", f"{re.sub(r' *[(][^)]*[)]', '', folds.title(k))}\n"
                         f"{folds.structure_code(k)}" if pictured else "")
            if V:                               # vertical: the enzyme names go in the caption (2026-11-14)
                sub = retext(f"fold|{k}|sub", folds.structure_code(k) if pictured else "")
                if args.no_ids:                  # the IDs too (2026-11-14)
                    sub = ""
                if args.residues and pictured:   # the catalytic residues instead, in the group colour
                    res_ = catalytic_of(k)
                    sub = " \u00b7 ".join(res_)
                    if len(sub) > 24:            # too long for one line: two even lines, broken between residues
                        h_ = (len(res_) + 1) // 2
                        sub = " \u00b7 ".join(res_[:h_]) + "\n" + " \u00b7 ".join(res_[h_:])
        else:
            title = retext(f"fold|{k}|name", folds.title(k))
            sub = retext(f"fold|{k}|sub", f"{folds.structure_code(k)}\n{folds.mechanism(k)}"
                         if pictured else "")
        if VFS:
            wrap_ = lambda t_: "\n".join(textwrap.fill(l_, 24) for l_ in t_.split("\n"))
            title, sub = wrap_(title), wrap_(sub)
        return title, sub

    cell_h = sz + 0.85
    if few:                                            # the grid spans the height of the tree
        cell_h = (total * ROW_H if V else 2 * r_out) / n_rows
        sz = min(cell_h - (0.5 * TS + 0.1 if TS > 1 else 0.85), 5.5)
    cell_w = sz + 0.45
    if VFS:          # room under each picture for its wrapped caption, and beside it for the widest line
        caps_ = [captions(gkey[g], gkey[g] in A and not args.no_structures) for g, _ in items]
        cap_h = max((t_.count("\n") + 1) * 1.25 + (s_.count("\n") + 1 if s_ else 0) * 1.3
                    for t_, s_ in caps_) * VFS / 72 + 0.45
        sz = min(cell_h - cap_h, 5.5)
        cell_w = max(sz + 0.45, max(max(text_w(t_, VFS, True), text_w(s_, VFS))
                                    for t_, s_ in caps_) + 0.5)
        if args.fill_cells:                    # the cell is set by the widest caption; fill it
            sz = min(cell_h - cap_h, cell_w - 0.45)
        if args.fit_content and V and few and args.panel:
            # rows from just under the column headers (drawn already) to the scale
            # bar's level; columns as wide as the widest protein needs to fill its
            # row's height, so no protein is held back by a narrow caption
            fit_top = cv.box[3] - 0.1
            fit_row_h = (fit_top - (ang(slices[-1][3]) - 1.0)) / n_rows
            need_w, fits_ = 0.0, []
            for g_, _ in items:
                k_ = gkey[g_]
                if k_ not in A or args.no_structures:
                    continue
                a_ = plt.imread(os.path.join(cards.OUTDIR, "full", f"{cards.slug(k_)}.png"))[..., 3]
                ys_, xs_ = np.where(a_ > 0.05)
                fx_, fy_ = (xs_.max() - xs_.min()) / a_.shape[1], (ys_.max() - ys_.min()) / a_.shape[0]
                need_w = max(need_w, (fit_row_h - cap_h - 0.15) / fy_ * fx_ + 0.3)
                fits_.append((fx_, fy_))
            # but the grid takes at most GRID_SHARE of the width: what the
            # structures gain the tree loses once the figure is printed at a fixed width
            rest_ = (r_out - bar_x) + BADGE / 72 + 0.9           # numbers, tree and plastic columns
            cap_w = args.grid_share / (1 - args.grid_share) * rest_ / cols
            cell_w = max(cell_w, min(need_w, cap_w))
            # how tall the rows have to be for the widest-fitting protein: the
            # caller can shorten the tree (--row-h) until there is no height to spare
            need_h = max(min((cell_w - 0.3) / fx_, (fit_row_h - cap_h - 0.15) / fy_) * fy_
                         for fx_, fy_ in fits_) + cap_h + 0.15
            print(f"fit rows: need {need_h:.2f} in, have {fit_row_h:.2f} in, "
                  f"tree {total * ROW_H:.2f} in")
    if V:
        grid_x0 = r_out + ((0.33 if TS > 1 else 0.5) + type_w if type_w else 1.1)   # left edge of the grid
        # left_side: the grid's right edge, as far left of the bars as the
        # grid sits right of them otherwise
        grid_x1 = bar_x + 0.07 - (0.21 + (0.33 if TS > 1 else 0.5) + type_w)
        grid_y0 = -total * ROW_H / 2 + n_rows * cell_h / 2 - sz / 2   # centre of the top row
        if few:
            grid_y0 = -sz / 2 - 0.1                    # a picture's top just under the tree's
    else:
        grid_x1 = -(left_ext if TS > 1 else r_out) - (0.33 if TS > 1 else 1.3)   # right edge of the grid
        grid_y0 = n_rows * cell_h / 2 - sz / 2
        if few:
            grid_y0 = r_out - sz / 2 - 0.1
    # a badge per fold beside its slice (outside the ring, or left of the
    # tree), spread so the small folds' badges do not overlap; a short tick
    # joins a badge that has moved to its slice
    if args.panel:
        bsep = (BADGE + 1.5) / 72 / (1.0 if V else r_out + 0.2)
        badge_at = spread_angles(want, bsep)
    INSET = 0.46                   # --insets: close up diameter, as a share of the picture
    INSET_FIT = 0.72               # --insets --fit-content: as a share of the column; may overlay the protein
    cell_of = {}                   # --fit-content: k -> (column left, column right, row top)
    even_r = []                    # --insets --fit-content: the one radius every close up gets
    pending_insets = []            # drawn after every structure is placed, so the size can be shared
    no_inset = set()               # --min-inset-mag: structures drawn without a close up
    alpha_box = {}

    def content_box(k, view):
        """The drawn extent of a render as fractions of its width (x0, y0, x1,
        y1, y down), from its alpha channel."""
        if (k, view) not in alpha_box:
            a_ = plt.imread(os.path.join(cards.OUTDIR, view, f"{cards.slug(k)}.png"))[..., 3]
            ys_, xs_ = np.where(a_ > 0.05)
            n_ = a_.shape[1]
            alpha_box[(k, view)] = (xs_.min() / n_, ys_.min() / n_, xs_.max() / n_, ys_.max() / n_)
        return alpha_box[(k, view)]

    views_ = {}
    if args.insets or args.label_residues:
        with open(os.path.join(cards.OUTDIR, "views.json")) as fh:
            views_ = json.load(fh)

    def inset_at(k, cx, cy, size_in):
        """Where the ring around the site (rx, ry, rr) and the round close up
        (ix, iy, ir) go, in inches. The ring is where the zoom window falls in
        the full render (views.json); the close up takes the corner of the
        picture farthest from it (the badge has the top left)."""
        ring_ = views_[cards.slug(k)]["full"]["ring"]
        px = size_in / views_[cards.slug(k)]["full"]["size"]          # inches per css px
        rx, ry = cx - size_in / 2 + ring_["x"] * px, cy + size_in / 2 - ring_["y"] * px
        rr = ring_["r"] * px
        ir = INSET * size_in / 2
        # the top right corner (the badge has the top left, the caption the
        # bottom), pushed up and out along the line from the site just far
        # enough to clear the ring
        ix, iy = cx + size_in / 2 - ir, cy + size_in / 2 - ir
        if k in cell_of:                 # --fit-content: a top corner of the column instead,
            # never above the row (the row above has its captions there) and never
            # over the ring: as large as fits clear of it, up to INSET_FIT of the
            # column, in whichever top corner allows the larger close up
            cl_, cr_, rt_ = cell_of[k]
            gap_ = 0.12
            best_ = None
            for side_ in (1, -1):        # right first, so a tie keeps the right
                at_ = lambda r_: ((cr_ - r_) if side_ == 1 else (cl_ + r_), rt_ - r_)
                clear_ = lambda r_: math.hypot(at_(r_)[0] - rx, at_(r_)[1] - ry) >= rr + r_ + gap_
                lo_, hi_ = 0.0, INSET_FIT * (cr_ - cl_) / 2
                if even_r:               # one size for every close up (the smallest of the maxima)
                    hi_ = min(hi_, even_r[0])
                if clear_(hi_):
                    lo_ = hi_
                else:
                    for _ in range(40):
                        mid_ = (lo_ + hi_) / 2
                        lo_, hi_ = (mid_, hi_) if clear_(mid_) else (lo_, mid_)
                if best_ is None or lo_ > best_[0] + 1e-6:
                    best_ = (lo_, *at_(lo_))
            ir, ix, iy = best_
            return rx, ry, rr, ix, iy, ir
        dx, dy = ix - rx, iy - ry
        dn = math.hypot(dx, dy) or 1
        short = rr + ir + 0.25 - dn
        if short > 0:
            ix, iy = ix + dx / dn * short, iy + dy / dn * short
        return rx, ry, rr, ix, iy, ir

    def badge_at_(k, cx, cy, size_in, view):
        """The structure's badge at the top left of the protein's visible extent,
        so it sits by its own structure rather than at the corner of the square
        picture, which can be nearer a neighbour. Vertical mechanism figure only."""
        a_ = plt.imread(os.path.join(cards.OUTDIR, view, f"{cards.slug(k)}.png"))[..., 3]
        ys_, xs_ = np.where(a_ > 0.05)
        n_ = a_.shape[1]
        left = cx - size_in / 2 + xs_.min() / n_ * size_in
        top = cy + size_in / 2 - ys_.min() / n_ * size_in
        br_ = BADGE / 144                                   # badge radius, inches
        # never left of its own cell, so it cannot drift towards a neighbour
        return max(left - 0.2 * br_, cx - cell_w / 2 + br_), top + 0.2 * br_

    def label_protein(k, g, cx, cy, size_in):
        """Name the catalytic residues on the whole protein: each name outward
        from the middle of the site, beyond its farthest residue, pushed clear of
        the other names, with a thin line in the group colour from the name's
        edge to the residue's tip (positions from views.json, full view)."""
        import matplotlib.patheffects as pe
        z_ = views_[cards.slug(k)]["full"]
        n_ = z_["size"]
        fs_ = F(6.6 * TS)
        bold_ = FontProperties(family=plt.rcParams["font.family"], weight="bold")
        labs = []
        for l_ in z_.get("labels", []):
            if l_["kind"] != "primary" or l_["text"].endswith("?"):
                continue
            t_ = l_["text"].split(":")[-1]               # B:T267 -> T267
            e_ = TextPath((0, 0), t_, size=fs_, prop=bold_).get_extents()
            labs.append({"t": t_, "ax": cx - size_in / 2 + l_["ax"] / n_ * size_in,
                         "ay": cy + size_in / 2 - l_["ay"] / n_ * size_in,
                         "w": e_.width / 72 + 0.08, "h": e_.height / 72 + 0.08})
        if not labs:
            return
        mx_ = sum(l_["ax"] for l_ in labs) / len(labs)
        my_ = sum(l_["ay"] for l_ in labs) / len(labs)
        reach_ = max(math.hypot(l_["ax"] - mx_, l_["ay"] - my_) for l_ in labs)
        for i_, l_ in enumerate(labs):                   # outward from the site's middle
            dx_, dy_ = l_["ax"] - mx_, l_["ay"] - my_
            dn_ = math.hypot(dx_, dy_)
            if dn_ < 1e-6:                               # one residue: to its upper left
                dx_, dy_, dn_ = -1.0, 0.6, math.hypot(1.0, 0.6)
            dx_, dy_ = dx_ / dn_, dy_ / dn_
            r_ = reach_ + 0.8 + max(l_["w"], l_["h"]) / 2      # past the sticks, not just the tips
            l_["x"], l_["y"] = mx_ + dx_ * r_, my_ + dy_ * r_
        for _ in range(300):                             # names apart from each other and off the residues
            for i_, a_ in enumerate(labs):
                for b_ in labs[i_ + 1:]:
                    ox_ = (a_["w"] + b_["w"]) / 2 + 0.05 - abs(a_["x"] - b_["x"])
                    oy_ = (a_["h"] + b_["h"]) / 2 + 0.05 - abs(a_["y"] - b_["y"])
                    if ox_ > 0 and oy_ > 0:
                        if oy_ < ox_:
                            sh_ = (1 if a_["y"] >= b_["y"] else -1) * oy_ / 2
                            a_["y"] += sh_; b_["y"] -= sh_
                        else:
                            sh_ = (1 if a_["x"] >= b_["x"] else -1) * ox_ / 2
                            a_["x"] += sh_; b_["x"] -= sh_
                for c_ in labs:                          # and clear of every residue tip
                    ox_ = a_["w"] / 2 + 0.45 - abs(a_["x"] - c_["ax"])
                    oy_ = a_["h"] / 2 + 0.45 - abs(a_["y"] - c_["ay"])
                    if ox_ > 0 and oy_ > 0:
                        dx_, dy_ = a_["x"] - mx_, a_["y"] - my_
                        dn_ = math.hypot(dx_, dy_) or 1
                        a_["x"] += dx_ / dn_ * 0.05; a_["y"] += dy_ / dn_ * 0.05
        col = folds.strong(k)
        for l_ in labs:
            # the line ends at the edge of the name's box, on the way to the residue
            dx_, dy_ = l_["ax"] - l_["x"], l_["ay"] - l_["y"]
            t_box = min(abs((l_["w"] / 2) / dx_) if abs(dx_) > 1e-9 else 1e9,
                        abs((l_["h"] / 2) / dy_) if abs(dy_) > 1e-9 else 1e9, 1.0)
            x0_, y0_ = l_["x"] + dx_ * t_box, l_["y"] + dy_ * t_box
            cv.tag(ax.add_line(Line2D([x0_, l_["ax"]], [y0_, l_["ay"]], color=col,
                                      lw=0.9 * args.line_scale, zorder=7, solid_capstyle="round",
                                      path_effects=[pe.withStroke(linewidth=0.9 * args.line_scale + 2.5,
                                                                  foreground="white")])),
                   f"struct|{g}", f"line to {l_['t']}")
            cv.tag(ax.text(l_["x"], l_["y"], l_["t"], fontsize=fs_, color=col, fontweight="bold",
                           ha="center", va="center", zorder=8,
                           path_effects=[pe.withStroke(linewidth=fs_ * 0.3, foreground="white")]),
                   f"struct|{g}", f"residue {l_['t']}")

    def label_inset(k, g, ix, iy, ir):
        """Name the catalytic residues in a close up, in the figure's text size
        with a white outline. Positions are the residues' tips in the zoom
        render (views.json); each name starts just outside its residue, away
        from the middle of the site, and the names are pushed apart and kept
        inside the circle."""
        import matplotlib.patheffects as pe
        z_ = views_[cards.slug(k)]["zoom"]
        n_ = z_["size"]
        fs_ = F(6.6 * TS)
        bold_ = FontProperties(family=plt.rcParams["font.family"], weight="bold")
        labs = []
        for l_ in z_["labels"]:
            if l_["kind"] != "primary":
                continue
            px_ = ix - ir + l_["ax"] / n_ * 2 * ir
            py_ = iy + ir - l_["ay"] / n_ * 2 * ir
            if math.hypot(px_ - ix, py_ - iy) > 0.92 * ir:
                continue
            t_ = l_["text"].split(":")[-1]               # B:T267 -> T267
            e_ = TextPath((0, 0), t_, size=fs_, prop=bold_).get_extents()
            labs.append({"t": t_, "ax": px_, "ay": py_, "w": e_.width / 72, "h": e_.height / 72})
        if not labs:
            return
        mx_ = sum(l_["ax"] for l_ in labs) / len(labs)
        my_ = sum(l_["ay"] for l_ in labs) / len(labs)
        for l_ in labs:                                  # start beside the residue, away from the site's middle
            dx_, dy_ = l_["ax"] - mx_, l_["ay"] - my_
            dn_ = math.hypot(dx_, dy_)
            dx_, dy_ = (dx_ / dn_, dy_ / dn_) if dn_ > 1e-6 else (0.0, 1.0)
            l_["x"] = l_["ax"] + dx_ * (l_["w"] / 2 + 0.12)
            l_["y"] = l_["ay"] + dy_ * (l_["h"] / 2 + 0.12)
        for _ in range(200):                             # apart, and inside the circle
            for i_, a_ in enumerate(labs):
                for b_ in labs[i_ + 1:]:
                    ox_ = (a_["w"] + b_["w"]) / 2 + 0.06 - abs(a_["x"] - b_["x"])
                    oy_ = (a_["h"] + b_["h"]) / 2 + 0.06 - abs(a_["y"] - b_["y"])
                    if ox_ > 0 and oy_ > 0:
                        if oy_ < ox_:
                            sh_ = (1 if a_["y"] >= b_["y"] else -1) * oy_ / 2
                            a_["y"] += sh_; b_["y"] -= sh_
                        else:
                            sh_ = (1 if a_["x"] >= b_["x"] else -1) * ox_ / 2
                            a_["x"] += sh_; b_["x"] -= sh_
            for l_ in labs:
                lim_ = ir * 0.95 - math.hypot(l_["w"], l_["h"]) / 2
                d_ = math.hypot(l_["x"] - ix, l_["y"] - iy)
                if d_ > lim_ > 0:
                    l_["x"] = ix + (l_["x"] - ix) * lim_ / d_
                    l_["y"] = iy + (l_["y"] - iy) * lim_ / d_
        for l_ in labs:
            cv.tag(ax.text(l_["x"], l_["y"], l_["t"], fontsize=fs_, color=folds.strong(k),
                           fontweight="bold", ha="center", va="center", zorder=7,
                           path_effects=[pe.withStroke(linewidth=fs_ * 0.28, foreground="white")]),
                   f"struct|{g}", f"residue {l_['t']}")

    def draw_inset(k, g, cx, cy, size_in):
        """The close up, its ring and the line joining them (inset_at)."""
        rx, ry, rr, ix, iy, ir = inset_at(k, cx, cy, size_in)
        zimg = plt.imread(os.path.join(cards.OUTDIR, "zoom", f"{cards.slug(k)}.png")).copy()
        n_ = zimg.shape[0]
        yy, xx = np.mgrid[:n_, :n_]
        zimg[..., 3] *= ((xx - n_ / 2) ** 2 + (yy - n_ / 2) ** 2 <= (n_ / 2) ** 2)
        col = folds.strong(k)
        lw_ = 0.12 * F(1)                                             # 4.5 pt at 37.5 pt text
        dx, dy = ix - rx, iy - ry
        dn = math.hypot(dx, dy) or 1
        if dn > rr + ir:                                              # the joining line, rim to rim
            cv.tag(ax.add_line(Line2D([rx + dx / dn * rr, ix - dx / dn * ir],
                                      [ry + dy / dn * rr, iy - dy / dn * ir],
                                      color=col, lw=lw_ * 0.6, zorder=5)), f"struct|{g}", "inset line")
        for x_, y_, r_, fc, z_ in ((rx, ry, rr, "none", 6), (ix, iy, ir, SURFACE, 5)):
            cv.tag(ax.add_patch(plt.Circle((x_, y_), r_, fc=fc, ec=col, lw=lw_, zorder=z_)),
                   f"struct|{g}", "inset ring")
        zi = OffsetImage(zimg, zoom=2 * ir * 72 / n_)
        cv.tag(ax.add_artist(AnnotationBbox(zi, (ix, iy), frameon=False, zorder=5.5)),
               f"struct|{g}", f"close up {folds.title(k)}")
        cv.tag(ax.add_patch(plt.Circle((ix, iy), ir, fc="none", ec=col, lw=lw_, zorder=6)),
               f"struct|{g}", "inset rim")
        if even_r:                       # --fit-content: the catalytic residues named, at print size
            label_inset(k, g, ix, iy, ir)

    for idx, ((g, a_mid), w_c, clockwise, lane) in enumerate(zip(items, want, got, lanes)):
        k = gkey[g]
        a_mid = -w_c if V else math.pi / 2 - w_c            # unwrapped, so arcs take the short way
        a = math.pi / 2 - clockwise
        number = sorted(range(len(want)), key=lambda i: want[i]).index(idx) + 1
        if args.panel:
            row, col = divmod(number - 1, cols)
            if V and not left_side:
                cx = grid_x0 + (col + 0.5) * cell_w
            else:
                cx = grid_x1 - (cols - col - 0.5) * cell_w
            cy = grid_y0 - row * cell_h
        else:
            cx, cy = xy(r_struct + (STAGGER if rank[idx] % 2 else 0.0), a)
        size_in, view = sz, "full"
        fit_bottom = None
        cap_x = cx                       # the caption's centre: the column's, also under --fit-content
        if args.fit_content and V and few and args.panel and k in A and not args.no_structures:
            # the row: from the column headers down to the scale bar, captions at its foot
            y_hi_ = fit_top - row * fit_row_h
            fit_bottom = y_hi_ - fit_row_h + cap_h
            x0_, y0_, x1_, y1_ = content_box(k, view)
            # height: the row less the caption, and less part of the close up at
            # the column's top right, so the close up covers little of the protein
            room_h = fit_row_h - cap_h - 0.15
            size_in = min((cell_w - 0.3) / (x1_ - x0_), room_h / (y1_ - y0_))
            if args.insets:
                # a close up only where it magnifies enough: its radius over the
                # radius of the ring it enlarges (sized for the reserved room)
                ring_r_ = views_[cards.slug(k)]["full"]["ring"]["r"] / views_[cards.slug(k)]["full"]["size"]
                ir_ = INSET_FIT * (cell_w - 0.2) / 2
                room_h = fit_row_h - cap_h - 0.15        # no room kept: the close up may overlay the protein
                size_r = min((cell_w - 0.3) / (x1_ - x0_), room_h / (y1_ - y0_))
                mag_ = ir_ / (ring_r_ * size_r)
                if mag_ >= args.min_inset_mag:
                    size_in = size_r
                else:
                    no_inset.add(k)
                if k in no_inset:
                    print(f"  inset {folds.type_name(k):<32} magnifies {mag_:.1f}x, left out")
            cx = cx - size_in * ((x0_ + x1_) / 2 - 0.5)          # drawn extent centred in the column
            cy = fit_bottom - size_in / 2 + size_in * y1_         # its foot on the row's foot
            cx_cell = cx + size_in * ((x0_ + x1_) / 2 - 0.5)
            cap_x = cx_cell
            cell_of[k] = (cx_cell - cell_w / 2 + 0.1, cx_cell + cell_w / 2 - 0.1, y_hi_)
        moved = layout.get("structures", {}).get(k)
        if moved:
            cx, cy = moved["x"], moved["y"]
            size_in, view = moved.get("size", sz), moved.get("view", "full")
        r_lane = r_out + LANE_IN + lane * LANE_GAP
        if args.panel:
            if left_side:
                bx, by = bar_x - 0.12 - BADGE / 144, name_y[folds.type_rank(k)]
                tick = []
                if args.no_group_names:
                    col_ = badge_col.get(folds.type_rank(k), 0)
                    bx -= BADGE / 72 * 1.08 * col_     # an outer column
                    if col_:                            # straight across to the bar, under any circle on the way
                        tick = [(bx, by), (bar_x + 0.035, by)]
                else:
                    mid_ = bar_mid[folds.type_rank(k)]
                    if abs(by - mid_) > 0.05:  # moved off its bar: a short tick back to the bar's middle
                        tick = [(bx + BADGE / 144, by), (bar_x - 0.02, mid_)]
            elif V and badge_on_bar:
                bx, by = bar_x + 0.07 + 0.12 + BADGE / 144, name_y[folds.type_rank(k)]
                tick = []
            elif not V and badge_on_bar:
                bx, by = xy(rb + 0.19 + BADGE / 144, name_a[folds.type_rank(k)])
                tick = []
            elif V:
                bx, by = INNER - 0.25, -badge_at[idx]
                tick = [(INNER - 0.02, a_mid), (bx + BADGE / 144, by)]
            else:
                ab = math.pi / 2 - badge_at[idx]
                bx, by = xy(r_out + 0.2, ab)
                tick = [xy(r_out + 0.04, a_mid), xy(r_out + 0.2 - BADGE / 144, ab)]
            lead = []
            if tick:
                cv.line(f"leader|{g}", (folds.strong(k), 1.2 * args.line_scale if args.no_group_names
                                        else 0.5, 3), tick)
                sc.line("leaders", k, tick, folds.strong(k), 0.5)
            sbx, sby = cx - size_in / 2 + 0.1, cy + size_in / 2 - 0.1
            if V and few and k in A and not args.no_structures:
                # vertical mechanism figure: the number just left of its title's
                # first line (2026-11-14), not on the picture
                t1_ = captions(k, True)[0].split("\n")[0]
                ty_ = fit_bottom - 0.08 if fit_bottom is not None else cy - size_in * 0.5 - 0.02
                sbx = cap_x - text_w(t1_, F(6.6 * TS), True) / 2 - BADGE / 144 - 0.15
                sby = ty_ - F(6.6 * TS) * 0.62 / 72
            for x_, y_, tag in ((bx, by, f"leader|{g}"), (sbx, sby, f"struct|{g}")):
                cv.tag(ax.scatter([x_], [y_], s=BADGE ** 2, c=folds.strong(k), zorder=10,
                                  linewidths=0), tag, f"badge {number}")
                # centred on the digit's own outline, not its text box (the box
                # includes room below the baseline, which sat the digits high)
                fs_ = F(5.2 * TS)
                ext_ = TextPath((0, 0), str(number), size=fs_,
                                prop=FontProperties(family=plt.rcParams["font.family"], weight="bold")).get_extents()
                cv.text(tag, x_ - (ext_.x0 + ext_.x1) / 2 / 72, y_ - (ext_.y0 + ext_.y1) / 2 / 72,
                        str(number), fontsize=fs_, color="white", ha="left", va="baseline",
                        fontweight="bold", zorder=11)
            sc.marker("leaders", k, bx, by, "o", folds.strong(k), BADGE ** 2, lw=0,
                      title=f"{folds.title(k)}: structure {number}")
            sc.text("leaders", k, bx, by, str(number), F(5.2 * TS), "white", weight="bold")
        else:
            lead = leader(a_mid, r_out + 0.05, r_lane, cx, cy,
                          struct_box(cx, cy, size_in, k in A and not args.no_structures), arc)
            # the leader in the colour of the fold it leads to
            cv.tag(ax.plot(*zip(*lead), color=folds.strong(k), lw=0.8, ls=(0, (1, 1.6)),
                           zorder=3)[0], f"leader|{g}", "leader to structure")
            ax.scatter([lead[0][0]], [lead[0][1]], s=4, c=folds.strong(k), zorder=3)
        pictured = k in A and not args.no_structures
        if pictured:
            png = os.path.join(cards.OUTDIR, view, f"{cards.slug(k)}.png")
            img = plt.imread(png)
            im = OffsetImage(img, zoom=size_in * 72 / img.shape[1])
            cv.tag(ax.add_artist(AnnotationBbox(im, (cx, cy), frameon=False, zorder=4)),
                   f"struct|{g}", f"structure {folds.title(k)}")
            if args.label_residues and view == "full":
                label_protein(k, g, cx, cy, size_in)
            if args.insets and view == "full" and k not in no_inset:
                if k in cell_of:                 # --fit-content: after the loop, at one shared size
                    pending_insets.append((k, g, cx, cy, size_in))
                else:
                    draw_inset(k, g, cx, cy, size_in)
            ty = cy - size_in * 0.5 - 0.02
            if fit_bottom is not None:
                ty = fit_bottom - 0.08
        else:
            ty = cy + 0.05              # the name alone, where the picture would be
        title, sub = captions(k, pictured)
        CS = TS if TS > 1 else 1.8 if few and args.panel else 1.0        # captions follow the larger pictures
        cv.text(f"label|{g}", cap_x, ty, title, fontsize=F(6.6 * CS), color=folds.strong(k),
                ha="center", va="top", fontweight="bold", zorder=9)
        cv.text(f"label|{g}", cap_x, (ty - (title.count("\n") + 1) * 1.25 * VFS / 72 - 0.08 if VFS else ty - 0.14 * CS), sub,
                fontsize=F(4.8 * CS), color=folds.strong(k) if args.residues else MUTED, ha="center",
                va="top", zorder=9, linespacing=1.3)
        if fit_bottom is None:
            cv.grow(cx - size_in / 2 - 0.3, cy - size_in / 2 - 0.65)
            cv.grow(cx + size_in / 2 + 0.3, cy + size_in / 2 + 0.1)
        else:                            # the drawn extent and the close up, not the square picture
            x0_, y0_, x1_, y1_ = content_box(k, view)
            cv.grow(cx - size_in / 2 + size_in * x0_ - 0.1, fit_bottom - 0.65)
            cv.grow(cx - size_in / 2 + size_in * x1_ + 0.1, cy + size_in / 2 - size_in * y0_ + 0.1)
        sc.structure(k, cx, cy, size_in,
                     None if args.panel else {"angle": a_mid, "r0": r_out + 0.05, "lane": r_lane},
                     title, sub, folds.strong(k), pictured)
        if args.panel:
            sc.structures[-1]["badge"] = number
            sc.structures[-1]["text_scale"] = CS
            sc.structures[-1]["badge_r"] = BADGE / 2
        sc.structures[-1]["contents"] = inside.get(k, [])
        sc.structures[-1]["leader"] = [list(p) for p in lead]
        if moved:
            sc.structures[-1]["view"] = view

    if pending_insets:
        # one size for all: the largest that clears the ring in every panel
        even_r.append(min(inset_at(*p_[:1], *p_[2:])[5] for p_ in pending_insets))
        for k_, g_, cx_, cy_, s_ in pending_insets:
            rx_, ry_, rr_, ix_, iy_, ir_ = inset_at(k_, cx_, cy_, s_)
            print(f"  inset {folds.type_name(k_):<32} magnifies {ir_ / rr_:.1f}x "
                  f"({2 * ir_ / (cell_w - 0.2):.0%} of the column, {'right' if ix_ > cx_ else 'left'})")
            draw_inset(k_, g_, cx_, cy_, s_)
            cv.grow(ix_ + ir_ + 0.05, iy_ + ir_ + 0.05)

    # --- vertical mechanism groups: the whole dataset as one band, from the
    # right of the rings up to a narrow tail at the top where the next pipeline
    # joins. Its height at the tree is the tree's height, so every record goes in.
    if left_side:
        y_top, y_bot = ang(slices[0][2]), ang(slices[-1][3])
        end_h = 1.0                                      # inches: the narrow tail
        fx0, flen, fstraight = r_out + 0.15, 6.0, 0.6   # start, taper, flat run at the end
        # Five nodes: the top edge, the flat tail, one cubic curve down to the
        # tree's bottom, and back up the left edge. Easy to adjust in Inkscape.
        xe_, yb_ = fx0 + flen + fstraight, y_top - end_h
        cx_ = fx0 + flen / 2                              # both handles halfway along
        verts_ = [(fx0, y_top), (xe_, y_top), (xe_, yb_), (fx0 + flen, yb_),
                  (cx_, yb_), (cx_, y_bot), (fx0, y_bot), (fx0, y_top)]
        codes_ = [MPath.MOVETO, MPath.LINETO, MPath.LINETO, MPath.LINETO,
                  MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
        s4_ = 3 * 0.4 ** 2 - 2 * 0.4 ** 3                # about where the curve sits under the label
        y_mid = (y_top + y_bot + (yb_ - y_bot) * s4_) / 2
        FLOW = "#f0f2f8"                                 # lavender, as in manuscript 2a
        if args.band:
            cv.tag(ax.add_patch(PathPatch(MPath(verts_, codes_), facecolor=FLOW, edgecolor="none",
                                          zorder=0)), "flow", "plasticome query")
            sc.poly("query flow", "", verts_[:4] + [verts_[4] + verts_[5] + verts_[6]], FLOW,
                    title="plasticome, the input query")
        n_ = "609 sequences" if SET == "all" else "411 centroids"
        fs_ = F(8 * TS)
        if args.band:     # right of the tail, where the band leaves the figure
            lx_, ly_, q_ = xe_ + 0.35, y_top - end_h / 2, "Input Query"
        else:             # above the tree at its left edge, as placed by hand (2026-11-14)
            lx_, ly_, q_ = bar_x + 0.25, y_top + 0.35 + F(5.5 * TS) * 1.25 / 72, "Query"
        lines_ = (("PETadex-PAZy", ly_ + 0.05, fs_, INK, "bold"),
                  (f"{n_} ({q_})", ly_ - 0.05, F(5.5 * TS), MUTED, "normal"))
        for s_txt, y_, size_, col_, wt_ in lines_:
            va_ = "bottom" if wt_ == "bold" else "top"
            t_txt = retext(f"flow|{wt_}", s_txt)    # ax.text: cv.text would widen the box now
            cv.tag(ax.text(lx_, y_, t_txt, fontsize=size_, color=col_, ha="left", va=va_,
                           fontweight=wt_, zorder=9), "flow", t_txt)
            sc.text("query flow", "", lx_, y_, s_txt, size_, col_, ha="left", va=va_, weight=wt_,
                    tid=f"flow|{wt_}")
        if args.band:
            cv.grow(fx0 + flen + fstraight, y_top)
        # the label's reach is added after the legends are placed: it sits above
        # them, so it must not push them right (the legend column starts at the box edge)
        label_reach_ = (lx_ + max(text_w(t_, z_, w_ == "bold") for t_, _, z_, _, w_ in lines_) + 0.2,
                        y_top)

    clashes = check_layout(sc)
    print("layout check: " + (", ".join(f"{n} {what}" for what, n in clashes.items())
                              if clashes else "no crossings or overlaps"))

    for sty, (col, lw) in LINES.items():
        d7.STYLE[sty] = (col, lw * args.line_scale, d7.STYLE[sty][2])
    cv.flush()

    # --- legends ---------------------------------------------------------------
    box = cv.box
    xl, yt = box[2] + (0.47 if TS > 1 else 0.3), box[3] - 0.2   # large text: equal gaps either side of the structures
    if left_side:
        cv.grow(*label_reach_)

    LS = TS if TS > 1 else 1.6 if few else 1.0                  # legend size, relative to the earlier figures

    def entry(h):
        """A matplotlib legend handle as a scene legend entry."""
        mk, ls = h.get_marker(), h.get_linestyle()
        e = {"label": h.get_label(), "color": ring_scene.css(h.get_color())}
        if mk == "s" and ls in ("None", "none", ""):
            fc = h.get_markerfacecolor()
            e.update(kind="square", color=ring_scene.css(fc) if fc != "none" else "white",
                     edge=ring_scene.css(h.get_markeredgecolor())
                     if h.get_markeredgecolor() != "none" else None)
        elif ls in ("None", "none", ""):
            hollow = h.get_markerfacecolor() == "none"
            e.update(kind="marker", shape=mk, size=h.get_markersize() ** 2, hollow=hollow,
                     color=ring_scene.css(h.get_markeredgecolor() if hollow
                                          else h.get_markerfacecolor()))
        else:
            e.update(kind="line", lw=h.get_linewidth(),
                     dash=[2.2, 1.6] if "--" in str(ls) or ls == (0, (2.2, 1.6)) or
                     h.get_label().startswith("backbone") else
                     [1, 1.6] if ls in (":", "dotted") else None)
        return e

    def legend(key, handles, title, y, ncol=1):
        off = layout.get("legends", {}).get(key, {})
        sc.legend(key, xl, y, title, [entry(h) for h in handles], ncol, scale=LS)
        title = retext(f"legend|{key}", title)
        lg = ax.legend(handles=handles, title=title, loc="upper left",
                       bbox_to_anchor=(xl + off.get("dx", 0), y + off.get("dy", 0)),
                       bbox_transform=ax.transData, frameon=False, fontsize=F(7 * LS), ncol=ncol,
                       handletextpad=0.5, columnspacing=1.2, labelspacing=0.5,
                       borderaxespad=0, alignment="left", handlelength=1.8)
        lg.get_title().set_fontsize(F(7.5 * LS))
        lg.get_title().set_color(INK)
        lg.get_title().set_fontweight("bold")
        ax.legend_ = None
        ax.add_artist(lg)
        cv.tag(lg, f"legend|{key}", f"{title} legend")

    sq = lambda c, lab, **kw: Line2D([], [], marker="s", ls="none", ms=7 * LS, mfc=c,
                                     mec=kw.get("mec", "none"), label=lab)
    dot = lambda c, lab, ms=6, **kw: Line2D([], [], marker="o", ls="none", ms=ms * LS, mfc=c,
                                            mec=kw.get("mec", "none"), mew=kw.get("mew", 0.5),
                                            label=lab)
    bar_ = lambda c, lab, lw, **kw: Line2D([], [], color=c, lw=lw * LS, label=lab, **kw)

    def height(handles, title, ncol):
        """A legend's height in inches (1 data unit = 1 inch), from its rows,
        so the blocks can be stacked with one even gap."""
        rows = [handles[n:n + ncol] for n in range(0, len(handles), ncol)]
        body = sum((max(h.get_label().count("\n") + 1 for h in r) * 1.2 + 0.5) * F(7 * LS) / 72
                   for r in rows)
        return (title.count("\n") + 1) * F(7.5 * LS) * 1.3 / 72 + body

    blocks = [
        ("lines", "Branches", [
            bar_(LINES["backbone"][0], "backbone (not to scale)", 1.0, ls=(0, (2.2, 1.6))),
            bar_(LINES["solid"][0], "support \u2265 70", 1.0),
            bar_(LINES["weak"][0], "support < 70", 1.0)], 1),
        ("plastic", "Plastic (rings, tip colour)",
         [sq(RING_COLOR[p], p) for p in PLASTICS] + [sq(SUBSTRATE_COLOR["Combi."], "mixed")], 2 if TS > 1 else 4),
        ("size", "Records per tip" if SET == "all" else "Records in cluster",
         [Line2D([], [], marker="o", ls="none", color="#9a9a9a",
                 ms=LS * math.sqrt(6 + 7 * (n - 1) ** 0.85 if n > 1 else 6), label=str(n))
          for n in ([1, 2, 5, 10] if SET == "all" else [1, 3, 10, 21])], 3),
    ]
    if TS > 1:           # the ring names across the top already key the plastic colours
        blocks = [b_ for b_ in blocks if b_[0] != "plastic"]
    if not args.no_structures and cards.OUTDIR.startswith(grey.OUTDIR):    # 17b renders, either view
        blocks.append(("structure", "Structure annotation", [
            bar_(MUTED, "mechanism (group colour)", 2.6),
            sq(grey.CARTOON, "protein", mec="#9aa0a8"),
            dot(grey.SUBSTRATE, "plastic", 8, mec=grey.OUTLINE, mew=0.8)], 1))
    elif not args.no_structures:
        blocks.append(("structure", "Structure annotation", [
            bar_(cards.CATALYTIC, "catalytic", 2.6), bar_(cards.SECONDARY, "secondary", 1.6),
            dot(cards.SUBSTRATE, "plastic", 8), bar_(cards.HEME, "heme", 1.8),
            dot(cards.METAL, "metal ion", 6, mec="#7d838a", mew=0.4)], 2))
    if left_side:        # vertical mechanism figure: no legend, the caption explains (2026-11-14)
        blocks = []
    GAP = 0.3
    y = yt
    if left_side:        # stacked in the bottom right corner, clear of the flow's tail at the top
        GAP = 0.15 * LS
        used = sum(height(h_, t_, n_) for _, t_, h_, n_ in blocks) + GAP * (len(blocks) - 1)
        y = box[1] + 0.45 + used
    elif TS > 1:         # the blocks spread down the whole height, on the empty side
        used = sum(height(h_, t_, n_) for _, t_, h_, n_ in blocks) + 0.2 * LS * 4
        GAP = max(GAP, (yt - box[1] - used) / len(blocks))
    for key, title, handles, ncol in blocks:
        legend(key, handles, title, y, ncol=ncol)
        y -= height(handles, title, ncol) + GAP
    bar = 0.5
    if left_side:        # the scale bar under the tree's bottom left, not beside the plastics
        xl, y = INNER, ang(slices[-1][3]) - 0.86
    cv.tag(ax.add_line(Line2D([xl, xl + bar * scale], [y, y], color=INK, lw=1.0 * args.line_scale)),
           "legend|scale", "scale bar")
    cv.text("legend|scale", xl, y + 0.06, f"{bar} substitutions/site", fontsize=F(7 * LS),
            color=INK, ha="left", va="bottom")
    if not left_side:
        cv.grow(xl + (1.8 * LS * F(7 * LS) / (7 * LS) if TS > 1 else 3.2 * LS + 0.8), y - 0.2)   # the widest legend line
    else:
        cv.grow(xl, y - 0.2)

    d8.OUT = OUT
    box_ = list(cv.box)
    w, h = d8.save(cv, STEM, vector_dpi=600)   # structures sharp in PDF/SVG

    # --- the editor ------------------------------------------------------------
    images = {}
    for k in {s_["key"] for s_ in sc.structures if s_["image"]}:
        images[k] = {v: ring_scene.png_data_uri(os.path.join(cards.OUTDIR, v, f"{cards.slug(k)}.png"))
                     for v in ("full", "zoom")}
    sc.write_html(os.path.join(OUT, f"{STEM}.html"),
                  f"Plasticome ring tree ({'all 609 records' if SET == 'all' else '411 centroids'})",
                  [b + d for b, d in zip(box_, (-0.25, -0.25, 0.25, 0.25))], images,
                  {"x": xl, "y": y, "length": bar * scale, "label": f"{bar} substitutions per site",
                   "scale": LS})

    # --- the order, for the heatmap ------------------------------------------
    ordered_tips.sort(key=lambda x: x[1])
    write_tsv(os.path.join(OUT, "ring_order.tsv"),
              ["tip", "order", "fold_group", "fold_key", "fold", "component", "enzyme_name",
               "position", "hue_hex", "size"],
              [[t["tip"], i, t["fold_group"], keyof(t), folds.name(keyof(t)), t["component"],
                t["enzyme_name"], f"{frac:.5f}",
                "#%02x%02x%02x" % tuple(round(c * 255) for c in rainbow(frac)), size[t["tip"]]]
               for i, (t, frac, a) in enumerate(ordered_tips)])
    c = collections.Counter(v[0] for v in cons.values())
    print(f"figure_ring: {len(ordered_tips)} tips, {len(slices)} folds, canvas {w:.1f} x {h:.1f} in")
    print("catalytic ring: " + ", ".join(f"{k} {v}" for k, v in c.most_common()))


if __name__ == "__main__":
    main()
