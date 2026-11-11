#!/usr/bin/env python3
"""
Build panel a-i: sequence similarity network of the plasticome v2 union set.

Encoding (matches the published PETadex/Logan workflow panel a-i):
  node colour  = substrate
  node shape   = Bacterial (diamond) / Eukaryotic (hexagon) / Metagenomic (circle)
  edge width   = aa identity %, scaled over 40-100
  grey blob    = connected component at the family threshold (30% id), labelled
                 with its dominant CATH domain family

Two thresholds do two different jobs:
  COMPONENT_ID  defines the grey blobs / family membership. At 30% the
                components become family-pure (94.3% -> 99.8% purity); below
                that the alpha/beta hydrolase component starts absorbing
                other families.
  DRAW_ID       the floor for an edge to be *drawn*, so the 40-100 width
                legend stays honest. The drawn set is thinned to a backbone
                (maximum spanning forest + top-K strongest edges per node)
                because a flat 40% cutoff on 609 nodes is 37,097 edges.

Usage:  python build_ssn.py [--component-id 30] [--draw-id 40] [--topk 3]
"""

import argparse
import collections
import csv
import math
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
import networkx as nx
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "plasticome.v2 - plasticome-v1-v260701-union-annotated (1).csv")
HITS = os.path.join(HERE, "work", "allvall.tsv")
SEEDS = os.path.join(HERE, "..", "..", "source-data", "plasticome-curated-seeds",
                     "plasticome-curated-seeds.tsv")
OUT = os.path.join(HERE, "out")

# --- palette -----------------------------------------------------------------
# Sampled from the published panel, then re-stepped on the pairs that failed
# scripts/validate_palette.js: PBAT was dE 2.0 from PET under deuteranopia
# (PET is 56% of nodes), PHA/PUR and PE/PET were below the normal-vision floor.
SUBSTRATE_COLOR = {
    "PET":    "#275ca0",
    "PE":     "#e8559e",
    "PA":     "#40b6c5",
    "PUR":    "#f0c95e",
    "PLA":    "#cc490b",
    "PBAT":   "#9a6fd4",
    "PBS":    "#23873d",
    "PHA":    "#8fbf5f",
    "PC":     "#8c6d46",
    "NR":     "#959595",
    "Combi.": "#262626",
}
SUBSTRATE_ORDER = ["PET", "PE", "PA", "PUR", "PLA", "PBAT", "PBS", "PHA", "PC", "NR", "Combi."]

SHAPE = {"Bacteria": "D", "Archaea": "D", "Eukaryota": "h", "Metagenomic": "o", "NA": "o"}
SHAPE_LABEL = [("Bacterial", "D"), ("Eukaryotic", "h"), ("Metagenomic", "o")]

INK = "#1a1a1a"
MUTED = "#6b6b6b"
BLOB = "#e8e8e6"
SURFACE = "#fcfcfb"
EDGE_COLOR = "#b4b4b0"
NODE_STROKE = "#3a3a3a"   # relief for the pale fills (PUR/PHA/PA contrast WARN)

NN_SPACING = 1.0    # target median nearest-neighbour node spacing, in data units
BLOB_PAD = 0.95     # grey halo radius around each node, in data units


def substrate_of(raw):
    s = (raw or "").strip()
    if s in ("", "NA"):
        return "NR"
    parts = [p.strip() for p in s.split(";") if p.strip()]
    if len(parts) > 1:
        return "Combi."
    return parts[0] if parts[0] in SUBSTRATE_COLOR else "Combi."


def family_of(domain):
    d = (domain or "").split(";")[0].strip().lower()
    if "alpha/beta" in d:
        return "A/B hydrolase"
    if "l-amino peptidase" in d:
        return "L-amino peptidase\nD-ALA esterase"
    if "amidase" in d:
        return "Amidase"
    if "trypsin" in d:
        return "Trypsin-like"
    if "peptidase" in d or "serralysin" in d:
        return "Peptidase"
    if "cytochrome" in d:
        return "Cytochrome C-like"
    if "oxygenase" in d or "monooxygenase" in d:
        return "Oxygenase-like"
    if "lactamase" in d:
        return "Beta-lactamase"
    if "arylesterase" in d:
        return "Arylesterase"
    if "cyclase" in d:
        return "Cyclase-like"
    if "uridase" in d or "amidohydrolase" in d:
        return "Uridase-like"
    if d in ("", "na", "-"):
        return "Unassigned"
    return "Other"


def shape_of(row):
    if (row.get("metagenomic") or "").strip().lower() == "yes":
        return "o"
    return SHAPE.get((row.get("superkingdom") or "NA").strip(), "o")


def load():
    fh = open(CSV, newline="")
    next(fh)  # provenance line above the header
    rows = {r["plasticome_id"]: r for r in csv.DictReader(fh)}
    fh.close()
    return rows


def slug(text):
    """Lowercase ascii slug for SVG ids."""
    out = "".join(c.lower() if c.isalnum() else "-" for c in text)
    while "--" in out:
        out = out.replace("--", "-")
    return out.strip("-")


def seed_nodes(rows):
    """Curated seed enzymes -> plasticome_id, matched on accession then sequence."""
    by_acc = {r["accession"].strip(): pid for pid, r in rows.items()}
    by_seq = {r["aa_sequence"].strip().upper(): pid for pid, r in rows.items()}
    found = {}
    with open(SEEDS, newline="") as fh:
        for s in csv.DictReader(fh, delimiter="\t"):
            acc = (s.get("accession") or "").strip()
            seq = (s.get("aa_sequence") or "").strip().upper()
            pid = by_acc.get(acc) or by_seq.get(seq)
            if pid is not None:
                found[s["enzyme_name"].strip()] = pid
    return found


def load_pairs(cov_min):
    """Best (highest-identity) alignment per unordered pair, coverage-filtered."""
    best = {}
    with open(HITS) as fh:
        for line in fh:
            q, s, pid, ln, ql, sl, bs, ev = line.rstrip("\n").split("\t")
            if q == s:
                continue
            pid = float(pid)
            cov = int(ln) / min(int(ql), int(sl))
            if cov < cov_min:
                continue
            k = (q, s) if q < s else (s, q)
            if k not in best or pid > best[k]:
                best[k] = pid
    return best


def max_spanning(G):
    """Maximum spanning edges, reproducibly.

    Many pairs here share an identical pident, and Kruskal breaks those ties in
    graph iteration order, which depends on the per process string hash seed.
    Two runs then pick different (equal weight) spanning trees and the whole
    layout shifts. Perturbing each weight by a canonical 1e-9 rank makes the
    maximum spanning tree unique, far below the 0.1 resolution of pident.
    """
    ordered = sorted(((a, b) if int(a) < int(b) else (b, a))
                     for a, b in G.edges())
    T = nx.Graph()
    T.add_nodes_from(sorted(G.nodes(), key=int))
    for i, (a, b) in enumerate(ordered):
        T.add_edge(a, b, _w=G[a][b]["pident"] + i * 1e-9)
    return sorted(nx.maximum_spanning_edges(T, weight="_w", data=False),
                  key=lambda e: (int(e[0]), int(e[1])))


def backbone(G, topk):
    """Maximum spanning forest (keeps each component in one piece) plus each
    node's topk strongest partners. Every kept edge is a real alignment."""
    keep = {(a, b) if a < b else (b, a) for a, b in max_spanning(G)}
    for n in G.nodes():
        nbrs = sorted(G[n].items(), key=lambda kv: -kv[1]["pident"])[:topk]
        for m, _ in nbrs:
            keep.add((n, m) if n < m else (m, n))
    return keep


def shelf_pack(boxes, aspect=1.15, slots=420):
    """Bottom-left skyline packing of (key, w, h) boxes, largest first.

    A plain shelf packer puts every box on a row whose height is the tallest
    box in it, so one tall block (the A/B hydrolase section, or its 329 node
    component) leaves a band of dead space beside everything shorter. The
    skyline tracks the filled height across x and drops each box at the lowest
    position it fits, which tucks the small sections alongside the tall one.

    Origin is top left and y grows downward; the caller flips y for the plot.
    Returns ({key: (x, y)}, total_w, total_h).
    """
    order = sorted(boxes, key=lambda b: (-(b[1] * b[2]), str(b[0])))
    total = sum(w * h for _, w, h in order)
    W = max(math.sqrt(total) * aspect, max(w for _, w, _ in order) * 1.02)
    sx = W / slots
    heights = np.zeros(slots)
    placed = {}
    for k, w, h in order:
        nw = min(max(int(math.ceil(w / sx)), 1), slots)
        best_y, best_i = None, 0
        for i in range(slots - nw + 1):
            y = float(heights[i:i + nw].max())
            if best_y is None or y < best_y - 1e-9:
                best_y, best_i = y, i
        placed[k] = (best_i * sx, best_y)
        heights[best_i:best_i + nw] = best_y + h
    return placed, W, float(heights.max())


def regroup_svg(path, groups):
    """Rewrite the matplotlib SVG so each component is one named <g>.

    Matplotlib emits every blob, edge collection and node collection as its own
    top level element, so in Inkscape a component cannot be selected or moved as
    a unit. Each artist is drawn with a gid, and this pass moves the artists
    belonging to one component into a single group carrying an inkscape:label,
    which is the name Inkscape shows in the Objects panel.
    """
    import xml.etree.ElementTree as ET

    SVG = "http://www.w3.org/2000/svg"
    INK = "http://www.inkscape.org/namespaces/inkscape"
    ET.register_namespace("", SVG)
    ET.register_namespace("inkscape", INK)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")

    tree = ET.parse(path)
    root = tree.getroot()
    parent_of = {c: p for p in root.iter() for c in p}
    by_id = {el.get("id"): el for el in root.iter() if el.get("id")}

    moved = 0
    for gid, label, members in groups:
        present = [by_id[m] for m in members if m in by_id]
        if not present:
            continue
        anchor = present[0]
        host = parent_of[anchor]
        g = ET.Element(f"{{{SVG}}}g")
        g.set("id", gid)
        g.set(f"{{{INK}}}label", label)
        host.insert(list(host).index(anchor), g)
        for el in present:
            parent_of[el].remove(el)
            g.append(el)
            moved += 1

    tree.write(path, encoding="utf-8", xml_declaration=True)
    return moved


def bridge_points(pts, step):
    """Points sampled along a Euclidean minimum spanning tree of `pts`.

    The halo has to read as ONE region per component, because a break in it
    reads as two components. A plain union of disks around the nodes splits
    whenever two members sit more than 2*pad apart, which happens in 5 of the
    35 multi-node components here. Adding the spanning tree guarantees the
    union is connected, and the tree's edges are the shortest bridges available
    so the necks stay short.
    """
    from scipy.sparse.csgraph import minimum_spanning_tree
    if len(pts) < 2:
        return np.empty((0, 2))
    d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    mst = minimum_spanning_tree(d).tocoo()
    out = []
    for a, b in zip(mst.row, mst.col):
        p, q = pts[a], pts[b]
        n = max(int(np.linalg.norm(q - p) / step) + 1, 2)
        out.append(p + (q - p) * np.linspace(0, 1, n)[:, None])
    return np.vstack(out) if out else np.empty((0, 2))


def draw_blob(ax, points, pad, color, zorder=0, res=0.10, gid=None,
              bridge_pad=None):
    """Grey halo that hugs the cluster and stays in one piece.

    Two radii: `pad` around every node, and a narrower `bridge_pad` along the
    spanning tree that links them, so the region is guaranteed connected
    without inflating into a blob that swallows the empty space. A convex hull
    was tried first and encloses huge empty regions on sparse components.
    """
    from scipy.spatial import cKDTree
    pts = np.asarray(points, dtype=float)
    if bridge_pad is None:
        bridge_pad = pad * 0.42
    lo = pts.min(axis=0) - pad * 1.6
    hi = pts.max(axis=0) + pad * 1.6
    nx_ = min(max(int((hi[0] - lo[0]) / res) + 2, 12), 1400)
    ny_ = min(max(int((hi[1] - lo[1]) / res) + 2, 12), 1400)
    gx = np.linspace(lo[0], hi[0], nx_)
    gy = np.linspace(lo[1], hi[1], ny_)
    X, Y = np.meshgrid(gx, gy)
    G = np.column_stack([X.ravel(), Y.ravel()])

    f = cKDTree(pts).query(G)[0] - pad
    br = bridge_points(pts, step=res * 0.5)
    if len(br):
        f = np.minimum(f, cKDTree(br).query(G)[0] - bridge_pad)
    f = f.reshape(X.shape)

    cs = ax.contourf(X, Y, f, levels=[f.min() - 1.0, 0.0], colors=[color],
                     zorder=zorder, antialiased=True)
    if gid:
        cs.set_gid(gid)
    return cs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--component-id", type=float, default=30.0)
    ap.add_argument("--draw-id", type=float, default=40.0)
    ap.add_argument("--coverage", type=float, default=0.5)
    ap.add_argument("--topk", type=int, default=3)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--tag", default="")
    ap.add_argument("--inner-aspect", type=float, default=1.00)
    ap.add_argument("--outer-aspect", type=float, default=1.10)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    rows = load()
    pairs = load_pairs(args.coverage)
    ids = list(rows)

    # --- components at the family threshold ---------------------------------
    C = nx.Graph()
    C.add_nodes_from(ids)
    for (a, b), p in pairs.items():
        if p >= args.component_id:
            C.add_edge(a, b, pident=p)
    comp_of = {}
    comps = []
    # Sort by size, then by smallest member id. Without the tiebreak, equal-sized
    # components come back in set-iteration order, which changes with the per
    # process string hash seed, and the figure is not reproducible run to run.
    cc_sorted = sorted((sorted(c, key=int) for c in nx.connected_components(C)),
                       key=lambda c: (-len(c), int(c[0])))
    for i, cc in enumerate(cc_sorted):
        comps.append(cc)
        for n in cc:
            comp_of[n] = i

    comp_family = {}
    for i, cc in enumerate(comps):
        cnt = collections.Counter(family_of(rows[n]["domain"]) for n in cc)
        comp_family[i] = cnt.most_common(1)[0][0]

    # --- drawn edges ---------------------------------------------------------
    D = nx.Graph()
    D.add_nodes_from(ids)
    for (a, b), p in pairs.items():
        if p >= args.draw_id:
            D.add_edge(a, b, pident=p)
    # sorted, not the raw set: set iteration order of string tuples varies with
    # the per process hash seed, which perturbs the force layout
    kept = sorted(backbone(D, args.topk), key=lambda e: (int(e[0]), int(e[1])))
    kept_set = set(kept)
    E = nx.Graph()
    E.add_nodes_from(ids)
    for a, b in kept:
        E.add_edge(a, b, pident=D[a][b]["pident"])

    # --- layout, per component, then packed ---------------------------------
    # Positions come from the COMPONENT graph (30%), not the drawn graph (40%).
    # A node whose only similarity is in the 30-40% band has no drawn edge; laid
    # out on the drawn graph it is a disconnected point and force layout flings
    # it out of its own family blob. The 30% subgraph is connected by
    # construction, so every node is placed inside its family.
    from scipy.spatial import cKDTree
    layouts = {}
    for i, cc in enumerate(comps):
        n_cc = len(cc)
        if n_cc == 1:
            layouts[i] = {cc[0]: (0.0, 0.0)}
            continue
        # Layout graph: the drawn backbone (real substructure) plus a maximum
        # spanning tree of the component's 30% graph (guarantees one piece, adds
        # the minimum extra structure). Using the full 30% subgraph instead
        # gives a featureless disc, because it is near-complete inside a family.
        L = nx.Graph()
        L.add_nodes_from(cc)
        for a, b in E.subgraph(cc).edges():
            L.add_edge(a, b, pident=E[a][b]["pident"])
        for a, b in max_spanning(C.subgraph(cc)):
            if not L.has_edge(a, b):
                L.add_edge(a, b, pident=C[a][b]["pident"])

        # Cubic identity weighting: a 90% pair pulls ~30x harder than a 40% pair,
        # so near-identical groups condense into visible subclusters while the
        # component as a whole stays a compact disc (r50 ~= sqrt(n/pi)).
        for a, b, d in L.edges(data=True):
            d["w"] = ((d["pident"] - args.component_id) /
                      (100.0 - args.component_id)) ** 3 * 10.0 + 0.05

        if n_cc >= 20:
            p = nx.forceatlas2_layout(L, max_iter=1200, jitter_tolerance=1.0,
                                      scaling_ratio=2.0, gravity=2.0,
                                      strong_gravity=True, weight="w",
                                      linlog=True, seed=args.seed)
        else:
            p = nx.spring_layout(L, k=2.5 / math.sqrt(n_cc), iterations=600,
                                 seed=args.seed, weight=None)
        arr = np.array([p[n] for n in cc], dtype=float)
        arr -= arr.mean(axis=0)
        # Normalise every component to the same node density: median
        # nearest-neighbour spacing becomes NN_SPACING data units, so one marker
        # size and one blob pad work everywhere.
        if n_cc > 1:
            nn = cKDTree(arr).query(arr, k=2)[0][:, 1]
            med = float(np.median(nn)) or 1.0
            arr *= NN_SPACING / med
        layouts[i] = {n: tuple(arr[j]) for j, n in enumerate(cc)}

    # --- sections: components grouped by catalytic family --------------------
    # Every component of a family is packed into one block, and the blocks are
    # packed on the page, so the reader sees all the A/B hydrolase components
    # together, all the Amidase components together, and so on.
    # Gap must exceed two blob radii, or neighbouring halos touch and two
    # components read as one.
    GAP = 2 * BLOB_PAD + 0.9
    SECTION_GAP = GAP * 2.2
    LABEL_H = 2.6

    sections = collections.defaultdict(list)
    for i in range(len(comps)):
        sections[comp_family[i]].append(i)
    # largest family first, ties by name so the page order is reproducible
    sec_order = sorted(sections, key=lambda f: (-sum(len(comps[i])
                                                     for i in sections[f]), f))

    comp_extent = {}
    for i, lay in layouts.items():
        xs = [p[0] for p in lay.values()]
        ys = [p[1] for p in lay.values()]
        comp_extent[i] = (min(xs), min(ys),
                          (max(xs) - min(xs)), (max(ys) - min(ys)))

    inner, sec_box = {}, {}
    for fam in sec_order:
        boxes = [(i, comp_extent[i][2] + GAP, comp_extent[i][3] + GAP)
                 for i in sections[fam]]
        placed, W, H = shelf_pack(boxes, aspect=args.inner_aspect)
        inner[fam] = placed
        sec_box[fam] = (W, H)

    # The section box must be wide enough for its own heading, or a short
    # section's label runs into the neighbour placed beside it.
    def label_width(fam):
        head = fam.replace("\n", " ")
        n_comp = len(sections[fam])
        sub = f"{n_comp} component{'s' if n_comp != 1 else ''}, ? seqs"
        return max(len(head) * 0.70, len(sub) * 0.50) + 1.0

    sec_placed, _, _ = shelf_pack(
        [(fam, max(sec_box[fam][0], label_width(fam)) + SECTION_GAP,
          sec_box[fam][1] + LABEL_H + SECTION_GAP) for fam in sec_order],
        aspect=args.outer_aspect)

    pos = {}
    section_anchor = {}
    for fam in sec_order:
        sx, sy = sec_placed[fam]
        section_anchor[fam] = (sx + SECTION_GAP / 2,
                               -(sy + SECTION_GAP / 2 + LABEL_H * 0.55))
        for i in sections[fam]:
            cx, cy = inner[fam][i]
            x0, y0, _, _ = comp_extent[i]
            dx = sx + SECTION_GAP / 2 + cx + GAP / 2 - x0
            dy = sy + SECTION_GAP / 2 + LABEL_H + cy + GAP / 2 - y0
            for n, (x, y) in layouts[i].items():
                # y is flipped: shelf packing grows downward, the plot grows up
                pos[n] = (x + dx, -(y + dy))

    # --- draw ----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13.5, 11.0))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # Every artist is tagged with the component it belongs to, so the SVG can be
    # regrouped afterwards into one movable, named group per component.
    from matplotlib.collections import LineCollection
    gids = collections.defaultdict(list)
    section_label_gids = {}

    def tag(artist, i, suffix):
        g = f"plc-{i}-{suffix}"
        artist.set_gid(g)
        gids[i].append(g)
        return artist

    for i, cc in enumerate(comps):
        if len(cc) < 2:
            continue
        draw_blob(ax, [pos[n] for n in cc], pad=BLOB_PAD, color=BLOB, zorder=0,
                  gid=f"plc-{i}-blob")
        gids[i].append(f"plc-{i}-blob")

    # A drawn edge is >= 40% identity, so both ends are always in the same 30%
    # component; edges group cleanly by component.
    edges_by_comp = collections.defaultdict(list)
    for a, b, d in E.edges(data=True):
        edges_by_comp[comp_of[a]].append((a, b, d["pident"]))
    for i, es in edges_by_comp.items():
        segs = [[pos[a], pos[b]] for a, b, _ in es]
        lws = [0.28 + 1.25 * min(max((p - 40.0) / 60.0, 0.0), 1.0)
               for _, _, p in es]
        lc = LineCollection(segs, colors=EDGE_COLOR, linewidths=lws,
                            capstyle="round", zorder=1, alpha=0.55)
        ax.add_collection(tag(lc, i, "edges"))

    for i, cc in enumerate(comps):
        by_style = collections.defaultdict(list)
        for n in cc:
            by_style[(shape_of(rows[n]),
                      SUBSTRATE_COLOR[substrate_of(rows[n]["substrate"])])].append(n)
        for k, ((mk, col), members) in enumerate(sorted(by_style.items())):
            sc = ax.scatter([pos[n][0] for n in members],
                            [pos[n][1] for n in members],
                            marker=mk, s=46, c=col, edgecolors=NODE_STROKE,
                            linewidths=0.5, zorder=3)
            tag(sc, i, f"nodes{k}")

    # One label per SECTION, drawn in the space reserved above each family
    # block. Labelling every component repeated "A/B hydrolase" 30 times, and
    # the section layout already groups them, so the family is named once.
    all_pts = np.array([pos[n] for n in ids])
    node_tree = cKDTree(all_pts)
    node_comp = np.array([comp_of[n] for n in ids])
    placed_boxes = []

    def foreign_blob_hits(box, own):
        """Nodes of OTHER components whose blob halo the box would sit on.

        A label overlapping a foreign halo means paint order decides which
        wins, and the SVG is regrouped by component, which changes that order.
        Keeping labels off foreign halos makes the grouping safe."""
        cxb, cyb = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
        r = max(box[2] - box[0], box[3] - box[1]) / 2 + BLOB_PAD + 0.2
        hits = 0
        for j in node_tree.query_ball_point([cxb, cyb], r):
            if node_comp[j] == own:
                continue
            px_, py_ = all_pts[j]
            dx = max(box[0] - px_, 0, px_ - box[2])
            dy = max(box[1] - py_, 0, py_ - box[3])
            if dx * dx + dy * dy < BLOB_PAD ** 2:
                hits += 1
        return hits

    for fam in sec_order:
        lx, ly = section_anchor[fam]
        n_comp = len(sections[fam])
        n_seq = sum(len(comps[i]) for i in sections[fam])
        txt = f"{fam.replace(chr(10), ' ')}"
        sub = f"{n_comp} component{'s' if n_comp != 1 else ''}, {n_seq} seqs"
        t = ax.text(lx, ly, txt, ha="left", va="center", fontsize=12.5,
                    color=INK, zorder=6)
        t.set_gid(f"section-label-{slug(fam)}")
        t2 = ax.text(lx, ly - 1.15, sub, ha="left", va="center", fontsize=8.5,
                     color=MUTED, zorder=6)
        t2.set_gid(f"section-sub-{slug(fam)}")
        section_label_gids[fam] = [f"section-label-{slug(fam)}",
                                   f"section-sub-{slug(fam)}"]
        w = max(len(txt), len(sub)) * 0.62
        placed_boxes.append((lx, ly - 1.8, lx + w, ly + 0.9))

    # Curated seed enzymes, called out by name with a leader line. Same
    # anchor-scoring as the family labels, so a callout does not land on a node
    # or on another label.
    seed_ids = seed_nodes(rows)
    for name, n in sorted(seed_ids.items(), key=lambda kv: kv[0].lower()):
        px, py = pos[n]
        w = len(name) * 0.60
        h = 1.15
        best = None
        for rad in (2.4, 3.6, 5.2, 7.0):
            for ang in range(0, 360, 30):
                a = math.radians(ang)
                lx, ly = px + rad * math.cos(a), py + rad * math.sin(a)
                ha = "left" if math.cos(a) >= 0 else "right"
                bx0 = lx if ha == "left" else lx - w
                box = (bx0, ly - h / 2, bx0 + w, ly + h / 2)
                idx = node_tree.query_ball_point([(box[0] + box[2]) / 2, ly],
                                                 max(w, h) / 2 + 0.4)
                covered = sum(1 for j in idx
                              if box[0] <= all_pts[j][0] <= box[2]
                              and box[1] <= all_pts[j][1] <= box[3])
                overlap = sum(1 for o in placed_boxes
                              if not (box[2] < o[0] or box[0] > o[2] or
                                      box[3] < o[1] or box[1] > o[3]))
                score = (12 * covered + 400 * overlap + rad * 0.8
                         + 600 * foreign_blob_hits(box, comp_of[n]))
                if best is None or score < best[0]:
                    best = (score, lx, ly, ha, box)
            if best and best[0] < 1.0 + rad * 0.8:
                break
        _, lx, ly, ha, box = best
        placed_boxes.append(box)
        ci = comp_of[n]
        safe = "".join(c if c.isalnum() else "-" for c in name)
        tag(ax.annotate(name, xy=(px, py), xytext=(lx, ly), fontsize=8.5,
                        color=INK, zorder=7, ha=ha, va="center",
                        arrowprops=dict(arrowstyle="-", color=MUTED,
                                        linewidth=0.7, shrinkA=0.5,
                                        shrinkB=3.5)),
            ci, f"seed-{safe}")
        tag(ax.scatter([px], [py], marker=shape_of(rows[n]), s=58,
                       facecolors="none", edgecolors=INK, linewidths=1.1,
                       zorder=5),
            ci, f"seedring-{safe}")

    ax.set_aspect("equal")
    ax.axis("off")
    # Include the label boxes in the extent, or an outer label is clipped.
    xs = [all_pts[:, 0].min() - BLOB_PAD, all_pts[:, 0].max() + BLOB_PAD]
    ys = [all_pts[:, 1].min() - BLOB_PAD, all_pts[:, 1].max() + BLOB_PAD]
    for bx0, by0, bx1, by1 in placed_boxes:
        xs += [bx0, bx1]
        ys += [by0, by1]
    padx = (max(xs) - min(xs)) * 0.015
    pady = (max(ys) - min(ys)) * 0.015
    ax.set_xlim(min(xs) - padx, max(xs) + padx)
    ax.set_ylim(min(ys) - pady, max(ys) + pady)

    # Size the canvas to the content. The axes are equal-aspect, so a fixed
    # figure size letterboxes the drawing with broad empty margins whenever the
    # packed layout is not the same shape as the figure.
    data_w = (max(xs) + padx) - (min(xs) - padx)
    data_h = (max(ys) + pady) - (min(ys) - pady)
    ax_frac_w, ax_frac_h = 0.835, 0.98
    scale = math.sqrt(108.0 / (data_w * data_h))
    fig.set_size_inches(min(max(data_w * scale / ax_frac_w, 7.0), 26.0),
                        min(max(data_h * scale / ax_frac_h, 6.0), 26.0))

    present = [s for s in SUBSTRATE_ORDER
               if any(substrate_of(rows[n]["substrate"]) == s for n in ids)]
    h_sub = [Line2D([], [], marker="s", linestyle="none", markersize=8,
                    markerfacecolor=SUBSTRATE_COLOR[s], markeredgecolor=NODE_STROKE,
                    markeredgewidth=0.4, label=s) for s in present]
    l1 = ax.legend(handles=h_sub, loc="upper left", bbox_to_anchor=(1.005, 1.0),
                   ncol=2, frameon=False, fontsize=9, handletextpad=0.5,
                   columnspacing=1.1, labelspacing=0.55, title="substrate")
    l1.get_title().set_fontsize(9)
    l1.get_title().set_color(MUTED)
    ax.add_artist(l1)

    h_shape = [Line2D([], [], marker=mk, linestyle="none", markersize=8,
                      markerfacecolor="none", markeredgecolor=INK,
                      markeredgewidth=0.9, label=lab) for lab, mk in SHAPE_LABEL]
    l2 = ax.legend(handles=h_shape, loc="upper left", bbox_to_anchor=(1.005, 0.70),
                   frameon=False, fontsize=9, handletextpad=0.5, labelspacing=0.55,
                   title="origin")
    l2.get_title().set_fontsize(9)
    l2.get_title().set_color(MUTED)
    ax.add_artist(l2)

    h_edge = [Line2D([], [], color=EDGE_COLOR, linewidth=0.28 + 1.25 * (v - 40) / 60,
                     label=str(v)) for v in (100, 80, 60, 40)]
    l3 = ax.legend(handles=h_edge, loc="upper left", bbox_to_anchor=(1.005, 0.48),
                   frameon=False, fontsize=9, handletextpad=0.8, labelspacing=0.55,
                   title="aa id%")
    l3.get_title().set_fontsize(9)
    l3.get_title().set_color(MUTED)

    fig.subplots_adjust(left=0.01, right=0.845, top=0.99, bottom=0.01)

    tag = args.tag or f"c{int(args.component_id)}_d{int(args.draw_id)}_k{args.topk}"
    for ext in ("svg", "pdf", "png"):
        md = {"Date": None} if ext == "svg" else {}
        fig.savefig(os.path.join(OUT, f"figure_ai_{tag}.{ext}"), dpi=400,
                    facecolor=SURFACE, metadata=md)
    plt.close(fig)

    # One named, movable group per component in the SVG.
    svg_path = os.path.join(OUT, f"figure_ai_{tag}.svg")
    comp_groups = []
    for i, cc in enumerate(comps):
        fam = comp_family[i].replace("\n", " ")
        kind = "singleton" if len(cc) == 1 else "cluster"
        comp_groups.append((f"component-{i:02d}",
                            f"{i:02d} {fam} ({kind}, n={len(cc)})", gids[i]))
    n_moved = regroup_svg(svg_path, comp_groups)

    # Second pass nests the component groups inside one group per family, so a
    # whole section can be moved as a unit and a single component still can.
    sec_groups = []
    for fam in sec_order:
        members = section_label_gids.get(fam, []) + [
            f"component-{i:02d}" for i in sorted(sections[fam],
                                                 key=lambda j: -len(comps[j]))]
        n_seq = sum(len(comps[i]) for i in sections[fam])
        n_c = len(sections[fam])
        sec_groups.append((f"section-{slug(fam)}",
                           f"{fam.replace(chr(10), ' ')} "
                           f"({n_c} component{'s' if n_c != 1 else ''}, "
                           f"{n_seq} seqs)",
                           members))
    n_sec = regroup_svg(svg_path, sec_groups)
    used = sum((comp_extent[i][2] + GAP) * (comp_extent[i][3] + GAP)
               for i in range(len(comps)))
    print(f"        canvas {data_w:.0f}x{data_h:.0f} (aspect {data_w/data_h:.2f}), "
          f"fill {100*used/(data_w*data_h):.0f}%")
    print(f"        svg: {len(sec_groups)} section groups containing "
          f"{len(comp_groups)} component groups "
          f"({n_moved}+{n_sec} elements grouped)")

    # --- tables + graphml ----------------------------------------------------
    with open(os.path.join(OUT, f"ssn_nodes_{tag}.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["plasticome_id", "enzyme_name", "accession", "source",
                    "substrate_group", "substrate_raw", "family_group", "domain",
                    "superkingdom", "metagenomic", "shape", "color",
                    "component", "component_family", "component_size",
                    "degree_drawn", "x", "y"])
        for n in ids:
            r = rows[n]
            ci = comp_of[n]
            w.writerow([n, r["enzyme_name"], r["accession"], r["source"],
                        substrate_of(r["substrate"]), r["substrate"],
                        family_of(r["domain"]), r["domain"], r["superkingdom"],
                        r["metagenomic"], shape_of(r),
                        SUBSTRATE_COLOR[substrate_of(r["substrate"])],
                        ci, comp_family[ci], len(comps[ci]), E.degree(n),
                        round(pos[n][0], 4), round(pos[n][1], 4)])

    with open(os.path.join(OUT, f"ssn_edges_{tag}.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["source", "target", "pident", "drawn"])
        for (a, b), p in sorted(pairs.items()):
            if p >= args.draw_id:
                w.writerow([a, b, f"{p:.1f}", int((a, b) in kept_set)])

    GX = nx.Graph()
    for n in ids:
        r = rows[n]
        ci = comp_of[n]
        GX.add_node(n, enzyme_name=r["enzyme_name"], accession=r["accession"],
                    source=r["source"], substrate=substrate_of(r["substrate"]),
                    family=family_of(r["domain"]), superkingdom=r["superkingdom"],
                    metagenomic=r["metagenomic"], shape=shape_of(r),
                    color=SUBSTRATE_COLOR[substrate_of(r["substrate"])],
                    component=ci, component_family=comp_family[ci],
                    x=float(pos[n][0]), y=float(pos[n][1]))
    for a, b in kept:
        GX.add_edge(a, b, pident=float(D[a][b]["pident"]))
    nx.write_graphml(GX, os.path.join(OUT, f"ssn_{tag}.graphml"))

    sizes = sorted((len(c) for c in comps), reverse=True)
    print(f"[{tag}] nodes={len(ids)}  components@{args.component_id:.0f}%={len(comps)} "
          f"(largest {sizes[0]}, singletons {sum(1 for s in sizes if s == 1)})")
    print(f"        edges >= {args.draw_id:.0f}%: {D.number_of_edges()}  -> drawn backbone: "
          f"{E.number_of_edges()}  (E/N = {E.number_of_edges()/len(ids):.2f})")
    print(f"        wrote out/figure_ai_{tag}.{{svg,pdf,png}}, ssn_{tag}.graphml, "
          f"ssn_nodes_{tag}.tsv, ssn_edges_{tag}.tsv")


if __name__ == "__main__":
    main()
