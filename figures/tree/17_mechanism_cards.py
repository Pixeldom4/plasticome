#!/usr/bin/env python3
"""
Step 17: one card per fold showing the catalytic mechanism on its structure.

Input is the hand-annotated set in ../structures/annotations: 15
co-folded models (protein, the plastic oligomer as chain L, and any heme or
metal cofactors) and petadex_domain_summary.csv, which lists the catalytic and
secondary residues of each in the model's own numbering (drawn in the
published numbering, see numbering()). Every listed residue
is checked against the model before drawing; the ones that do not match, or
that the sheet flags, are drawn faded with a "?" and reported.

The structures are drawn by Mol*, in the style of the iGEM wiki's structure
viewer, with one representation per kind of thing so they cannot be mistaken
for each other (molstar/page.js has the details):

  protein             Mol*'s cartoon in the fold's colour; a second protein
                      chain (the SGNH dimer, the NylC beta subunit) lighter
  catalytic residues  pink sticks, labelled in bold
  secondary residues  thinner amber sticks (oxyanion hole, metal and heme
                      ligands, binding residues)
  modelled plastic    teal balls at 70% of van der Waals size (sticks for
                      residues, balls for the substrate)
  cofactors           heme as dark slate sticks, metal ions as silver balls

molstar/render.mjs runs it headless in Chrome with an orthographic camera.
Each structure is turned so the mechanism is in clear view (best_view): of
about 800 directions, the camera takes the one where the fewest atoms of the
catalytic residues, secondary residues, plastic and cofactors sit behind the
drawn backbone, with a slight preference for looking in from outside the
protein. The long axis of the fold then lies across the picture. Where the
site is buried and no direction is clear, both views cut away the residues
of cartoon that lie between the camera and the mechanism. A render
whose content touches the edge of the image is redone with a wider frame.

Writes runs/mechanisms/
  scenes.json           what the renderer draws, one entry per fold
  full/<key>.png        the whole structure, transparent
  zoom/<key>.png        the site alone, labelled
  views.json            where the zoom window and labels landed
  cards/<key>.{png,svg} full, zoom and the text of the mechanism together
  mechanism_sheet.{png,pdf}  all 15 cards
  sites.tsv             every residue drawn, with its check against the model

Usage: python 17_mechanism_cards.py [--no-render]   (reuse the last renders)
"""

import argparse
import importlib
import json
import os
import re
import subprocess
import textwrap

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["svg.hashsalt"] = "plasticome-mechanisms"
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

import folds
from common import INK, MUTED, SURFACE, write_tsv

r12 = importlib.import_module("12_render_structures")

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "runs", "mechanisms")
MOLSTAR = os.path.join(HERE, "molstar")
REFERENCE_PDB = os.path.join(HERE, "reference_pdb")      # the cited structures, from RCSB

# the wiki viewer's colours (StructureViewer.js), for the legend
CATALYTIC, SECONDARY, SUBSTRATE = "#e8559a", "#e0a526", "#1f9e89"
HEME, METAL = "#4a5560", "#b4b9bf"     # no brown: slate heme, silver metals
MECH_OFF = 0.2        # the mechanism may sit this share of the half width off centre
ZOOM_R = 11.0          # angstrom, half width of the zoom window
ZOOM_MAX = 17.0        # the window may grow to this to keep the site, plastic and cofactor in frame
METALS = {"FE": "Fe", "CU": "Cu", "ZN": "Zn", "MN": "Mn", "MG": "Mg", "CA": "Ca", "NI": "Ni"}
DPI = 450

# Cards are ordered as the tree is, largest fold first.
ORDER = ["3.40.50.1820", "3.60.70.12", "3.90.1300.10",
         "domain:oxygenase MpaB family protein", "3.40.50.200", "2.40.10.10",
         "3.40.710.10", "1.10.760.10", "2.60.40.420", "3.40.50.1110",
         "1.10.420.10;1.10.520.10", "1.10.630.10", "2.40.70.10", "3.20.20.370",
         "3.60.21.10"]

# PLASTICOME_GROUPING=mechanism: the cards of the mechanism groups, one fold per
# group in the group's colour, in runs/mechanisms-mechanism (ring tree, step 18)
GROUPED = os.environ.get("PLASTICOME_GROUPING")
if GROUPED:
    folds.use_grouping(GROUPED)
    OUTDIR = os.path.join(HERE, "runs", f"mechanisms-{GROUPED}")
    ORDER = sorted((k for k in ORDER if k in folds.GROUP_PICTURED), key=folds.type_rank)


def slug(k):
    return folds.name(k).lower().replace("/", "-").replace(" ", "_").replace("(", "").replace(")", "")


def rotation_to_z(v):
    """Rotation matrix taking unit vector v onto +z."""
    v = v / np.linalg.norm(v)
    z = np.array([0.0, 0.0, 1.0])
    c = float(v @ z)
    if c > 0.999999:
        return np.eye(3)
    if c < -0.999999:
        return np.diag([1.0, -1.0, -1.0])
    ax = np.cross(v, z)
    s = np.linalg.norm(ax)
    ax = ax / s
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K


def load(a):
    """Protein residues per chain, ligand atoms, cofactor residues."""
    model = next(iter(PDBParser(QUIET=True).get_structure("m", a["pdb"])))
    protein, ligand, cofactors = [], [], []
    for ch in model:
        for r in ch:
            het = r.id[0] != " "
            if not het and "CA" in r:
                protein.append(r)
            elif het and ch.id == "L":
                ligand.extend(r.get_atoms())
            elif het and r.id[0] != "W":
                cofactors.append(r)
    return model, protein, ligand, cofactors


def cofactor_name(r):
    """A hetero group by its atoms: the models call every one of them LIG."""
    atoms = list(r.get_atoms())
    el = [at.element.upper() for at in atoms]
    if len(atoms) == 1:
        return f"{METALS.get(el[0], el[0].capitalize())} ion", True
    if "FE" in el and len(atoms) >= 30:
        return "heme", False
    return f"cofactor ({len(atoms)} atoms)", False


def ligand_facts(a):
    """(name, formula, ipTM) of the modelled plastic, from the sheet."""
    lig, _ = folds.ligand_name(a)
    m = re.search(r"\(([A-Z][A-Za-z0-9]+)\)", a["ligand"])
    i = re.search(r"ligand-ipTM\s*([0-9.]+)", a["ligand"])
    return lig, m.group(1) if m else "", float(i.group(1)) if i else None


def numbering(k, model):
    """{(chain, model number): number as published}. The sheet lists residues
    in each model's own numbering, which for about half the folds starts at
    the mature chain or the first modelled residue (IsPETase S128 is S160 in
    the literature). The literature uses the numbering of the experimental
    structure, which is not always the full-length one either (subtilisin
    S221 and pepsin D32 count from the mature protein). So each model residue
    goes to its position in the full-length source protein, and from there to
    the author number the reference PDB entry gives that position; where the
    entry has no residue there, the entry's typical offset is used. For an
    AlphaFold reference the full-length (UniProt) position is the number."""
    src = folds.source_sequence(k)
    chains = [(ch.id, [r for r in ch if r.id[0] == " " and "CA" in r]) for ch in model]
    to_src = {}
    for cid, res in chains:
        if not res:
            continue
        aln = TAG_ALIGNER.align("".join(seq1(r.get_resname()) for r in res), src)[0]
        to_src.update({(cid, res[i].id[1]): j for i, j in zip(*aln.indices) if i >= 0 and j >= 0})
    path = os.path.join(REFERENCE_PDB, folds.structure_code(k) + ".cif")
    if not os.path.exists(path):
        return {key: str(j + 1) for key, j in to_src.items()}
    ref = next(iter(MMCIFParser(QUIET=True).get_structure("r", path)))
    best = {}
    for ch in ref:
        res = [r for r in ch if r.id[0] == " " and "CA" in r]
        if not res:
            continue
        aln = TAG_ALIGNER.align("".join(seq1(r.get_resname()) for r in res), src)[0]
        got = {int(j): f"{res[i].id[1]}{res[i].id[2].strip()}"
               for i, j in zip(*aln.indices) if i >= 0 and j >= 0}
        if len(got) > len(best):
            best = got
    off = int(np.median([int(re.match(r"-?\d+", v).group()) - j for j, v in best.items()]))
    return {key: best.get(j, str(j + off)) for key, j in to_src.items()}


def site_residues(model, a, published=None):
    """[(kind, group, residue or None, label, ok, got, flag)] for every
    annotated residue; labels in the published numbering when `published`
    (from numbering()) is given."""
    published = published or {}
    out = []
    for kind in ("primary", "secondary"):
        for group, res in a[kind]:
            for ch, aa, num, flag in res:
                try:
                    r = model[ch][num]
                    got = seq1(r.get_resname())
                except KeyError:
                    r, got = None, "?"
                n = published.get((ch, num), num)
                lab = f"{aa}{n}" if ch == "A" else f"{ch}:{aa}{n}"
                out.append((kind, group, r, lab, got == aa and not flag, got, flag))
    return out


TAG_ALIGNER = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1,
                                    end_gap_score=0)
TAG_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")


def not_tag(res, src):
    """The chain without terminal stretches that are not in the source
    protein: expression tags and linkers the model was built with (the
    IsPETase model ends in 23 residues of vector sequence)."""
    if not src:
        return res
    s = "".join(seq1(r.get_resname()) for r in res)
    aln = TAG_ALIGNER.align(s, src)[0]
    on = [int(i) for i, j in zip(*aln.indices) if i >= 0 and j >= 0]
    return res[min(on):max(on) + 1] if on else res


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], 1)


TUBE = 2.6        # angstrom: how close to a backbone atom a line of sight counts as blocked


def hidden(d, groups, occ):
    """Weighted share of the mechanism hidden behind the backbone, looking
    along -d (the camera is out along +d)."""
    total = 0.0
    for pts, w in groups:
        v = occ[None, :, :] - pts[:, None, :]            # site atom -> occluder
        t = v @ d                                         # distance towards the camera
        perp2 = (v ** 2).sum(-1) - t ** 2
        blocked = ((t > 1.5) & (perp2 < TUBE ** 2)).any(axis=1)
        total += w * blocked.mean()
    return total


def best_view(groups, occ, outward, offset=None, n=800):
    """The direction that shows the mechanism best, and the hidden share
    before (straight down the outward axis) and after. `offset` is the
    vector from the protein's centre to the mechanism's, over the protein's
    radius of gyration: directions that put the mechanism off to the side of
    the picture are penalised, so it sits in the middle of the protein."""
    outward = outward / np.linalg.norm(outward)
    D = fibonacci_sphere(n)
    score = np.array([hidden(d, groups, occ) for d in D])
    # a slight pull towards the outward view, so ties go to the familiar one
    total = score + 0.15 * (1 - D @ outward) / 2
    if offset is not None:
        side = np.sqrt(np.maximum((offset ** 2).sum() - (D @ offset) ** 2, 0))
        # weighted to win over a little occlusion: what still covers the
        # mechanism from a centred view is cut away (zoom_keep)
        total = total + 3.0 * side
    i = int(np.argmin(total))
    return D[i], hidden(outward, groups, occ), score[i]


def scene(k, a):
    """What the renderer needs for one fold, in the model's own coordinates."""
    model, protein, ligand, cofactors = load(a)
    published = numbering(k, model)
    sites = site_residues(model, a, published)
    ca = np.array([r["CA"].coord for r in protein], float)
    sc = [at.coord for kind, _, r, *_ in sites if r is not None and kind == "primary"
          for at in r if at.get_id() not in ("N", "C", "O")]
    centre = np.mean(sc, axis=0) if sc else ca.mean(axis=0)

    # what has to stay in sight, weighted: catalytic side chains most
    site_ids = {(r.get_parent().id, r.id[1]) for _, _, r, *_ in sites if r is not None}
    side = lambda kind: np.array([at.coord for kd, _, r, *_ in sites if r is not None and kd == kind
                                  for at in r if at.get_id() not in ("N", "C", "O", "CA")]
                                 or np.zeros((0, 3)), float)
    groups = [(side("primary"), 3.0), (side("secondary"), 1.0),
              (np.array([at.coord for at in ligand], float), 1.0),
              (np.array([at.coord for r in cofactors for at in r], float), 0.5)]
    groups = [(p, w) for p, w in groups if len(p)]
    # what can hide it: the backbone of every drawn residue that is not a site
    trimmed = {}
    for r in protein:
        trimmed.setdefault(r.get_parent().id, []).append(r)
    src = folds.source_sequence(k)
    drawn = [r for res in trimmed.values()
             for r in r12.trim_predicted(not_tag(res, src), "model")[0]]
    occ = np.array([r[n].coord for r in drawn if (r.get_parent().id, r.id[1]) not in site_ids
                    for n in ("N", "CA", "C", "O") if n in r], float)
    # the mechanism's centre: catalytic side chains and the plastic, equally
    lig = np.array([at.coord for at in ligand], float)
    mech_c = np.mean([centre, lig.mean(axis=0)], axis=0) if len(lig) else centre
    rg = np.sqrt(((ca - ca.mean(axis=0)) ** 2).sum(1).mean())
    face, before, after = best_view(groups, occ, centre - ca.mean(axis=0),
                                    offset=(mech_c - ca.mean(axis=0)) / rg)
    # the close up cuts away the stretches of cartoon still in front of the
    # mechanism from that direction
    occ_res = [r for r in drawn if (r.get_parent().id, r.id[1]) not in site_ids]
    mech = np.concatenate([p for p, _ in groups])
    near = mech[np.linalg.norm(mech - centre, axis=1) < ZOOM_R * 1.6]
    cut = set()
    for r in occ_res:
        bb = np.array([r[n].coord for n in ("N", "CA", "C", "O") if n in r], float)
        v = bb[None, :, :] - near[:, None, :]
        t = v @ face
        perp2 = (v ** 2).sum(-1) - t ** 2
        if ((t > 1.5) & (perp2 < TUBE ** 2)).any():
            cut.add((r.get_parent().id, r.id[1]))
    zoom_keep = []
    for r in drawn:
        key = (r.get_parent().id, r.id[1])
        if key in cut:
            continue
        last = zoom_keep[-1] if zoom_keep else None
        if last and last["chain"] == key[0] and last["end"] == key[1] - 1:
            last["end"] = key[1]
        else:
            zoom_keep.append({"chain": key[0], "beg": key[1], "end": key[1]})
    print(f"  {folds.name(k):<34} mechanism hidden {before:.2f} -> {after:.2f} "
          f"(weighted, of {sum(w for _, w in groups):.1f}); close up cuts away "
          f"{len(cut)} residues in front", flush=True)
    R = rotation_to_z(face)
    # in plane: long axis of the rotated CA cloud across the picture
    p = (ca - ca.mean(axis=0)) @ R.T
    _, _, vt = np.linalg.svd(p[:, :2] - p[:, :2].mean(axis=0), full_matrices=False)
    t = np.arctan2(vt[0, 1], vt[0, 0])
    Rz = np.array([[np.cos(-t), -np.sin(-t), 0], [np.sin(-t), np.cos(-t), 0], [0, 0, 1]])
    M = Rz @ R                      # rows: screen x, y, z as model directions

    # predicted models: drop dangling low confidence termini (step 12 rule),
    # so a signal peptide does not take over the picture
    chains = {}
    for r in protein:
        chains.setdefault(r.get_parent().id, []).append(r)
    keep, shown = [], []
    for i, (ch, res) in enumerate(chains.items()):
        n = len(res)
        res = r12.trim_predicted(not_tag(res, src), "model")[0]
        # the first chain in the fold's colour, any partner chain lighter
        keep.append({"chain": ch, "beg": int(res[0].id[1]), "end": int(res[-1].id[1]),
                     "color": folds.label(k) if i == 0 else folds.band(k), "length": n})
        shown += [r["CA"].coord for r in res]
    # the frame: near the mechanism, wide enough for everything drawn
    # (core, plastic, cofactors) as seen down the viewing axis
    shown += [at.coord for at in ligand] + [at.coord for r in cofactors for at in r]
    o = ca.mean(axis=0)
    q = (np.array(shown, float) - o) @ M.T
    # The frame centre may move from the mechanism toward the middle of the
    # drawn content, as far as keeps the mechanism within MECH_OFF of the
    # frame's half width from the centre: an active site off to one side
    # (the alpha/beta hydrolase) otherwise leaves the protein small in a
    # mostly empty frame.
    mid = (mech_c - o) @ M.T
    box = (q[:, :2].min(axis=0) + q[:, :2].max(axis=0)) / 2
    half, c = None, mid[:2]
    for t in np.linspace(0, 1, 41):
        ct = mid[:2] + t * (box - mid[:2])
        h = np.abs(q[:, :2] - ct).max() + 3.0
        if half is None or (h < half and np.linalg.norm(ct - mid[:2]) <= MECH_OFF * h):
            half, c = h, ct
    frame_centre = o + M[0] * c[0] + M[1] * c[1] + M[2] * mid[2]

    # The close up: ZOOM_R around the catalytic residues, widened and re-centred
    # on what it shows (site side chains, plastic, cofactors) when any of it
    # would fall outside, so nothing is cut by the frame edge.
    tips = [np.array(at.coord, float) for at in ligand]
    tips += [np.array(at.coord, float) for r in cofactors for at in r]
    for kind, group, r, lab, ok, got, flag in sites:
        if r is not None:
            tips += [np.array(at.coord, float) for at in r
                     if at.get_id() not in ("N", "C", "O", "CA")] or [np.array(r["CA"].coord, float)]
    qz = (np.array(tips, float) - centre) @ M.T
    zoom_c, zoom_r = centre, ZOOM_R
    if (np.abs(qz[:, :2]).max(axis=1) * 1.1 + 2.5 > ZOOM_R).any():
        zb = (qz[:, :2].min(axis=0) + qz[:, :2].max(axis=0)) / 2
        zoom_r = float(min(ZOOM_MAX, max(ZOOM_R, np.abs(qz[:, :2] - zb).max() * 1.1 + 2.5)))
        zoom_c = centre + M[0] * zb[0] + M[1] * zb[1]

    out_sites = []
    for kind, group, r, lab, ok, got, flag in sites:
        if r is None:
            continue
        side = [at for at in r if at.get_id() not in ("N", "C", "O", "CA")] or [r["CA"]]
        tip = max(side, key=lambda at: np.linalg.norm(at.coord - r["CA"].coord))
        out_sites.append({"chain": r.get_parent().id, "num": int(r.id[1]), "kind": kind,
                          "ok": bool(ok), "label": lab, "anchor": tip.coord.tolist()})
    cof = []
    for r in cofactors:
        name, metal = cofactor_name(r)
        cof.append({"chain": r.get_parent().id, "num": int(r.id[1]), "name": name,
                    "metal": metal,
                    "anchor": np.mean([at.coord for at in r], axis=0).tolist()})
    lig_name, formula, iptm = ligand_facts(a)
    return {
        "key": k, "slug": slug(k), "label": folds.name(k),
        "pdb_path": os.path.abspath(a["pdb"]),
        "color": folds.label(k),
        "keep": keep, "sites": out_sites, "cofactors": cof,
        "numbering": {f"{c}:{n}": v for (c, n), v in published.items()},
        "ligand": bool(ligand), "ligand_name": lig_name, "ligand_formula": formula,
        "ligand_short": re.split(r" / |, ", lig_name)[0],
        "ligand_iptm": iptm,
        "ligand_anchor": (np.mean([at.coord for at in ligand], axis=0).tolist()
                          if ligand else None),
        "zoom_keep": zoom_keep, "zoom_cut": len(cut),
        "frame_scale": 1.0, "hidden_before": round(before, 3), "hidden_after": round(after, 3),
        "centre": centre.tolist(), "face": M[2].tolist(), "up": M[1].tolist(),
        "frame_centre": frame_centre.tolist(), "frame_radius": float(half),
        "scene_radius": float(np.linalg.norm(np.array(shown, float) - frame_centre, axis=1).max() + 5),
        "zoom_radius": zoom_r, "zoom_centre": zoom_c.tolist(),
    }, sites


def tidy(group):
    return (group.replace("Stabillizing", "Stabilising").replace("Nucelophile", "Nucleophile")
            .replace("Proximinal", "Proximal"))


def contents(sc, a):
    """What is in one structure, as (role, how it is drawn, what it is) rows.
    The cards print these, and the ring tree editor shows them on hover."""
    code = folds.structure_code(sc["key"])
    rows = [("representative structure", "", "co-folded model of "
             + folds.representative(sc["key"]) + "; residues numbered as in "
             + (f"PDB {code}" if not code.startswith("AF-") else "UniProt (full length)"))]
    prot = [f"chain {c['chain']}, {c['length']} aa" for c in sc["keep"]]
    rows.append(("protein", "cartoon",
                 prot[0] + (f"; {', '.join(prot[1:])} (second chain, drawn lighter)"
                            if len(prot) > 1 else "")))
    pub = lambda ch, n: sc["numbering"].get(f"{ch}:{n}", n)
    res = lambda kind: "; ".join(
        tidy(g) + " " + " ".join((f"{aa}{pub(ch, n)}" if ch == "A" else f"{ch}:{aa}{pub(ch, n)}")
                                 for ch, aa, n, _ in r) for g, r in a[kind])
    rows.append(("catalytic residues", "pink sticks", res("primary")))
    if a["secondary"]:
        rows.append(("secondary residues", "thinner amber sticks", res("secondary")))
    if sc["ligand"]:
        _, plastic = folds.ligand_name(a)
        rows.append(("modelled plastic", "teal balls",
                     f"{sc['ligand_name']} ({sc['ligand_formula']}), standing in for {plastic}, "
                     f"co-folded with the protein as chain L" + (f"; ligand ipTM {sc['ligand_iptm']:.2f}"
                                                 if sc["ligand_iptm"] else "")))
    if sc["cofactors"]:
        n = {}
        for c in sc["cofactors"]:
            n[c["name"]] = n.get(c["name"], 0) + 1
        rows.append(("cofactors", "heme as slate sticks, metal ions as silver balls",
                     ", ".join(f"{v} x {k}" if v > 1 else k for k, v in n.items())))
    return rows


SWATCH = {"protein": None, "catalytic residues": CATALYTIC, "secondary residues": SECONDARY,
          "modelled plastic": SUBSTRATE}


def card(k, a, sc, view, path_stem):
    """Full view, zoom, and the key to what is in them, for one fold. Laid out
    in inches; the card is as tall as its key needs."""
    W, HEAD, IMG, CAP = 6.4, 0.62, 3.0, 0.24
    rows = []
    for role, how, what in contents(sc, a):
        body = textwrap.fill(what, 118)
        rows.append((role, how, body, 0.19 + 0.125 * (body.count("\n") + 1) + 0.07))
    flagged = any(not st["ok"] for st in sc["sites"])
    H = HEAD + IMG + CAP + sum(r[3] for r in rows) + (0.2 if flagged else 0) + 0.1
    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor(SURFACE)
    fx, fy = (lambda x: x / W), (lambda y: 1 - y / H)      # inches from the top left
    s = slug(k)
    for n, sub in enumerate(("full", "zoom")):
        x0 = 0.1 + n * (W / 2)
        ax = fig.add_axes([fx(x0), fy(HEAD + IMG), IMG / W, IMG / H])
        size = view[sub]["size"]
        ax.imshow(plt.imread(os.path.join(OUTDIR, sub, f"{s}.png")),
                  extent=[0, size, size, 0], interpolation="lanczos")
        if sub == "full":
            ring = view["full"]["ring"]
            ax.add_patch(plt.Circle((ring["x"], ring["y"]), ring["r"], fill=False, ec=MUTED,
                                    lw=0.8, ls=(0, (2, 1.5))))
        else:
            ax.add_patch(plt.Rectangle((8, 8), size - 16, size - 16, fill=False, ec=MUTED,
                                       lw=0.8, ls=(0, (2, 1.5))))
        ax.set_xlim(0, size)
        ax.set_ylim(size, 0)
        ax.set_aspect("equal")
        ax.axis("off")
        fig.text(fx(x0 + IMG / 2), fy(HEAD + IMG + 0.05),
                 "whole model (dashed circle: the close up)" if sub == "full"
                 else "active site, close up", fontsize=5.8, color=MUTED, ha="center", va="top")
    fig.text(fx(0.12), fy(0.1), folds.name(k), fontsize=11, color=folds.strong(k),
             fontweight="bold", va="top")
    fig.text(fx(0.12), fy(0.36), f"{folds.representative(k, short=True)}   |   "
             f"{folds.mechanism(k)}   |   CATH "
             f"{k.split('|')[0] if not k.startswith('domain') else 'none'}"
             f"   |   reference {a['reference']}", fontsize=6.4, color=MUTED, va="top")
    # the key: one row per kind of thing in the picture
    y = HEAD + IMG + CAP
    for role, how, body, h in rows:
        col = SWATCH.get(role)
        if role == "protein":
            fig.add_artist(plt.Line2D([fx(0.14), fx(0.32)], [fy(y + 0.07)] * 2, lw=3.2,
                                      color=folds.label(k), solid_capstyle="round"))
        elif role == "cofactors":
            fig.add_artist(plt.Line2D([fx(0.12), fx(0.22)], [fy(y + 0.07)] * 2, lw=2.4,
                                      color=HEME, solid_capstyle="round"))
            fig.add_artist(plt.Circle((fx(0.29), fy(y + 0.07)), 0.05 / W, color=METAL,
                                      ec="#7d838a", lw=0.4, transform=fig.transFigure))
        elif col:
            fig.add_artist(plt.Circle((fx(0.23), fy(y + 0.07)), 0.055 / W, color=col,
                                      transform=fig.transFigure))
        fig.text(fx(0.42), fy(y), role.capitalize() + (f" ({how})" if how else ""),
                 fontsize=6.6, fontweight="bold", color=INK, va="top")
        fig.text(fx(0.42), fy(y + 0.19), body, fontsize=6.3, color=INK, va="top",
                 linespacing=1.35)
        y += h
    if flagged:
        fig.text(fx(0.42), fy(y + 0.02), "Faded, with \"?\": flagged in the annotation sheet, "
                 "needs a check.", fontsize=6.0, color=MUTED, va="top", style="italic")
    for ext in ("png", "svg"):
        fig.savefig(f"{path_stem}.{ext}", dpi=DPI, facecolor=SURFACE,
                    **({"metadata": {"Date": None}} if ext == "svg" else {}))
    plt.close(fig)


def margin(path):
    """Pixels between the drawn content and the nearest edge of a render."""
    from PIL import Image
    alpha = np.array(Image.open(path))[:, :, 3] > 10
    ys, xs = np.where(alpha)
    if not len(xs):
        return 0
    return int(min(xs.min(), ys.min(), alpha.shape[1] - 1 - xs.max(), alpha.shape[0] - 1 - ys.max()))


def render(scenes, path):
    """Run Mol* on scenes, then widen the frame of any full view whose content
    touches the edge and run those again, up to three times."""
    if not os.path.isdir(os.path.join(MOLSTAR, "node_modules")):
        subprocess.run(["npm", "install", "--no-audit", "--no-fund"], cwd=MOLSTAR, check=True)
    todo = scenes
    for attempt in range(4):
        with open(path, "w") as fh:
            json.dump(todo, fh, indent=1)
        subprocess.run(["node", os.path.join(MOLSTAR, "render.mjs"), path, OUTDIR], check=True)
        todo = [sc for sc in todo if margin(os.path.join(OUTDIR, "full", f"{sc['slug']}.png")) < 24]
        if not todo or attempt == 3:
            break
        for sc in todo:
            sc["frame_scale"] *= 1.15
        print(f"  widening the frame of {', '.join(sc['slug'] for sc in todo)}")
    with open(path, "w") as fh:
        json.dump(scenes, fh, indent=1)
    for sc in todo:
        print(f"  still touching the edge: {sc['slug']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-render", action="store_true", help="reuse the last Mol* renders")
    args = ap.parse_args()
    os.makedirs(os.path.join(OUTDIR, "cards"), exist_ok=True)
    A = folds.annotations()
    table, scenes, all_sites = [], [], {}
    for k in ORDER:
        sc, sites = scene(k, A[k])
        sc["contents"] = contents(sc, A[k])
        scenes.append(sc)
        all_sites[k] = sites
    path = os.path.join(OUTDIR, "scenes.json")
    if args.no_render:
        with open(path, "w") as fh:
            json.dump(scenes, fh, indent=1)
    else:
        render(scenes, path)
    with open(os.path.join(OUTDIR, "views.json")) as fh:
        views = json.load(fh)

    for k, sc in zip(ORDER, scenes):
        s = slug(k)
        card(k, A[k], sc, views[s], os.path.join(OUTDIR, "cards", s))
        sites = all_sites[k]
        for kind, group, r, lab, ok, got, flag in sites:
            table.append([k, folds.name(k), kind, group, lab,
                          f"{r.get_parent().id}:{r.id[1]}" if r is not None else "", got,
                          "match" if got == lab.split(":")[-1][0] else "MISMATCH", flag])
        bad = [lab for kind, group, r, lab, ok, got, flag in sites if not ok]
        print(f"  {folds.name(k):<34} {len(sites):>2} residues, margin "
              f"{margin(os.path.join(OUTDIR, 'full', s + '.png'))} px"
              + (f"  not verified: {', '.join(bad)}" if bad else ""), flush=True)
    write_tsv(os.path.join(OUTDIR, "sites.tsv"),
              ["fold_key", "fold", "kind", "group", "residue", "model_number", "residue_in_model",
               "check", "flag_in_sheet"], table)

    # all cards on one sheet, three columns, each column stacked
    from PIL import Image
    ims = [Image.open(os.path.join(OUTDIR, "cards", f"{slug(k)}.png")).convert("RGB")
           for k in ORDER]
    cols = 3
    w = ims[0].width
    heights = [sum(im.height for im in ims[c::cols]) for c in range(cols)]
    sheet = Image.new("RGB", (w * cols, max(heights)), SURFACE)
    for c in range(cols):
        y = 0
        for im in ims[c::cols]:
            sheet.paste(im, (c * w, y))
            y += im.height
    sheet = sheet.resize((sheet.width // 2, sheet.height // 2), Image.LANCZOS)
    sheet.save(os.path.join(OUTDIR, "mechanism_sheet.png"))
    sheet.save(os.path.join(OUTDIR, "mechanism_sheet.pdf"), resolution=150)
    print(f"{len(ORDER)} cards -> {os.path.relpath(OUTDIR, HERE)}")


if __name__ == "__main__":
    main()
