#!/usr/bin/env python3
"""
Step 17b: the six mechanism-group structures redrawn with three colour roles,
ranked mechanism > structure > substrate (adopted for the vertical tree on
2026-11-14):

  mechanism   catalytic and secondary residues, heme and metal ions, all in
              the group's lighter shade (folds.label; the strong shade of the
              badges sank into the grey for groups 2 and 5). Catalytic
              sticks thicker, secondary thinner.
  structure   the protein cartoon in one light grey for every group
  substrate   the modelled plastic as thin off-white sticks

Mol*'s dark outline is on so the off-white plastic reads against the grey.
Reuses step 17's scenes for the mechanism grouping (same sites), so run that
first: PLASTICOME_GROUPING=mechanism python 17_mechanism_cards.py

--spread (2026-11-14) also picks a new camera for each structure, so the
catalytic residues are easy to tell apart:
  centred      the frame is centred on the catalytic side chains, so the
               mechanism sits in the middle of the picture
  apart        of about 1500 directions, the camera takes the one where the
               catalytic residues overlap least on screen (every pair of side
               chains at least SEP apart in projection, scored by how far
               short they fall), with a small reward for spreading them out
               and a small penalty for the share hidden behind the backbone. With one catalytic residue (NylC), the
               secondary residues and cofactors are spread around it instead.
  tilt         a structure in TILT is turned that many degrees further, towards
               the face of its cofactor (panel 6: the heme was seen edge on)
  see-through  backbone still in front of the mechanism from that direction
               is drawn at SEE_ALPHA opacity instead of removed (step 17
               removes it), and so are the plastic's atoms in front of it

Writes runs/mechanisms-mechanism-grey/ (step 17's views) or
runs/mechanisms-mechanism-grey-spread/ (--spread):
  scenes.json, views.json   the scenes as rendered
  full/<slug>.png           whole structure, transparent (used by step 18 --vertical)
  zoom/<slug>.png           the active site close up (unlabelled with --spread, for step 18 --insets)

Usage: python 17b_mechanism_colour_renders.py [--spread] [--no-render]
"""

import argparse
import importlib
import itertools
import json
import os
import subprocess

import numpy as np
from Bio.PDB import PDBParser

import folds

s17 = importlib.import_module("17_mechanism_cards")

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.join(HERE, "runs", "mechanisms-mechanism")
OUTDIR = os.path.join(HERE, "runs", "mechanisms-mechanism-grey")
CARTOON = "#e3e5e9"      # protein: light grey (#c8ccd2 rendered too dark under Mol* shading)
SUBSTRATE = "#f2efe8"    # plastic: off-white, outlined
OUTLINE = "#2b2f36"      # Mol*'s outline colour (molstar/page.js)

N_DIRS = 1500
SEP = 3.5                # angstrom: side chains this far apart on screen read as separate sticks
W_HIDDEN = 0.3           # weight of the hidden share (what still hides it is cut away anyway)
W_SPREAD = 0.1           # reward for side chains spread out on screen, per 10 angstrom
PAD = 3.0                # angstrom: margin around the drawn content
TILT = {"oxygenase_mpab_lcp": 20}   # degrees towards the heme's face; the best spread view sees it edge on
UPRIGHT = False          # long axis vertical; tried 2026-11-14 and rejected, keep the views as they were
SUBSTRATE_ALPHA = 0.5    # the plastic's atoms in front of the mechanism (0.25 left it almost invisible)
SEE_ALPHA = 0.25         # opacity of what sits in front of the mechanism (backbone, plastic)
ZOOM_PAD, ZOOM_MIN = 2.0, 6.0   # angstrom: close up margin past the farthest mechanism atom, smallest radius


def colour(s):
    mech = folds.label(s["key"])     # the group's lighter shade: the dark ones (2, 5) sank into the grey
    s["cartoon_override"] = CARTOON
    s["palette"] = {"catalytic": mech, "secondary": mech, "heme": mech,
                    "metal": mech, "substrate": SUBSTRATE}
    s["catalytic_size"], s["secondary_size"] = 1.3, 0.7
    s["substrate_style"] = "sticks"
    s["outline"] = True


def basis(d):
    """Two unit vectors across the screen when looking along -d."""
    a = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(d, a)
    u /= np.linalg.norm(u)
    return u, np.cross(d, u)


def overlap(d, items, spread=False):
    """How far short of SEP the closest atoms of each pair of items fall on
    screen, squared and weighted, summed over pairs."""
    u, v = basis(d)
    flat = [(np.stack([p @ u, p @ v], 1), w) for p, w in items]
    out = 0.0
    for (a, wa), (b, wb) in itertools.combinations(flat, 2):
        if spread:
            out -= W_SPREAD * wa * wb * np.linalg.norm(a.mean(0) - b.mean(0)) / 10
            continue
        dmin = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1).min())
        out += wa * wb * max(0.0, SEP - dmin) ** 2 / SEP ** 2
    return out


def spread_view(s):
    """A new camera for scene s: centred on the catalytic side chains, from the
    direction where they overlap least and are least hidden."""
    st = PDBParser(QUIET=True).get_structure("m", s["pdb_path"])
    model = next(iter(st))
    res = {(r.get_parent().id, r.id[1]): r for r in model.get_residues()}
    side = lambda r: np.array([a.coord for a in r if a.get_id() not in ("N", "C", "O")], float)
    keys = lambda kind: list(dict.fromkeys((x["chain"], x["num"]) for x in s["sites"]
                                           if x["kind"] == kind and x["ok"]))
    prim, sec = keys("primary"), [k for k in keys("secondary") if k not in keys("primary")]
    cof = [(c["chain"], c["num"]) for c in s["cofactors"]]
    items = [(side(res[k]), 1.0) for k in prim]
    if len(items) < 2:                               # one catalytic residue: spread its partners
        items += [(side(res[k]), 0.5) for k in sec]
    items += [(np.array([a.coord for a in res[k]], float), 0.5) for k in cof]
    centre = np.concatenate([p for p, w in items if w == 1.0]).mean(axis=0)

    # drawn residues (step 17's trim) and what can hide the mechanism
    drawn = [res[(k["chain"], n)] for k in s["keep"] for n in range(k["beg"], k["end"] + 1)
             if (k["chain"], n) in res]
    site_ids = {(x["chain"], x["num"]) for x in s["sites"]}
    occ_res = [r for r in drawn if (r.get_parent().id, r.id[1]) not in site_ids]
    occ = np.array([r[n].coord for r in occ_res for n in ("N", "CA", "C", "O") if n in r], float)
    groups = [(np.concatenate([p for p, w in items if w == 1.0]), 3.0)]
    if sec:
        groups.append((np.concatenate([side(res[k]) for k in sec]), 1.0))

    D = s17.fibonacci_sphere(N_DIRS)
    score = np.array([overlap(d, items) + overlap(d, items, True) + W_HIDDEN * s17.hidden(d, groups, occ)
                      for d in D])
    face = D[int(np.argmin(score))]
    if s["slug"] in TILT and cof:
        # turn the camera TILT degrees towards the face of the cofactor's plane
        # (the heme seen edge on shows only its rim)
        het = np.concatenate([np.array([a.coord for a in res[k]], float) for k in cof])
        n = np.linalg.svd(het - het.mean(0), full_matrices=False)[2][2]
        n = n if n @ face >= 0 else -n
        side_ = n - (n @ face) * face
        if np.linalg.norm(side_) > 1e-6:
            side_ /= np.linalg.norm(side_)
            before_ = np.degrees(np.arccos(np.clip(abs(n @ face), 0, 1)))
            t_ = np.radians(TILT[s["slug"]])
            face = np.cos(t_) * face + np.sin(t_) * side_
            print(f"  {s['slug']:<22} tilted {TILT[s['slug']]} deg: cofactor plane "
                  f"{90 - before_:.0f} -> {90 - before_ + TILT[s['slug']]:.0f} deg to the screen")
    old = np.array(s["face"])

    # cut away backbone in front of the mechanism, as step 17 does
    mech = np.concatenate([p for p, _ in groups])
    cut = set()
    for r in occ_res:
        bb = np.array([r[n].coord for n in ("N", "CA", "C", "O") if n in r], float)
        v = bb[None, :, :] - mech[:, None, :]
        t = v @ face
        if ((t > 1.5) & ((v ** 2).sum(-1) - t ** 2 < s17.TUBE ** 2)).any():
            cut.add((r.get_parent().id, r.id[1]))
    # --spread: the backbone in front is drawn see-through rather than removed
    # (one residue of overlap at each end so the cartoon stays joined), and so
    # are the plastic's atoms in front of the mechanism
    see = []
    for r in drawn:
        key = (r.get_parent().id, r.id[1])
        if not any((key[0], key[1] + o) in cut for o in (-1, 0, 1)):
            continue
        if see and see[-1]["chain"] == key[0] and see[-1]["end"] == key[1] - 1:
            see[-1]["end"] = key[1]
        else:
            see.append({"chain": key[0], "beg": key[1], "end": key[1]})
    front = []
    if s["ligand"]:
        for a in model["L"].get_atoms():
            v = a.coord - mech
            t = v @ face
            if ((t > 1.5) & ((v ** 2).sum(-1) - t ** 2 < s17.TUBE ** 2)).any():
                front.append(a.get_id())
    keep = []
    for r in drawn:
        key = (r.get_parent().id, r.id[1])
        if key in cut:
            continue
        if keep and keep[-1]["chain"] == key[0] and keep[-1]["end"] == key[1] - 1:
            keep[-1]["end"] = key[1]
        else:
            keep.append({"chain": key[0], "beg": key[1], "end": key[1]})

    # the protein's long axis up (UPRIGHT) or across the picture, the frame centred on the mechanism
    u, v = basis(face)
    ca = np.array([r["CA"].coord for r in drawn if "CA" in r], float)
    q = np.stack([(ca - ca.mean(0)) @ u, (ca - ca.mean(0)) @ v], 1)
    long_ = np.linalg.svd(q, full_matrices=False)[2][0]
    if UPRIGHT:          # long axis up the picture: the grid columns are narrow, the rows tall
        up = long_[0] * u + long_[1] * v
        right = np.cross(up, face)
    else:
        right = long_[0] * u + long_[1] * v
        up = np.cross(face, right)
    shown = [ca]
    if s["ligand"]:
        shown.append(np.array([a.coord for a in model["L"].get_atoms()], float))
    shown += [p for p, w in items if w < 1.0]
    pts = np.concatenate(shown) - centre
    half = max(np.abs(pts @ right).max(), np.abs(pts @ up).max()) + PAD

    # the close up: every mechanism atom (catalytic and secondary side chains,
    # heme, metal ions) inside the round inset of step 18. Centred on the
    # middle of their extent on screen, radius to the farthest atom plus a margin
    allmech = [side(res[k]) for k in prim + sec] + [np.array([a.coord for a in res[k]], float)
                                                   for k in cof]
    flat = np.concatenate(allmech) - centre
    xy = np.stack([flat @ right, flat @ up], 1)
    mid = (xy.min(0) + xy.max(0)) / 2
    zr = float(max(ZOOM_MIN, np.linalg.norm(xy - mid, axis=1).max() + ZOOM_PAD))
    zoom_centre = centre + right * mid[0] + up * mid[1]
    s.update({"zoom_radius": zr})
    s.update({"face": face.tolist(), "up": up.tolist(), "centre": centre.tolist(),
              "zoom_centre": zoom_centre.tolist(), "frame_centre": centre.tolist(),
              "frame_radius": float(half), "frame_scale": 1.0, "zoom_keep": keep,
              "scene_radius": float(np.linalg.norm(pts, axis=1).max() + 5),
              "zoom_labels": False, "full_label_positions": True, "see_through": see, "see_through_alpha": SEE_ALPHA, "substrate_alpha": SUBSTRATE_ALPHA,
              "substrate_front": front,
              "substrate_back": [a.get_id() for a in model["L"].get_atoms() if a.get_id() not in front]
              if s["ligand"] else []})          # the close up is an inset in the tree, too small for labels
    print(f"  {s['slug']:<22} overlap {overlap(old, items):.2f} -> {overlap(face, items):.2f}, "
          f"hidden {s17.hidden(old, groups, occ):.2f} -> {s17.hidden(face, groups, occ):.2f}, "
          f"turned {np.degrees(np.arccos(np.clip(face @ old, -1, 1))):.0f} deg, "
          f"see-through {len(cut)} residues and {len(front)} plastic atoms, close up {zr:.1f} A",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spread", action="store_true")
    ap.add_argument("--no-render", action="store_true")
    args = ap.parse_args()
    out = OUTDIR + ("-spread" if args.spread else "")
    folds.use_grouping("mechanism")
    with open(os.path.join(SOURCE, "scenes.json")) as fh:
        scenes = json.load(fh)
    for s in scenes:
        colour(s)
        if args.spread:
            spread_view(s)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "scenes.json")
    with open(path, "w") as fh:
        json.dump(scenes, fh, indent=1)
    if not args.no_render:
        subprocess.run(["node", os.path.join(HERE, "molstar", "render.mjs"), path, out], check=True)


if __name__ == "__main__":
    main()
