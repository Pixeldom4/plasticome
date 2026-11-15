#!/usr/bin/env python3
"""
Step 12: the fold representatives drawn as structures, one SVG and PNG each.

Every structure is tinted with its own fold colour, the same hex the slice
carries in the pie (08_draw_pie.FOLD_TINT, the `fill` column of
cath_folds.tsv), so a panel of structures and the pie read as one figure.

Where the structures come from, in the order tried:

  rcsb    the representative's PDB entry, when the union table gives one
  afdb    an AlphaFold DB model whose sequence matches exactly
  esmfold the ESM Atlas API, for a sequence under 400 residues with neither

Of the ten, eight have a PDB entry, LMC021 takes an AlphaFold DB model and
PLAase I is folded by ESMFold. The source of every image is recorded in
structures/sources.tsv, so a predicted model is never mistaken for a
deposition.

Drawing. There is no PyMOL or ChimeraX on this machine, so the cartoon is
drawn here: the CA trace, smoothed to a spline, as a tube whose width follows
the secondary structure and whose lightness follows depth. Segments are drawn
back to front, which is what gives the crossings. Secondary structure comes
from the HELIX and SHEET records when the file has them, otherwise from the
CA(i) to CA(i+3) and CA(i+4) distances, which separates helix from strand
well enough for a cartoon.

The view is the one that shows the most: coordinates are centred and rotated
onto their principal axes, then projected down the axis of least variance.

Writes into the bundle, structures/{pdb,svg,png}/ plus sources.tsv and the
contact sheet folds_contact_sheet.{svg,png}.
"""

import argparse
import collections
import importlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from scipy.interpolate import splev, splprep

from common import OUT, load_rows, read_tsv, write_tsv

st = importlib.import_module("03_structures")
pie = importlib.import_module("08_draw_pie")

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.path.join(HERE, "..", "structures", "fold-representatives")
SDIR = os.path.join(BUNDLE, "structures")
RCSB = "https://files.rcsb.org/download/{acc}.pdb"
ESM_MAX = 400

# tube width in points per secondary structure, and the spline smoothing
WIDTH = {"H": 10.0, "E": 7.5, "L": 3.6}
DEPTH = (0.74, 1.10)        # lightness factor back to front
SAMPLES = 8                 # spline points per residue
BANDS = 26                  # depth bands, see draw()
DPI = 400

ALIGNER = Align.PairwiseAligner(mode="global", open_gap_score=-11,
                                extend_gap_score=-1, end_gap_score=0)
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")


# --- getting a structure -----------------------------------------------------

def structure_for(rep, row, seq):
    """(pdb text, source, note) for one representative."""
    acc = rep["rep_pdb"].strip()
    if rep["rep_has_structure"].strip().lower() == "yes" and acc not in ("", "NA", "-"):
        code = acc.split("_")[0][:4]
        txt = st.fetch(RCSB.format(acc=code), f"rcsb_{code}.pdb")
        if txt:
            return txt, "rcsb", code
    for a in st.accessions(row, seq):
        js = st.fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{a}", f"afdb_{a}.json")
        if not js:
            continue
        for m in json.loads(js):
            if st.match(seq, m.get("sequence", "")) == "exact" and m.get("pdbUrl"):
                return st.fetch(m["pdbUrl"], f"afdb_{a}.pdb"), "afdb", a
    if len(seq) <= ESM_MAX:
        pdb, plddt = st.esmfold(seq, f"esmfold_{rep['rep_plasticome_id']}.pdb")
        return pdb, "esmfold", f"mean pLDDT {plddt:.1f}"
    raise SystemExit(f"{rep['slice']}: no structure and {len(seq)} aa is over the "
                     f"ESMFold limit of {ESM_MAX}")


def trim_predicted(res, source, cut=70.0, keep=0.5):
    """Trim dangling termini from a predicted model.

    A predicted model often carries a signal peptide or a propeptide as an
    extended arm off the folded core. It is real output, but it takes over the
    picture and it is not the domain the fold is about, so the leading and
    trailing dangling stretches come off. Crystal structures are left alone,
    since the disordered parts are already missing from the deposition.

    Two rules, applied from each end. A residue goes if its pLDDT, in the B
    factor column, is under `cut`, or if it is dangling: fewer than two CA
    atoms within 10 A once its own sequence neighbours are discounted, which
    is what an extended arm looks like and a folded terminus does not. ESMFold
    can be confident about such an arm, so the contact rule is the one that
    catches it. Nothing inside the chain is touched, and if trimming would
    leave under half the chain it is abandoned.
    """
    if source == "rcsb":
        return res, 0, 0
    ca = np.array([r["CA"].coord for r in res], dtype=float)
    conf = [float(np.mean([a.get_bfactor() for a in r])) for r in res]
    conf = [c * 100 if c <= 1.0 else c for c in conf]
    d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=-1)
    near = (d < 10.0)
    for k in range(-4, 5):                       # discount sequence neighbours
        np.fill_diagonal(near[max(0, -k):, max(0, k):], False)
    contacts = near.sum(axis=1)

    def dangling(i):
        return conf[i] < cut or contacts[i] < 2

    i, j = 0, len(res)
    while i < j and dangling(i):
        i += 1
    while j > i and dangling(j - 1):
        j -= 1
    if j - i < keep * len(res):
        return res, 0, 0
    return res[i:j], i, len(res) - j


def best_chain(pdb_text, seq, tag):
    """The chain whose sequence matches the representative best."""
    path = os.path.join(SDIR, "pdb", f"{tag}.pdb")
    with open(path, "w") as fh:
        fh.write(pdb_text)
    model = next(iter(PDBParser(QUIET=True).get_structure(tag, path)))
    best = None
    for ch in model:
        res = [r for r in ch if is_aa(r, standard=True) and "CA" in r]
        if len(res) < 30:
            continue
        s = "".join(seq1(r.get_resname()) for r in res)
        aln = ALIGNER.align(seq, s)[0]
        same = sum(1 for a, b in zip(aln[0], aln[1]) if a == b and a != "-")
        ident = same / max(len(seq), len(s))
        if best is None or ident > best[0]:
            best = (ident, ch.id, res)
    if best is None:
        raise SystemExit(f"{tag}: no protein chain of 30 or more residues")
    return best[1], best[2], best[0], path


# --- secondary structure -----------------------------------------------------

def ss_from_records(path, chain, res):
    """HELIX and SHEET records, when the file carries them."""
    want = {"H": set(), "E": set()}
    for line in open(path):
        if line.startswith("HELIX") and line[19] == chain:
            want["H"] |= set(range(int(line[21:25]), int(line[33:37]) + 1))
        elif line.startswith("SHEET") and line[21] == chain:
            want["E"] |= set(range(int(line[22:26]), int(line[33:37]) + 1))
    if not want["H"] and not want["E"]:
        return None
    out = []
    for r in res:
        n = r.id[1]
        out.append("H" if n in want["H"] else "E" if n in want["E"] else "L")
    return "".join(out)


def ss_from_ca(ca):
    """Helix and strand from CA(i)->CA(i+3) and CA(i+4) distances."""
    n = len(ca)
    ss = ["L"] * n
    d3 = [np.linalg.norm(ca[i + 3] - ca[i]) if i + 3 < n else 99 for i in range(n)]
    d4 = [np.linalg.norm(ca[i + 4] - ca[i]) if i + 4 < n else 99 for i in range(n)]
    for i in range(n):
        if 4.7 <= d3[i] <= 6.3 and 5.0 <= d4[i] <= 7.6:
            for j in range(i, min(i + 4, n)):
                ss[j] = "H"
    for i in range(n):
        if ss[i] == "L" and 9.2 <= d3[i] <= 11.2:
            for j in range(i, min(i + 3, n)):
                if ss[j] == "L":
                    ss[j] = "E"
    # drop runs of one or two, they are noise at this scale
    out, i = list(ss), 0
    while i < n:
        j = i
        while j < n and out[j] == out[i]:
            j += 1
        if out[i] != "L" and j - i < 3:
            for k in range(i, j):
                out[k] = "L"
        i = j
    return "".join(out)


# --- drawing -----------------------------------------------------------------

def project(ca):
    """Centre, rotate onto principal axes, and return x, y, depth."""
    x = ca - ca.mean(axis=0)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    p = x @ vt.T                      # columns ordered by variance
    return p[:, 0], p[:, 1], p[:, 2]


def draw(ax, ca, ss, tint, lw_scale=1.0):
    """The CA trace as a depth shaded tube, back to front.

    Drawing one path per spline segment gives a megabyte of SVG and tens of
    thousands of elements, which no illustration program enjoys. Depth is
    banded into `BANDS` steps instead, and each run of consecutive segments
    that shares a band and a secondary structure becomes one polyline. The
    picture is the same to the eye and the file is two orders of magnitude
    smaller.
    """
    t = np.arange(len(ca))
    k = min(3, len(ca) - 1)
    tck, _ = splprep([ca[:, 0], ca[:, 1], ca[:, 2]], u=t, s=0, k=k)
    fine = np.linspace(0, len(ca) - 1, max(2, int(len(ca) * SAMPLES)))
    sx, sy, sz = splev(fine, tck)
    x, y, z = project(np.column_stack([sx, sy, sz]))

    zmid = (z[:-1] + z[1:]) / 2
    lo, span = zmid.min(), (zmid.max() - zmid.min()) or 1.0
    band = np.clip(((zmid - lo) / span * BANDS).astype(int), 0, BANDS - 1)
    kind = np.array([ss[min(int(round(f)), len(ss) - 1)] for f in fine[:-1]])

    # maximal runs of consecutive segments sharing a band and a structure
    runs = []
    i = 0
    while i < len(band):
        j = i + 1
        while j < len(band) and band[j] == band[i] and kind[j] == kind[i]:
            j += 1
        runs.append((band[i], kind[i], i, j))
        i = j

    for b in range(BANDS):
        here = [r for r in runs if r[0] == b]
        if not here:
            continue
        f = (b + 0.5) / BANDS
        light = DEPTH[0] + f * (DEPTH[1] - DEPTH[0])
        for pass_, (color, extra, zo) in enumerate((
                (pie.shade(tint, light * 0.52), 2.0, 2 * b + 1),
                (pie.shade(tint, light), 0.0, 2 * b + 2))):
            for _, sstype, i, j in here:
                w = WIDTH[sstype] * lw_scale + extra
                ax.plot(x[i:j + 1], y[i:j + 1], color=color, lw=w, zorder=zo,
                        solid_capstyle="round", solid_joinstyle="round",
                        snap=False)

    # a square window around the centre, so every panel is the same box and
    # the structures stay to one scale against each other
    half = 0.53 * max(np.ptp(x), np.ptp(y))
    cx, cy = (x.min() + x.max()) / 2, (y.min() + y.max()) / 2
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
    ax.axis("off")
    return np.ptp(x), np.ptp(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-esmfold", action="store_true")
    ap.add_argument("--size", type=float, default=3.2, help="inches per panel")
    args = ap.parse_args()
    if os.environ.get("PLASTICOME_TREE_SET") != "centroids":
        raise SystemExit("run with PLASTICOME_TREE_SET=centroids")
    for sub in ("pdb", "svg", "png"):
        os.makedirs(os.path.join(SDIR, sub), exist_ok=True)
    os.makedirs(st.CACHE, exist_ok=True)

    rows = load_rows()
    reps = read_tsv(os.path.join(OUT, "fold_representatives.tsv"))
    fill = {r["slice"]: r["fill"] for r in read_tsv(os.path.join(OUT, "cath_folds.tsv"))}

    drawn, table = [], []
    for rep in reps:
        g = rep["slice"]
        row = rows[rep["rep_plasticome_id"]]
        seq = row["aa_sequence"].strip().upper()
        tag = f"{g}_{rep['representative'].replace(' ', '_').replace('/', '-')}"
        pdb_text, source, note = structure_for(rep, row, seq)
        chain, res, ident, path = best_chain(pdb_text, seq, tag)
        res, cut_n, cut_c = trim_predicted(res, source)
        ca = np.array([r["CA"].coord for r in res], dtype=float)
        ss = ss_from_records(path, chain, res) or ss_from_ca(ca)
        tint = fill[g]

        fig, ax = plt.subplots(figsize=(args.size, args.size))
        draw(ax, ca, ss, tint)
        fig.subplots_adjust(0, 0, 1, 1)
        for ext, dpi in (("svg", None), ("png", DPI)):
            fig.savefig(os.path.join(SDIR, ext, f"{tag}.{ext}"),
                        transparent=True, **({"dpi": dpi} if dpi else {}))
        plt.close(fig)

        n = collections.Counter(ss)
        drawn.append((rep, tag, ca, ss, tint))
        table.append([g, rep["representative"], rep["fold"], rep["catalytic_cath"],
                      tint, source, note, chain, f"{100 * ident:.1f}", len(res),
                      len(seq), n["H"], n["E"], n["L"],
                      "records" if ss_from_records(path, chain, res) else "ca_geometry",
                      cut_n, cut_c])
        print(f"  {g} {rep['representative']:<20} {source:<8} {note:<18} "
              f"chain {chain}  {len(res):>4} res  id {100 * ident:5.1f}%  "
              f"H{n['H']:>4} E{n['E']:>4}"
              + (f"  trimmed {cut_n}+{cut_c} dangling" if cut_n or cut_c else ""),
              flush=True)

    write_tsv(os.path.join(SDIR, "sources.tsv"),
              ["slice", "representative", "fold", "catalytic_cath", "fill",
               "source", "source_note", "chain", "percent_identity_to_rep",
               "residues_drawn", "rep_length_aa", "helix", "strand", "loop",
               "secondary_structure_from", "trimmed_n_term", "trimmed_c_term"], table)

    # --- contact sheet ------------------------------------------------------
    cols, size = 5, args.size
    rowsn = (len(drawn) + cols - 1) // cols
    fig, axes = plt.subplots(rowsn, cols, figsize=(cols * size, rowsn * size * 1.1))
    for ax in np.ravel(axes):
        ax.axis("off")
    for ax, (rep, tag, ca, ss, tint) in zip(np.ravel(axes), drawn):
        draw(ax, ca, ss, tint, lw_scale=0.72)
        ax.set_title(f"{rep['slice']}  {rep['representative']}\n"
                     f"{rep['fold'][:38]}\n{rep['catalytic_cath']}",
                     fontsize=7.5, color="#2b2b2b", linespacing=1.5, pad=6)
    fig.subplots_adjust(0.01, 0.01, 0.99, 0.95, wspace=0.05, hspace=0.3)
    for ext, dpi in (("svg", None), ("png", DPI)):
        fig.savefig(os.path.join(SDIR, f"folds_contact_sheet.{ext}"),
                    transparent=True, **({"dpi": dpi} if dpi else {}))
    plt.close(fig)

    print(f"\n{len(drawn)} structures -> {os.path.normpath(SDIR)}")
    print(f"  svg/ and png/ one file each, plus folds_contact_sheet.{{svg,png}}")
    print("  " + ", ".join(f"{k} {v}" for k, v in
                           collections.Counter(r[5] for r in table).items()))


if __name__ == "__main__":
    main()
