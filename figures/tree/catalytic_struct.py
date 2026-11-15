"""Catalytic residues placed by structure, for sequences too far from their
fold's annotated model to align by sequence (step 18's catalytic ring).

A sequence under 25 percent identity to the annotated model cannot be aligned
to it reliably, so its catalytic positions are found in two steps instead:

  1. the representative structure of the sequence's SSN component (step 3,
     intermediate/structures/cNN.pdb) is superposed on the annotated model
     with Foldseek in TM-align mode, which aligns by shape and holds up far
     below the identity where sequence alignment fails;
  2. the representative is aligned to the sequence by sequence. Members of a
     component are at least 30 percent identical to each other, so this step
     is the easy one.

The catalytic residue of the model maps through both alignments to a position
in the sequence. A superposition can slip by a residue or two in a loop, so a
residue counts as kept when the expected amino acid is at that position or
within two of it, and as lost when it is not. Where the position does not
align at all (a loop the superposition left unpaired) the residue is unknown,
and a sequence with an unknown residue gets no call unless another residue is
clearly lost. The call is only made when
the superposition is good (TM-score 0.5 or more, normalised by the model) and
the sequence is at least 25 percent identical to its representative (over
the shorter of the two);
otherwise the sequence is "not assessed".
"""

import os
import subprocess
import tempfile

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.SeqUtils import seq1

FOLDSEEK = "/opt/homebrew/Caskroom/miniconda/base/envs/plasticome-tree/bin/foldseek"
TM_MIN, ID_MIN = 0.5, 0.25

ALIGNER = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1,
                                end_gap_score=0)
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")


def residues(path, chain=None):
    """[(chain, number, one letter)] of the protein residues with a CA, in
    file order, the residues Foldseek reads."""
    model = next(iter(PDBParser(QUIET=True).get_structure("s", path)))
    out = []
    for ch in model:
        if chain and ch.id != chain:
            continue
        for r in ch:
            if r.id[0] == " " and "CA" in r:
                out.append((ch.id, r.id[1], seq1(r.get_resname())))
    return out


class _Chain(Select):
    def __init__(self, cid):
        self.cid = cid

    def accept_chain(self, chain):
        return chain.id == self.cid

    def accept_residue(self, r):
        return r.id[0] == " "


def _write_chain(src, cid, dst):
    s = PDBParser(QUIET=True).get_structure("s", src)
    io = PDBIO()
    io.set_structure(s)
    io.save(dst, _Chain(cid))


def superpose(pairs, workdir):
    """pairs: [(query pdb, target pdb, target chain)]. Returns
    {(query, target, chain): (tm, {target index: query index})}, indices into
    residues(query) and residues(target, chain)."""
    qdir, tdir = os.path.join(workdir, "q"), os.path.join(workdir, "t")
    os.makedirs(qdir, exist_ok=True)
    os.makedirs(tdir, exist_ok=True)
    names = {}
    for q, t, c in pairs:
        qn = os.path.basename(q).replace(".pdb", "")
        tn = f"{os.path.basename(t).replace('.pdb', '')}__{c}"
        if not os.path.exists(os.path.join(qdir, qn + ".pdb")):
            _write_chain(q, residues(q)[0][0], os.path.join(qdir, qn + ".pdb"))
        if not os.path.exists(os.path.join(tdir, tn + ".pdb")):
            _write_chain(t, c, os.path.join(tdir, tn + ".pdb"))
        names[(qn, tn)] = (q, t, c)
    out_m8 = os.path.join(workdir, "aln.m8")
    subprocess.run([FOLDSEEK, "easy-search", qdir, tdir, out_m8, os.path.join(workdir, "tmp"),
                    "--alignment-type", "1", "--exhaustive-search", "1", "-e", "inf",
                    "--format-output", "query,target,alntmscore,ttmscore,qstart,tstart,qaln,taln",
                    "-v", "1"], check=True)
    out = {}
    with open(out_m8) as fh:
        for line in fh:
            qn, tn, _, ttm, qs, ts, qa, ta = line.rstrip("\n").split("\t")
            # Foldseek names an entry by its file, sometimes with a chain suffix
            key = next(((q, t) for q, t in names if qn.startswith(q) and tn.startswith(t)), None)
            if key is None or names[key] in out:
                continue
            i, j, m = int(qs) - 1, int(ts) - 1, {}
            for a, b in zip(qa, ta):
                if a != "-" and b != "-":
                    m[j] = i
                if a != "-":
                    i += 1
                if b != "-":
                    j += 1
            out[names[key]] = (float(ttm), m)
    return out


def assess(items, workdir=None):
    """items: [(tip, tip sequence, representative pdb, reference pdb,
    [(chain, number, one letter)] wanted)]. Returns {tip: (category, tm,
    identity to representative, kept, total)}; category is one of "all
    catalytic residues kept", "some kept", "none kept", "not assessed"."""
    own = workdir is None
    workdir = workdir or tempfile.mkdtemp(prefix="catalytic_struct_")
    pairs = sorted({(rep, ref, c) for _, _, rep, ref, want in items if rep
                    for c, _, _ in want})
    sup = superpose(pairs, workdir) if pairs else {}
    out = {}
    for tip, seq, rep, ref, want in items:
        total = len(want)
        if not rep:
            out[tip] = ("not assessed", 0.0, 0.0, 0, total)
            continue
        rep_res = residues(rep)
        rep_seq = "".join(a for _, _, a in rep_res)
        aln = ALIGNER.align(rep_seq, seq)[0]
        pairs_ = [(int(i), int(j)) for i, j in zip(*aln.indices) if i >= 0 and j >= 0]
        ident = sum(rep_seq[i] == seq[j] for i, j in pairs_) / max(min(len(rep_seq), len(seq)), 1)
        to_tip = dict(pairs_)
        kept = lost = unknown = 0
        tms = []
        for c, num, aa in want:
            tm, m = sup.get((rep, ref, c), (0.0, {}))
            tms.append(tm)
            ref_res = residues(ref, c)
            idx = next((n for n, (_, rn, _) in enumerate(ref_res) if rn == num), None)
            q = m.get(idx) if idx is not None else None
            j = to_tip.get(q) if q is not None else None
            if j is None:
                unknown += 1
            elif any(0 <= j + d < len(seq) and seq[j + d] == aa for d in (0, -1, 1, -2, 2)):
                kept += 1
            else:
                lost += 1
        tm = min(tms) if tms else 0.0
        if tm < TM_MIN or ident < ID_MIN or (unknown and not lost):
            out[tip] = ("not assessed", tm, ident, kept, total)
        else:
            cat = ("all catalytic residues kept" if kept == total else
                   "some kept" if kept else "none kept")
            out[tip] = (cat, tm, ident, kept, total)
    if own:
        subprocess.run(["rm", "-rf", workdir])
    return out
