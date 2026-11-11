#!/usr/bin/env python3
"""
Step 1: amino acid identity for every pair of sequences in the plasticome.

The SSN in ../ssn was built from a DIAMOND search, which only reports the
pairs it finds. A median over the pairs of two components needs every pair, so
this step aligns all of them: Smith-Waterman, BLOSUM62, gap open 11 and extend
1, the BLAST defaults.

Sequences are the distinct sequences the trees kept
(../tree/runs/all-sequences/intermediate/tips.tsv, status kept), so identical records count once and the one fragment is out.
The centroid set is a subset of these, so one run serves both sets.

Per pair:
  identity  identical columns over alignment columns, in percent. This is the
            pident of BLAST and DIAMOND, so it is on the scale of the SSN.
  coverage  alignment columns over the length of the shorter sequence, the
            definition ../ssn/build_ssn.py uses. Step 2 applies the SSN rule
            to it, 50 percent or more.

Writes work/pairs.tsv, one row per unordered pair.
"""

import itertools
import multiprocessing
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tree"))
from common import WORK_ALL, load_rows, read_tsv, write_tsv  # noqa: E402

WORK = os.path.join(HERE, "work")
GAP_OPEN, GAP_EXTEND = 11, 1

_aligner = None
_seq = None


def aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.mode = "local"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    # BLAST charges open + extend for the first gap column
    al.open_gap_score = -(GAP_OPEN + GAP_EXTEND)
    al.extend_gap_score = -GAP_EXTEND
    return al


def init(seq):
    global _aligner, _seq
    _aligner, _seq = aligner(), seq


def align(pair):
    a, b = pair
    sa, sb = _seq[a], _seq[b]
    aln = _aligner.align(sa, sb)[0]
    co = aln.coordinates
    ident = cols = 0
    for k in range(co.shape[1] - 1):
        a0, a1, b0, b1 = co[0, k], co[0, k + 1], co[1, k], co[1, k + 1]
        la, lb = a1 - a0, b1 - b0
        cols += max(la, lb)
        if la and lb:
            ident += sum(x == y for x, y in zip(sa[a0:a1], sb[b0:b1]))
    cov = cols / min(len(sa), len(sb))
    return [a, b, f"{100 * ident / cols:.1f}", ident, cols, f"{cov:.3f}",
            int(aln.score), len(sa), len(sb)]


def main():
    rows = load_rows()
    tips = sorted((t["tip"] for t in read_tsv(os.path.join(WORK_ALL, "tips.tsv"))
                   if t["status"] == "kept"), key=int)
    seq = {t: rows[t]["aa_sequence"].strip().upper() for t in tips}
    pairs = list(itertools.combinations(tips, 2))
    os.makedirs(WORK, exist_ok=True)
    with multiprocessing.Pool(os.cpu_count(), init, (seq,)) as pool:
        out = pool.map(align, pairs, chunksize=500)
    write_tsv(os.path.join(WORK, "pairs.tsv"),
              ["a", "b", "identity", "identical", "columns", "coverage_shorter",
               "score", "len_a", "len_b"], out)
    print(f"{len(tips)} sequences, {len(out)} pairs -> work/pairs.tsv")


if __name__ == "__main__":
    main()
