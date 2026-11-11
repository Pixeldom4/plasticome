#!/usr/bin/env python3
"""
Step 1: Pfam domains on every PETai sequence that carries an activity label.

Runs in the venv of ../../../petadex-adjacency-landscape, which has pyhmmer and
the Pfam-A HMMs (data/Pfam-A.hmm.gz); the plasticome env has neither.

    ~/Documents/projects/petadex-adjacency-landscape/.venv/bin/python 01_scan_pfam.py

Every domain hit at independent E-value 1e-5 or better is kept. Step 2 decides
which hit is the catalytic domain and which are accessory.

Writes work/petai_pfam.tsv (seq_md5, length, pfam_acc, pfam_name, start, end,
i_evalue, score).
"""

import os
import time

import pyhmmer

HERE = os.path.dirname(os.path.abspath(__file__))
LAND = os.path.join(HERE, "..", "..", "..", "petadex-adjacency-landscape", "data")
SEQS = os.path.join(LAND, "petai_sequences.fasta")
PFAM = os.path.join(LAND, "Pfam-A.hmm.gz")
OUT = os.path.join(HERE, "work", "petai_pfam.tsv")
EVALUE = 1e-5


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(SEQS, digital=True, alphabet=alpha) as fh:
        seqs = fh.read_block()
    length = {s.name.decode() if isinstance(s.name, bytes) else s.name: len(s) for s in seqs}
    t0 = time.time()
    n = 0
    with pyhmmer.plan7.HMMFile(PFAM) as hmms, open(OUT, "w") as out:
        out.write("seq_md5\tlength\tpfam_acc\tpfam_name\tstart\tend\ti_evalue\tscore\n")
        for hits in pyhmmer.hmmsearch(hmms, seqs, cpus=8, E=1e-3, domE=EVALUE):
            q = hits.query
            acc = q.accession.decode() if isinstance(q.accession, bytes) else q.accession
            name = q.name.decode() if isinstance(q.name, bytes) else q.name
            for h in hits:
                if not h.included:
                    continue
                sid = h.name.decode() if isinstance(h.name, bytes) else h.name
                for d in h.domains:
                    if d.i_evalue > EVALUE:
                        continue
                    out.write(f"{sid}\t{length[sid]}\t{acc}\t{name}\t{d.env_from}\t{d.env_to}\t"
                              f"{d.i_evalue:.2e}\t{d.score:.1f}\n")
                    n += 1
    print(f"{len(seqs)} sequences, {n} domain hits, {time.time() - t0:.0f} s -> {OUT}")


if __name__ == "__main__":
    main()
