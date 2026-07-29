#!/usr/bin/env python3
"""Adapter: centroids FASTA -> curated-TSV shape that step1_nodes.py consumes.

The v11 input (`plasticome.v1.1_centroids.fasta`) is already a de-duplicated set of
cluster centroids, so there is no md5-collapse to do -- but rather than reimplement
the frozen node-assembly / overlay logic, we simply re-express each centroid as one
TSV row with the columns step1_nodes.py reads, then run the unchanged pipeline.

FASTA header format (pipe-delimited, 4 fields):
    >PL<n>|<accession>|<enzyme_name>|<pazy_id>
Any field may be empty (e.g. `>PL116|||134`, `>PL504|||`). `organism` is not carried
in the FASTA, so it is emitted blank.
"""
import argparse
import csv
import re
from pathlib import Path

COLS = ["identifier", "enzyme_name", "pazy_id", "accession", "organism", "aa_sequence"]


def read_fasta(path):
    header, seq = None, []
    for line in Path(path).open():
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def parse_header(h):
    f = (h.split("|") + ["", "", "", ""])[:4]
    ident, accession, enzyme_name, pazy_id = (x.strip() for x in f)
    return {"identifier": ident, "accession": accession,
            "enzyme_name": enzyme_name, "pazy_id": pazy_id, "organism": ""}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for header, seq in read_fasta(args.fasta):
            row = parse_header(header)
            row["aa_sequence"] = re.sub(r"\s", "", seq)
            w.writerow(row)
            n += 1
    print(f"[centroids_to_tsv] {n} centroids -> {args.out}")


if __name__ == "__main__":
    main()
