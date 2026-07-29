#!/usr/bin/env python3
"""Convert the cleaned PAZy TSV into a FASTA using the Plasticome (PL) header format.

Header format (``|`` = delimiter, ``;`` = sub-delimiter for multi-valued fields):

    >PL<n>|<accession(s)>|<genbank_accessions>|<pazy_id>

Field mapping:
  * PL<n>             : the ``identifier`` column (already prefixed, e.g. PL1)
  * accession(s)      : the ``accession`` column (input accessions), ``,`` -> ``;``
  * genbank_accessions: left empty (accessions mapped via sequence search go here)
  * pazy_id           : the ``pazy_id`` column
"""

import argparse
import csv
import sys

SUBDELIM = ";"


def norm_accessions(value):
    """Split a comma-separated accession list, trim, and re-join with ';'."""
    parts = [p.strip() for p in (value or "").replace(";", ",").split(",")]
    return SUBDELIM.join(p for p in parts if p)


def build_header(identifier, accession, genbank, pazy_id):
    return "{ident}|{acc}|{gb}|{pazy}".format(
        ident=(identifier or "").strip(),
        acc=norm_accessions(accession),
        gb=genbank,  # intentionally empty
        pazy=(pazy_id or "").strip(),
    )


def tsv_to_fasta(in_path, out_path, id_col, acc_col, pazy_col, seq_col):
    written = 0
    with open(in_path, newline="") as fin, open(out_path, "w") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        for field in (id_col, acc_col, pazy_col, seq_col):
            if field not in reader.fieldnames:
                sys.exit(
                    f"error: column '{field}' not found. "
                    f"available: {reader.fieldnames}"
                )
        for i, row in enumerate(reader, start=2):  # line 1 is the header
            seq = (row[seq_col] or "").strip()
            header = build_header(row[id_col], row[acc_col], "", row[pazy_col])
            if not (row[id_col] or "").strip():
                print(f"warning: skipping line {i} (no identifier)", file=sys.stderr)
                continue
            if not seq:
                print(f"warning: skipping line {i} (empty sequence) [{header}]",
                      file=sys.stderr)
                continue
            fout.write(f">{header}\n{seq}\n")
            written += 1
    return written


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="input TSV file")
    p.add_argument("output", help="output FASTA file")
    p.add_argument("--id-col", default="identifier",
                   help="identifier column, e.g. PL1 (default: identifier)")
    p.add_argument("--acc-col", default="accession",
                   help="input accession column (default: accession)")
    p.add_argument("--pazy-col", default="pazy_id",
                   help="column for pazy_id (default: pazy_id)")
    p.add_argument("--seq-col", default="aa_sequence",
                   help="sequence column (default: aa_sequence)")
    args = p.parse_args()

    n = tsv_to_fasta(args.input, args.output, args.id_col, args.acc_col,
                     args.pazy_col, args.seq_col)
    print(f"wrote {n} records to {args.output}")


if __name__ == "__main__":
    main()
