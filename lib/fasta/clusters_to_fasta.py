#!/usr/bin/env python3
"""Emit the cluster-representative FASTA from an annotated clusters table.

Input is a Step-03 `alignment.tsv` (one row per cluster, carrying the centroid's
accession/pazy_id/sequence plus the connected-component id). Output is one record
per cluster with a five-field pipe header:

    >identifier|accession|genbank_accessions|pazy_id|component

  * `identifier`   -- a clean 1..N re-index over the sorted rows (see below),
                      prefixed (default `PL`). It is *not* carried over from the
                      union table's `plasticome_id`: the row set here is the 411
                      centroids, not the 607 union rows, so the numbering differs.
  * `accession`    -- the representative centroid's accession, i.e. the one the
                      sequence was pulled from. Blank for centroids with no
                      accession (placeholders like `jmPE13`).
  * `genbank_accessions` -- other accessions mapping to the same sequence.
                      Always emitted empty here; the pipe is kept so the field
                      count is fixed and the column can be filled in later.
  * `pazy_id`      -- blank on centroids that have none.
  * `component`    -- the alignment-graph component id (`C001`...).

Ordering (mirrors `lib/01 union/build_union.py:sort_key`, with the amino-acid
sequence replacing input order as the final tiebreak):

    1. pazy_id, numerically; rows *without* a pazy_id sort after every row that
       has one -- same rule build_union.py applies to the union table.
    2. component, numerically (`C012` -> 12).
    3. aa_sequence, lexically.

Because pazy_id is unique across the rows that have one, keys 2 and 3 only
discriminate among the blank-pazy_id tail.

Example
-------
  python lib/fasta/clusters_to_fasta.py \\
      runs/2026-07-30.plasticome.v1.212.union-spec/03-alignment.tsv \\
      -o runs/2026-07-30.plasticome.v1.212.union-spec/04-plasticome.v1.212.union-spec.fasta
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tsv_to_fasta import duplicate_columns, normalize, wrap  # noqa: E402

INF = float("inf")


def num(x: str) -> float:
    """Leading integer in a key, or +inf so blank/non-numeric sorts last.

    Handles both bare pazy ids (`16`) and prefixed component ids (`C012`); the
    union table's `num()` only ever saw the former.
    """
    m = re.search(r"\d+", x or "")
    return int(m.group()) if m else INF


def sort_key(row: dict, cols: dict) -> tuple:
    pazy = (row.get(cols["pazy"]) or "").strip()
    comp = (row.get(cols["component"]) or "").strip()
    seq = normalize(row.get(cols["seq"], ""))
    return (
        0 if pazy else 1,        # pazy_id rows first, blanks last
        num(pazy), pazy,
        num(comp), comp,
        seq,
    )


def records(rows, cols, prefix, pad):
    """Yield (header, sequence) in re-indexed order; skip blank-sequence rows."""
    keep = [r for r in rows if normalize(r.get(cols["seq"], ""))]
    keep.sort(key=lambda r: sort_key(r, cols))
    for i, r in enumerate(keep, 1):
        get = lambda c: (r.get(cols[c]) or "").strip()  # noqa: E731
        ident = f"{prefix}{i:0{pad}d}" if pad else f"{prefix}{i}"
        # genbank_accessions is intentionally empty -- the pipe is a placeholder.
        header = f"{ident}|{get('accession')}||{get('pazy')}|{get('component')}"
        yield header, normalize(r[cols["seq"]])


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", type=Path, help="Step-03 alignment/clusters TSV")
    ap.add_argument("-o", "--out", default=None,
                    help="output FASTA path, or '-' for stdout "
                         "(default: input name with a .fasta suffix)")
    ap.add_argument("--prefix", default="PL",
                    help="identifier prefix (default: PL)")
    ap.add_argument("--pad", type=int, default=0,
                    help="zero-pad the identifier to this width (default: 0, unpadded)")
    ap.add_argument("--width", type=int, default=0,
                    help="sequence line-wrap width; 0 = one line per record "
                         "(default: 0, matching the prior centroid FASTA)")
    ap.add_argument("--accession-col", default="rep_accession")
    ap.add_argument("--pazy-col", default="rep_pazy_id")
    ap.add_argument("--component-col", default="component_id")
    ap.add_argument("--seq-col", default="rep_aa_sequence")
    args = ap.parse_args()

    if not args.tsv.exists():
        print(f"error: no such file: {args.tsv}", file=sys.stderr)
        return 1

    cols = {"accession": args.accession_col, "pazy": args.pazy_col,
            "component": args.component_col, "seq": args.seq_col}

    with args.tsv.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        have = reader.fieldnames or []
        dupes = duplicate_columns(have)
        if dupes:
            print(f"error: duplicate column name(s) {dupes} in {args.tsv.name}; "
                  f"only the last of each is readable, which silently drops data.",
                  file=sys.stderr)
            return 1
        missing = [c for c in cols.values() if c not in have]
        if missing:
            print(f"error: column(s) {missing} not in {args.tsv.name}; "
                  f"available: {have}", file=sys.stderr)
            return 1
        rows = list(reader)

    if args.out == "-":
        out_fh, out_name, close = sys.stdout, "<stdout>", False
    else:
        out_path = Path(args.out) if args.out else args.tsv.with_suffix(".fasta")
        out_fh, out_name, close = out_path.open("w"), str(out_path), True

    n = 0
    try:
        for header, seq in records(rows, cols, args.prefix, args.pad):
            out_fh.write(f">{header}\n")
            if args.width > 0:
                for line in wrap(seq, args.width):
                    out_fh.write(line + "\n")
            else:
                out_fh.write(seq + "\n")
            n += 1
    finally:
        if close:
            out_fh.close()

    print(f"{out_name}: wrote {n} records from {len(rows)} rows "
          f"({len(rows) - n} blank seq)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(0)
