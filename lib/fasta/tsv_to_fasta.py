#!/usr/bin/env python3
"""Generate a protein FASTA from a TSV sheet.

Every plasticome sheet (accession validation, clustering, the usearch
input) is a tab-separated table carrying an amino-acid column and one or more
identifier columns. The alignment/clustering tools downstream all speak FASTA,
so this is the one conversion that keeps getting re-implemented per pipeline
(step0_sanity.py, step1_nodes.py, ...). This is the shared, table-agnostic
version: point it at any TSV, name the id/sequence columns, get a FASTA.

Conventions match the rest of the repo:
  * normalize() strips every non-letter and upper-cases (same as the usearch
    steps), so whitespace/gaps/`*` stops never leak into a record;
  * sequences wrap at 60 columns;
  * rows whose sequence is blank after normalization are dropped, not emitted
    as empty records (they would break blastp/usearch downstream).

Defaults target runs/2026-07-20_collapse-v9-v10/accession-validation/cleaned_pazy_final-accession_validation.tsv
(id=identifier, seq=aa_sequence) but every column name is a flag.

Examples
--------
  # default columns, write beside the input
  python tsv_to_fasta.py runs/2026-07-20_collapse-v9-v10/accession-validation/cleaned_pazy_final-accession_validation.tsv

  # header ">accession pazy_id | organism", dedup identical sequences
  python tsv_to_fasta.py in.tsv -o out.fasta \\
      --id-col accession --desc-cols pazy_id organism --dedup

  # pipe to stdout
  python tsv_to_fasta.py in.tsv -o -
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import Counter
from pathlib import Path


def duplicate_columns(fieldnames) -> list[str]:
    """Column names that appear more than once in the header row.

    csv.DictReader keeps only the *last* column of a repeated name, so a stray
    empty duplicate silently blanks every row: the v260701 sheet carries a
    second, near-empty `aa_sequence` at column 12, and reading it as the
    sequence yields 4 records out of 473 with a zero exit status. Callers
    refuse to run rather than emit a truncated FASTA.
    """
    return [c for c, n in Counter(fieldnames or []).items() if n > 1]


def normalize(seq: str) -> str:
    """Strip non-letters and upper-case — the repo-wide sequence normalizer."""
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def wrap(seq: str, width: int = 60):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


class _Blank(dict):
    """Missing keys format to '' so a template never raises on an absent column."""
    def __missing__(self, key):
        return ""


def build_header(row: dict, id_col: str, desc_cols: list[str], sep: str,
                 template: str | None) -> str:
    if template is not None:
        # Fixed layout: every placeholder is emitted verbatim, empty or not, so
        # literal separators (e.g. the always-empty genbank pipe in
        # "{identifier}|{accession}||{pazy_id}") are preserved.
        vals = _Blank({k: (v or "").strip() for k, v in row.items()})
        return template.format_map(vals)
    ident = (row.get(id_col) or "").strip()
    if not ident:
        raise KeyError(id_col)
    extras = [(row.get(c) or "").strip() for c in desc_cols]
    extras = [e for e in extras if e]
    return f"{ident}{(sep + sep.join(extras)) if extras else ''}"


def convert(rows, out_fh, id_col, desc_cols, seq_col, sep, width, dedup,
            template=None) -> dict:
    """Write FASTA records to out_fh; return skip/dup/write counts."""
    stats = {"total": 0, "written": 0, "blank": 0, "dup": 0, "no_id": 0}
    seen: set[str] = set()  # md5(seq), only populated when dedup
    for row in rows:
        stats["total"] += 1
        seq = normalize(row.get(seq_col, ""))
        if not seq:
            stats["blank"] += 1
            continue
        try:
            header = build_header(row, id_col, desc_cols, sep, template)
        except KeyError:
            stats["no_id"] += 1
            continue
        if dedup:
            digest = hashlib.md5(seq.encode()).hexdigest()
            if digest in seen:
                stats["dup"] += 1
                continue
            seen.add(digest)
        out_fh.write(f">{header}\n")
        for line in wrap(seq, width):
            out_fh.write(line + "\n")
        stats["written"] += 1
    return stats


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", type=Path, help="input tab-separated table")
    ap.add_argument("-o", "--out", default=None,
                    help="output FASTA path, or '-' for stdout "
                         "(default: input name with a .fasta suffix)")
    ap.add_argument("--id-col", default="identifier",
                    help="column used as the FASTA record id (default: identifier)")
    ap.add_argument("--desc-cols", nargs="*", default=[],
                    help="extra columns appended to the header after --sep")
    ap.add_argument("--seq-col", default="aa_sequence",
                    help="amino-acid column (default: aa_sequence)")
    ap.add_argument("--sep", default=" ",
                    help="separator between header fields (default: single space)")
    ap.add_argument("--width", type=int, default=60,
                    help="sequence line-wrap width (default: 60)")
    ap.add_argument("--dedup", action="store_true",
                    help="emit only the first record per identical sequence (md5)")
    ap.add_argument("--header-template", default=None,
                    help="Python str.format template over the row columns, e.g. "
                         "'{identifier}|{accession}||{pazy_id}'. Placeholders are "
                         "emitted verbatim (empty fields kept), so literal "
                         "separators like an always-empty pipe are preserved. "
                         "Overrides --id-col/--desc-cols/--sep.")
    args = ap.parse_args()

    if not args.tsv.exists():
        print(f"error: no such file: {args.tsv}", file=sys.stderr)
        return 1

    with args.tsv.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        dupes = duplicate_columns(cols)
        if dupes:
            print(f"error: duplicate column name(s) {dupes} in {args.tsv.name}; "
                  f"only the last of each is readable, which silently drops "
                  f"data. Rename or remove the duplicate column(s).",
                  file=sys.stderr)
            return 1
        if args.header_template is not None:
            needed = re.findall(r"\{([^}:!]+)", args.header_template)
            missing = [c for c in [args.seq_col, *needed] if c not in cols]
        else:
            missing = [c for c in [args.id_col, args.seq_col, *args.desc_cols]
                       if c not in cols]
        if missing:
            print(f"error: column(s) {missing} not in {args.tsv.name}; "
                  f"available: {cols}", file=sys.stderr)
            return 1
        rows = list(reader)

    if args.out == "-":
        out_fh, out_name, close = sys.stdout, "<stdout>", False
    else:
        out_path = Path(args.out) if args.out else args.tsv.with_suffix(".fasta")
        out_fh, out_name, close = out_path.open("w"), str(out_path), True

    try:
        s = convert(rows, out_fh, args.id_col, args.desc_cols, args.seq_col,
                    args.sep, args.width, args.dedup, args.header_template)
    finally:
        if close:
            out_fh.close()

    print(f"{out_name}: wrote {s['written']} records from {s['total']} rows "
          f"({s['blank']} blank seq, {s['no_id']} missing id"
          f"{f', {s['dup']} dup collapsed' if args.dedup else ''})",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        # downstream closed early (e.g. `... | head`); exit quietly
        sys.stderr.close()
        raise SystemExit(0)
