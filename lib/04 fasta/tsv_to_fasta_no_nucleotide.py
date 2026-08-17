#!/usr/bin/env python3
"""Generate the plasticome FASTA, dropping nucleotide accessions from the header.

Same record format as tsv_to_fasta.py's `--header-template`:

    >petadex_id|accession_1;accession_2|<genbank always empty>|pazy_id

but each row's accession list is filtered to *protein / stable* identifiers only
-- GenBank/RefSeq **nucleotide** accessions are removed. The `accession` column in
cleaned_pazy_final.tsv mixes protein records (RefSeq WP_/NP_/XP_, UniProt, PDB,
MGnify, GenBank protein AAA#####/MBX#######) with the *nucleotide* entry the
protein was translated from (GenBank AB302136, OK558825, ...). For a protein
FASTA those nucleotide accessions are noise, so this script strips them.

Nucleotide classification (see is_nucleotide) is deliberately conservative:
  * RefSeq nucleotide prefixes  NC_/NG_/NM_/NR_/NT_/NW_/NZ_/XM_/XR_/AC_
  * GenBank nucleotide shapes   1 letter + 5 digits, 2 letters + 6 digits,
    2 letters + 8 digits (WGS) -- e.g. J01415, AB302136, OK558825
UniProt is checked first, because a Swiss-Prot id like P00590 has the same
1-letter+5-digit shape as a GenBank nucleotide id but is a *protein*. MGnify
(MGYP...) and PDB (leading digit) never match a nucleotide shape and are kept.

Sequence handling matches the rest of the repo: normalize() strips non-letters
and upper-cases, sequences wrap at 60 columns, blank-sequence rows are dropped.

Usage
-----
  python tsv_to_fasta_no_nucleotide.py cleaned_pazy_final.tsv
  python tsv_to_fasta_no_nucleotide.py in.tsv -o out.fasta
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


def duplicate_columns(fieldnames) -> list[str]:
    """Column names that appear more than once in the header row.

    csv.DictReader keeps only the *last* column of a repeated name, so a stray
    empty duplicate silently blanks every row: the v260701 sheet carries a
    second, near-empty `aa_sequence` at column 12, and reading it as the
    sequence yields 4 records out of 473 with a zero exit status. main()
    refuses to run rather than emit a truncated FASTA.
    """
    return [c for c, n in Counter(fieldnames or []).items() if n > 1]

# --- accession classification ------------------------------------------------

_UNIPROT = (re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]$'),
            re.compile(r'[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$'))
_REFSEQ_NT = re.compile(r'^(AC|NC|NG|NM|NR|NT|NW|NZ|XM|XR)_')
_GENBANK_NT = (re.compile(r'^[A-Z][0-9]{5}$'),      # 1 letter + 5 digits  (J01415)
               re.compile(r'^[A-Z]{2}[0-9]{6}$'),   # 2 letters + 6 digits (AB302136)
               re.compile(r'^[A-Z]{2}[0-9]{8}$'))   # 2 letters + 8 digits (WGS scaffold)


def is_uniprot(acc: str) -> bool:
    a = acc.split('.')[0].upper()
    return any(p.match(a) for p in _UNIPROT)


def is_nucleotide(acc: str) -> bool:
    """True for GenBank/RefSeq nucleotide accessions; False for protein/PDB/etc."""
    a = acc.split('.')[0].upper()          # drop the .version suffix
    if _REFSEQ_NT.match(a):
        return True
    if a.startswith('MGYP'):               # MGnify protein, keep
        return False
    if is_uniprot(a):                      # Swiss-Prot P00590 etc. shadow 1+5, keep
        return False
    return any(p.match(a) for p in _GENBANK_NT)


# --- sequence helpers (repo conventions) -------------------------------------

def normalize(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def wrap(seq: str, width: int = 60):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def protein_accessions(cell: str) -> str:
    """Return the accession cell with nucleotide tokens removed, ';'-joined."""
    kept = [t.strip() for t in (cell or "").split(";")
            if t.strip() and not is_nucleotide(t.strip())]
    return ";".join(kept)


class _Blank(dict):
    """Missing template keys format to '' instead of raising."""
    def __missing__(self, key):
        return ""


DEFAULT_TEMPLATE = "{identifier}|{accession}||{pazy_id}"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", type=Path, help="input tab-separated table")
    ap.add_argument("-o", "--out", default=None,
                    help="output FASTA path, or '-' for stdout "
                         "(default: <input>.protein.fasta)")
    ap.add_argument("--width", type=int, default=60, help="line-wrap width")
    ap.add_argument("--acc-col", default="accession",
                    help="column holding the ';'-separated accession list "
                         "(default: accession). The v1.1 sheet has no column "
                         "by that name -- point this at 'retrieved'. Whatever "
                         "the column is called, the filtered list is exposed "
                         "to the template as '{accession}'.")
    ap.add_argument("--header-template", default=DEFAULT_TEMPLATE,
                    help="Python str.format template over the row columns "
                         f"(default: '{DEFAULT_TEMPLATE}'). '{{accession}}' "
                         "resolves to the nucleotide-filtered accession list; "
                         "empty placeholders are kept, so literal separators "
                         "like the always-empty genbank pipe are preserved.")
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
        placeholders = re.findall(r"\{([^}:!]+)", args.header_template)
        # '{accession}' is synthesized from --acc-col, so it is not required to
        # exist under that name in the sheet.
        need = [args.acc_col, "aa_sequence",
                *(p for p in placeholders if p != "accession")]
        missing = [c for c in need if c not in cols]
        if missing:
            print(f"error: column(s) {missing} not in {args.tsv.name}; "
                  f"available: {cols}", file=sys.stderr)
            return 1
        rows = list(reader)

    if args.out == "-":
        out_fh, out_name, close = sys.stdout, "<stdout>", False
    else:
        out_path = (Path(args.out) if args.out
                    else args.tsv.with_suffix(".protein.fasta"))
        out_fh, out_name, close = out_path.open("w"), str(out_path), True

    written = blank = dropped_tokens = rows_emptied = 0
    try:
        for r in rows:
            seq = normalize(r.get("aa_sequence", ""))
            if not seq:
                blank += 1
                continue
            raw = [t.strip() for t in (r.get(args.acc_col) or "").split(";") if t.strip()]
            acc = protein_accessions(r.get(args.acc_col, ""))
            dropped_tokens += len(raw) - (len(acc.split(";")) if acc else 0)
            if raw and not acc:
                rows_emptied += 1
            # feed the template the filtered accession in place of the raw column
            vals = _Blank({k: (v or "").strip() for k, v in r.items()})
            vals["accession"] = acc
            header = args.header_template.format_map(vals)
            out_fh.write(f">{header}\n")
            for line in wrap(seq, args.width):
                out_fh.write(line + "\n")
            written += 1
    finally:
        if close:
            out_fh.close()

    print(f"{out_name}: wrote {written} records from {len(rows)} rows "
          f"({blank} blank seq; dropped {dropped_tokens} nucleotide accessions; "
          f"{rows_emptied} rows left with an empty accession field)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(0)
