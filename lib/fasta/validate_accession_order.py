#!/usr/bin/env python3
"""Validate the ordering convention of the semicolon-separated `accession` column.

Convention (see "Plasticome v2 Overview"): each accession field lists identifiers
in *identifier-type priority* order, not alphabetically:

    1. primary protein-sequence accession   -- NCBI protein preferred:
                                                RefSeq (WP_/NP_/XP_/YP_/AP_) else
                                                GenBank protein (e.g. CCK74972.1)
    2. cross-references                      -- UniProt acc / entry name, nucleotide
                                                (DNA) accession, MGnify, locus tags
    3. PDB structure codes                   -- always last (e.g. 5XFY, 7Z6B)

The check assigns each token a type-rank (protein=0, cross-ref=1, PDB=2) and
verifies the rank sequence is non-decreasing. It also flags rows whose first
token is not the primary protein accession when a protein accession exists later.

Usage:
    python validate_accession_order.py [TSV ...]        # defaults to cleaned_pazy_final.tsv
    python validate_accession_order.py --acc-col accession --strict
Exit status is non-zero if any exceptions are found.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter

# --- token classification --------------------------------------------------

REFSEQ_PREFIXES = ("WP_", "NP_", "XP_", "YP_", "AP_")

# Official UniProtKB accession pattern.
UNIPROT_ACC = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)


def classify(tok: str) -> str:
    """Return a coarse identifier type for one accession token."""
    t = tok.strip()
    if not t:
        return "empty"
    up = t.upper()

    # PDB: 4-char code, starts with a digit (5XFY, 7Z6B, 501l), optionally with a
    # chain suffix (5KXVA, 7CUV_A, 5ZOA_A, 7NEI_A/B).
    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}([_]?[A-Za-z](/[A-Za-z])?)?", t):
        return "PDB"
    # RefSeq protein.
    if up.startswith(REFSEQ_PREFIXES):
        return "RefSeq"
    # MGnify predicted protein.
    if up.startswith("MGYP"):
        return "MGnify"
    # UniProt accession (P19833, A0A0K8P6T7, Q6A0I3, ...).
    if UNIPROT_ACC.fullmatch(t):
        return "UniProt"
    # UniProt entry name or NCBI locus tag (has an underscore: R4YKL9_OLEAN, A4W93_05950).
    if "_" in t and not up.startswith(REFSEQ_PREFIXES):
        return "UniProt/locus"
    # GenBank protein accession: 3 letters + >=5 digits, optional .version (CCK74972.1).
    if re.fullmatch(r"[A-Z]{3}[0-9]{5,}(\.[0-9]+)?", t):
        return "GenBank-protein"
    # GenBank nucleotide / other short accession (OK558825, AJ810119).
    if re.fullmatch(r"[A-Z]{1,2}[0-9]{5,6}(\.[0-9]+)?", t):
        return "Nucleotide"
    # Bare locus tag: letters + a short run of digits, no underscore (RPA1511).
    # (Real protein/nucleotide accessions have >=5 digits and are caught above.)
    if re.fullmatch(r"[A-Z]{2,5}[0-9]{3,4}", up):
        return "Locus"
    return "UNKNOWN"


# type-rank: protein source (0) -> cross-ref (1) -> structure (2)
RANK = {
    "RefSeq": 0,
    "GenBank-protein": 0,
    "MGnify": 0,          # protein-sequence source when no NCBI protein exists
    "UniProt": 1,
    "UniProt/locus": 1,
    "Nucleotide": 1,
    "Locus": 1,
    "PDB": 2,
}
PROTEIN_TYPES = {"RefSeq", "GenBank-protein"}


def validate_row(field: str) -> list[str]:
    """Return a list of human-readable problems with one accession field (empty = ok)."""
    toks = [x.strip() for x in field.split(";") if x.strip()]
    if not toks:
        return []  # blank accession is not an ordering problem (handled elsewhere)

    kinds = [classify(t) for t in toks]
    problems: list[str] = []

    # (a) unrecognized tokens
    for tok, kind in zip(toks, kinds):
        if kind == "UNKNOWN":
            problems.append(f"unrecognized token '{tok}'")

    # (b) type-rank must be non-decreasing (this enforces "PDB last")
    ranks = [RANK.get(k) for k in kinds]
    known = [(t, k, r) for t, k, r in zip(toks, kinds, ranks) if r is not None]
    for (t0, k0, r0), (t1, k1, r1) in zip(known, known[1:]):
        if r1 < r0:
            problems.append(
                f"out-of-order: {k1} '{t1}' follows {k0} '{t0}' "
                f"(rank {r1} < {r0})"
            )

    # (c) first token should be the primary protein accession when one exists
    if kinds[0] not in PROTEIN_TYPES:
        later_protein = next(
            (t for t, k in zip(toks[1:], kinds[1:]) if k in PROTEIN_TYPES), None
        )
        if later_protein:
            problems.append(
                f"primary misplaced: field starts with {kinds[0]} '{toks[0]}' "
                f"but protein accession '{later_protein}' appears later"
            )

    return problems


def run(path: str, acc_col: str) -> tuple[int, int, list[tuple]]:
    total = flagged = 0
    exceptions: list[tuple] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if acc_col not in reader.fieldnames:
            sys.exit(f"error: column '{acc_col}' not in {path} "
                     f"(have: {', '.join(reader.fieldnames)})")
        id_col = "identifier" if "identifier" in reader.fieldnames else reader.fieldnames[0]
        for row in reader:
            field = row[acc_col]
            if not field.strip():
                continue
            total += 1
            problems = validate_row(field)
            if problems:
                flagged += 1
                exceptions.append((row.get(id_col, "?"), field, problems))
    return total, flagged, exceptions


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", nargs="*", default=["cleaned_pazy_final.tsv"],
                    help="TSV file(s) to check (default: cleaned_pazy_final.tsv)")
    ap.add_argument("--acc-col", default="accession",
                    help="name of the accession column (default: accession)")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any exceptions are found")
    args = ap.parse_args()

    any_flagged = False
    for path in args.tsv:
        total, flagged, exceptions = run(path, args.acc_col)
        print(f"\n=== {path} ===")
        print(f"checked {total} non-empty '{args.acc_col}' fields; "
              f"{flagged} exception(s)")
        by_type = Counter()
        for ident, field, problems in exceptions:
            print(f"\n  row {ident}: {field}")
            for p in problems:
                print(f"      - {p}")
                by_type[p.split(':')[0].split('(')[0].strip()] += 1
        if exceptions:
            print("\n  summary by kind:")
            for kind, n in by_type.most_common():
                print(f"      {n:4d}  {kind}")
        any_flagged = any_flagged or bool(flagged)

    if args.strict and any_flagged:
        sys.exit(1)


if __name__ == "__main__":
    main()
