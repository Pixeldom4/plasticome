#!/usr/bin/env python3
"""Apply the accession-ordering fixes flagged by validate_accession_order.py.

Two kinds of fix:

  * REORDER  -- reorder tokens by type-rank (protein -> cross-ref -> PDB),
                stable within a rank. No identifier is changed, added, or removed.
  * EDIT     -- an explicit, hand-verified content change (typo / dup / junk token).

Every EDIT is verified and documented inline below. Rows that need a human
decision (row 458) are deliberately left untouched and reported at the end.

Writes <tsv>.bak once, rewrites the file in place, prints a before->after diff,
then re-runs the validator.

    python fix_accession_order.py            # fix cleaned_pazy_final.tsv
    python fix_accession_order.py --dry-run  # show changes, write nothing
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys

from validate_accession_order import classify, RANK

# Rows to reorder by type-rank (no identifier changes).
REORDER_IDS = {"19", "32", "101", "118", "130", "131", "132"}

# Explicit content edits, keyed by identifier: (before, after, reason).
EDITS = {
    "58": ("BAC67195.1;Q83VDO", "BAC67195.1;Q83VD0",
           "typo O->0: UniProt Q83VD0's EMBL cross-ref is exactly BAC67195.1"),
    "76": ("BAC67242.1;AB102945.1;Q874E;Q874E9", "BAC67242.1;AB102945.1;Q874E9",
           "dropped 'Q874E' -- a truncated duplicate of 'Q874E9' in the same field"),
    "479": ("Chen, 2024;G0EX11", "G0EX11",
            "dropped 'Chen, 2024' -- a literal citation, not an accession"),
}

# Rows intentionally left for human review (reported, never edited).
NEEDS_REVIEW = {
    "458": ("N20M5AZM016;A4Y035;2FX5",
            "A4Y035 resolves to Ectopseudomonas mendocina (=ABP86951.1), but the row "
            "is Halopseudomonas aestusnigri PE-H -- organism mismatch; and "
            "'N20M5AZM016' matches no known accession format. Needs the correct "
            "primary NCBI protein accession for PE-H."),
}


def reorder(field: str) -> str:
    """Stable-sort tokens by type-rank; unknown-rank tokens keep their position."""
    toks = [t.strip() for t in field.split(";") if t.strip()]
    ranked = []
    for i, t in enumerate(toks):
        r = RANK.get(classify(t))
        # Unknown -> use a fractional rank tied to current index so it doesn't jump.
        key = r if r is not None else (toks_rank_fallback(toks, i))
        ranked.append((key, i, t))
    ranked.sort(key=lambda x: (x[0], x[1]))
    return ";".join(t for _, _, t in ranked)


def toks_rank_fallback(toks, i):
    """Keep an unrecognized token adjacent to its current neighbours."""
    return i / max(len(toks), 1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", nargs="?", default="cleaned_pazy_final.tsv")
    ap.add_argument("--acc-col", default="accession")
    ap.add_argument("--dry-run", action="store_true", help="print changes, write nothing")
    args = ap.parse_args()

    with open(args.tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)

    changes: list[tuple[str, str, str, str]] = []  # id, before, after, reason
    for row in rows:
        ident = row.get("identifier", "?")
        before = row[args.acc_col]

        if ident in EDITS:
            exp_before, after, reason = EDITS[ident]
            if before != exp_before:
                sys.exit(f"ABORT: row {ident} accession is {before!r}, "
                         f"expected {exp_before!r} -- data changed since verification.")
            row[args.acc_col] = after
            changes.append((ident, before, after, reason))
        elif ident in REORDER_IDS:
            after = reorder(before)
            if after != before:
                changes.append((ident, before, after, "reordered by type-rank"))
                row[args.acc_col] = after

    # --- report ---
    print(f"{'DRY-RUN: ' if args.dry_run else ''}{len(changes)} row(s) changed in {args.tsv}\n")
    for ident, before, after, reason in changes:
        print(f"  row {ident}: {reason}")
        print(f"      - {before}")
        print(f"      + {after}\n")

    if NEEDS_REVIEW:
        print("Left for manual review (NOT changed):")
        for ident, (field, why) in NEEDS_REVIEW.items():
            print(f"  row {ident}: {field}")
            print(f"      -> {why}\n")

    if args.dry_run:
        print("dry-run: no files written.")
        return

    backup = args.tsv + ".bak"
    shutil.copyfile(args.tsv, backup)
    print(f"backup written: {backup}")

    with open(args.tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"rewrote: {args.tsv}\n")

    print("re-running validator:")
    subprocess.run([sys.executable, "validate_accession_order.py", args.tsv])


if __name__ == "__main__":
    main()
