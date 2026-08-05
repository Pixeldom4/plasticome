#!/usr/bin/env python3
"""Build plasticome_v1.v260701-union.tsv by folding v1.1 into v260701.

Join key: v260701 `accession` matched against v1.1 `retrieved` ONLY. No other v1
accession column (pazy_accession / pdb_accession / uniprot_accession) is used,
either for matching or for the output accession. v260701 is the base table;
unmatched v1.1 rows are appended.

Sequence rule
-------------
The sequence retrieved from the accession wins wherever retrieval is possible.
The v1.1 table encodes that policy in its own data: it carries an aa_sequence
only on rows Stage 2 could never fill -- those whose `retrieved` is blank, or is
a placeholder rather than a database accession. So the rule here is just:

  v260701-only  -> v260701 aa_sequence, per --v2-seq (below)
  both (merged) -> same
  v1.1-only     -> v1.1's aa_sequence if it has one, else BLANK for
                   fetch_sequences.py (Stage 2) to fill from the accession

The `both` case cannot collide with a carried sequence: `retrieved` is the join
key, and a v1.1 row that carries a sequence has none, so it can never merge. A
merged row that does carry one means the input table has broken that invariant,
so it is reported rather than silently dropped.

This was previously a name test -- is_manual(): `Enzyme <n>` (Erickson 2022) plus
jmPE13 -- which kept v1.1's hand-curated sequence even for rows with a working
accession. It was dropped because the name was never a sound proxy for
provenance: the `Enzyme <n>` numbering is not Erickson-exclusive (cf.
`CtPL; Enzyme 504`, Avilan 2023), and only some Erickson enzymes are absent from
NCBI -- the ones that are present are better taken from their accession. Curate
the input table instead: blank a sequence to have it fetched, keep one to
protect it.

--v2-seq picks what the v260701 aa_sequence column means:

  as-given        the column already holds sequences retrieved from accession,
                  so take it verbatim. Correct for
                  cleaned_pazy-260701_retrieving_from_accession.tsv.
  from-accession  the column holds the sequence attached to the PAZy record,
                  which is authoritative only where there is no accession to
                  retrieve from. Keep it on accession-less rows; blank it
                  everywhere else so Stage 2 fetches from the accession.
                  Correct for cleaned_pazy-260701-singletons.tsv.

A blanked sequence is not discarded: it is written to v2_attached_sequences.tsv
beside the union, and Stage 2 restores it for any accession the databases cannot
resolve. That keeps "retrieved from accession wins" true wherever retrieval is
actually possible, without silently dropping a sequence curated by hand from a
paper for an accession that nothing serves (e.g. MGYP000321434903).

"Accession-less" means blank after norm_accession(), so a placeholder accession
that is dropped from the output (PLACEHOLDER_ACCESSIONS) also keeps its attached
sequence -- Stage 2 would have nothing to look it up by.

Ordering
--------
`plasticome_id` is 1..N over rows sorted on pazy_id (primary) then component
(secondary). pazy_id exists only on v260701-derived rows and component only on
v1.1-derived rows, so in practice: pazy_id rows first in numeric pazy_id order,
then the rest in numeric component order, ties broken by input order.

Deterministic: no network, no external state. Re-running regenerates the output
identically.

Usage
-----
  python3 build_union.py <v1.1.tsv> <v260701.tsv> [-o out.tsv] [--v2-seq MODE]

e.g.
  python3 "lib/01 union/build_union.py" source-data/plasticome.v1.1/plasticome.v1.1.no-seq.tsv \\
      source-data/plasticome.v260701/cleaned_pazy-260701-singletons.tsv \\
      --v2-seq from-accession -o runs/<run>/01-union.tsv
"""
import argparse
import csv
import os
from collections import Counter

OUT_COLS = ["plasticome_id", "enzyme_name", "accession", "pazy_id", "aa_sequence", "source"]
DEFAULT_OUT = "plasticome_v1.v260701-union.tsv"
ATTACHED_SIDECAR = "v2_attached_sequences.tsv"

SRC_V1 = "v1.1"
SRC_V2 = "v260701"
SRC_BOTH = "both"

# v1.1 `retrieved` values that are labels, not database accessions. jmPE13 has
# no database record at all (sequence curated from the paper supplement), so its
# name was parked in `retrieved` as a placeholder.
PLACEHOLDER_ACCESSIONS = {"jmPE13"}


def norm(x):
    x = (x or "").strip()
    return "" if x in ("", "-") else x


def norm_accession(x):
    """Accession as it should appear in the output: placeholders become blank."""
    x = norm(x)
    return "" if x in PLACEHOLDER_ACCESSIONS else x


def read_rows(path):
    """Return (header, rows) with rows as raw lists, preserving duplicate columns."""
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        return header, [row for row in reader]


def col(header, row, name, occurrence=0):
    """Value of the `occurrence`-th column named `name` (handles dup headers)."""
    seen = -1
    for i, h in enumerate(header):
        if h == name:
            seen += 1
            if seen == occurrence:
                return row[i] if i < len(row) else ""
    return ""


def num(x):
    """Numeric sort value for pazy_id / component; non-numeric and blank sort last."""
    x = (x or "").strip()
    return int(x) if x.isdigit() else float("inf")


def sort_key(entry):
    """pazy_id primary, component secondary, input order as the final tiebreak."""
    pazy = (entry["pazy_id"] or "").strip()
    comp = (entry["_component"] or "").strip()
    return (
        0 if pazy else 1,      # pazy_id rows first
        num(pazy), pazy,
        num(comp), comp,
        entry["_idx"],
    )


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("v1", help="v1.1 table (joined on its `retrieved` column)")
    p.add_argument("v2", help="v260701 table, the base (joined on its `accession` column)")
    p.add_argument("-o", "--out", default=DEFAULT_OUT, help=f"output TSV (default: ./{DEFAULT_OUT})")
    p.add_argument("--v2-seq", choices=("as-given", "from-accession"), default="as-given",
                   help="what the v260701 aa_sequence column holds: sequences already "
                        "retrieved from accession (as-given, the default), or the "
                        "PAZy-attached sequence, kept only where the row has no "
                        "accession and blanked for Stage 2 otherwise (from-accession)")
    return p.parse_args()


def main():
    args = parse_args()
    v1_header, v1_rows = read_rows(args.v1)
    v2_header, v2_rows = read_rows(args.v2)

    # Index v1 rows by their `retrieved` accession (the ONLY join field).
    v1_by_ret = {}
    for r in v1_rows:
        k = norm(col(v1_header, r, "retrieved"))
        if k:
            v1_by_ret.setdefault(k, []).append(r)

    consumed = set()  # id() of v1 rows folded into a v260701 row
    out = []
    merged_with_seq = []  # v1.1 rows that merged yet carry a sequence -- see docstring
    attached = []         # (accession, enzyme_name, seq) blanked under from-accession

    # 1. Emit every v260701 row (base table).
    for idx, r in enumerate(v2_rows):
        a = norm(col(v2_header, r, "accession"))
        out_acc = norm_accession(col(v2_header, r, "accession"))
        seq = col(v2_header, r, "aa_sequence", occurrence=0).strip()
        # Under from-accession the column is the PAZy-attached sequence: it stands
        # only where there is no accession for Stage 2 to retrieve from. What is
        # blanked here is handed to Stage 2 as a fallback, not thrown away.
        if args.v2_seq == "from-accession" and out_acc:
            if seq:
                attached.append((out_acc, col(v2_header, r, "enzyme_name"), seq))
            seq = ""
        source = SRC_V2
        component = ""
        if a and a in v1_by_ret:
            source = SRC_BOTH
            for cand in v1_by_ret[a]:
                if id(cand) not in consumed:
                    consumed.add(id(cand))
                    component = col(v1_header, cand, "component")
                    # A merged v1.1 row should never carry a sequence: it merged,
                    # so it has an accession, so Stage 2 can fetch it. Report the
                    # breach instead of dropping the sequence without a word.
                    v1_seq = col(v1_header, cand, "aa_sequence").strip()
                    if v1_seq:
                        merged_with_seq.append(
                            (col(v1_header, cand, "enzyme_name"), a, len(v1_seq)))
                    break
            # else: accession present in v1 but its row was already consumed by an
            # earlier duplicate v260701 accession -- still `both`, no component.
        out.append({
            "enzyme_name": col(v2_header, r, "enzyme_name"),
            "accession": out_acc,
            "pazy_id": col(v2_header, r, "pazy_id"),
            "aa_sequence": seq,
            "source": source,
            "_component": component,
            "_idx": idx,
        })

    # 2. Append unmatched v1.1 rows.
    #    accession = `retrieved` only (may be blank); pazy_id blank.
    #    aa_sequence is taken as given: v1.1 carries one only where Stage 2 could
    #    not fetch it (blank or placeholder `retrieved`), and blank everywhere
    #    else is exactly what hands the row to the Stage 2 accession lookup.
    for idx, r in enumerate(v1_rows):
        if id(r) in consumed:
            continue
        out.append({
            "enzyme_name": col(v1_header, r, "enzyme_name"),
            "accession": norm_accession(col(v1_header, r, "retrieved")),
            "pazy_id": "",
            "aa_sequence": col(v1_header, r, "aa_sequence").strip(),
            "source": SRC_V1,
            "_component": col(v1_header, r, "component"),
            "_idx": len(v2_rows) + idx,
        })

    # 3. Order and number.
    out.sort(key=sort_key)
    for i, row in enumerate(out, start=1):
        row["plasticome_id"] = i

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    # Sequences blanked above, for Stage 2 to fall back on. Written even when
    # empty so a stale sidecar from a previous build can never be picked up.
    side = os.path.join(os.path.dirname(os.path.abspath(args.out)), ATTACHED_SIDECAR)
    with open(side, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "enzyme_name", "aa_sequence"])
        w.writerows(attached)

    print(f"wrote   : {args.out}")
    print(f"v2 seq  : {args.v2_seq}")
    counts = Counter(row["source"] for row in out)
    print(f"{SRC_V2:<8}: {counts[SRC_V2]}")
    print(f"{SRC_BOTH:<8}: {counts[SRC_BOTH]}")
    print(f"{SRC_V1:<8}: {counts[SRC_V1]}")
    print(f"total   : {len(out)}")

    carried = [r for r in out if r["source"] == SRC_V1 and r["aa_sequence"].strip()]
    print(f"\nv1.1 sequences carried through (rows Stage 2 cannot fetch): {len(carried)}")
    if merged_with_seq:
        print(f"WARNING: {len(merged_with_seq)} merged v1.1 row(s) carry a sequence, which is "
              f"dropped in favour of the accession -- v1.1 should not hold one here:")
        for name, acc, n in merged_with_seq:
            print(f"  {name:<20} {acc:<18} v1.1 len={n}")

    if attached:
        print(f"\nv260701 attached sequences held for Stage 2 fallback: {len(attached)} "
              f"-> {os.path.basename(side)}")
        for acc, name, seq in attached:
            print(f"    {name:<20} {acc:<18} len={len(seq)}")

    blank = [r for r in out if not r["aa_sequence"].strip()]
    no_acc = [r for r in blank if not r["accession"].strip()]
    print(f"\nrows needing Stage 2 (blank aa_sequence): {len(blank)}")
    print(f"  of which blank accession (Stage 2 CANNOT fill): {len(no_acc)}")
    for r in no_acc:
        print(f"    plasticome_id={r['plasticome_id']} {r['enzyme_name']} [{r['source']}]")


if __name__ == "__main__":
    main()
