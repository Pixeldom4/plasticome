#!/usr/bin/env python3
"""Build plasticome_v1.v260701-union.tsv by folding v1.1 into v260701.

Join key: v260701 `accession` matched against v1.1 `retrieved` ONLY. No other v1
accession column (pazy_accession / pdb_accession / uniprot_accession) is used,
either for matching or for the output accession. v260701 is the base table;
unmatched v1.1 rows are appended.

Sequence rule
-------------
Sequences are always the ones retrieved from accession, never carried over from
v1.1 -- EXCEPT for the manual-assignment rows selected by is_manual(), whose
aa_sequence was set by hand in the v1 cleaning step and must survive.

  v260701-only  -> v260701 aa_sequence (already retrieved from accession)
  both (merged) -> v260701 aa_sequence, unless manually assigned -> v1.1's
  v1.1-only     -> BLANK, for fetch_sequences.py (Stage 2) to fill from
                   accession, unless manually assigned -> v1.1's

Leaving v1.1-only sequences blank here is deliberate: it is what hands those
rows to the Stage 2 accession lookup. Do not "helpfully" carry v1.1's sequence
forward, or Stage 2 becomes a no-op and the never-carried-over rule is silently
violated.

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
  python3 build_union.py <v1.1.tsv> <v260701.tsv> [-o out.tsv]

e.g.
  python3 "lib/01 union/build_union.py" source-data/plasticome.v1.1/plasticome.v1.1.tsv \\
      source-data/plasticome.v260701/cleaned_pazy-260701_retrieving_from_accession.tsv \\
      -o runs/<run>/01-union.tsv
"""
import argparse
import csv
from collections import Counter

OUT_COLS = ["plasticome_id", "enzyme_name", "accession", "pazy_id", "aa_sequence", "source"]
DEFAULT_OUT = "plasticome_v1.v260701-union.tsv"

SRC_V1 = "v1.1"
SRC_V2 = "v260701"
SRC_BOTH = "both"

# v1.1 `retrieved` values that are labels, not database accessions. jmPE13 has
# no database record at all (sequence curated from the paper supplement), so its
# name was parked in `retrieved` as a placeholder.
PLACEHOLDER_ACCESSIONS = {"jmPE13"}


def is_manual(enzyme_name):
    """True for rows whose aa_sequence was assigned by hand during v1 cleaning.

    Two groups, per the v1.1 cleaning step:
      - the Erickson primary/only enzymes, named `Enzyme <n>` / `Enzyme <n> like`,
        whose sequences were set to Erickson supplementary table D1 verbatim
        (they differ from the database record by signal sequence and/or His-tag,
        or are outright divergent);
      - jmPE13, carried forward from its paper supplement.

    Rows where Erickson is a *secondary* reference keep their primary name
    (`Est1; Enzyme 708`, `MtCut; Enzyme 606`, `RgCut-II`) and are NOT manual --
    hence the startswith test rather than a substring match.
    """
    n = (enzyme_name or "").strip()
    return n.startswith("Enzyme ") or n == "jmPE13"


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
    kept_manual = []  # merged rows where the manual v1.1 sequence beat v260701's

    # 1. Emit every v260701 row (base table).
    for idx, r in enumerate(v2_rows):
        a = norm(col(v2_header, r, "accession"))
        seq = col(v2_header, r, "aa_sequence", occurrence=0).strip()
        source = SRC_V2
        component = ""
        manual = False
        if a and a in v1_by_ret:
            source = SRC_BOTH
            for cand in v1_by_ret[a]:
                if id(cand) not in consumed:
                    consumed.add(id(cand))
                    component = col(v1_header, cand, "component")
                    v1_name = col(v1_header, cand, "enzyme_name")
                    v1_seq = col(v1_header, cand, "aa_sequence").strip()
                    # Manual assignments outrank the database sequence.
                    if is_manual(v1_name) and v1_seq:
                        if v1_seq != seq.strip():
                            kept_manual.append((v1_name, a, len(v1_seq), len(seq.strip())))
                        seq = v1_seq
                        manual = True
                    break
            # else: accession present in v1 but its row was already consumed by an
            # earlier duplicate v260701 accession -- still `both`, no component.
        out.append({
            "enzyme_name": col(v2_header, r, "enzyme_name"),
            "accession": norm_accession(col(v2_header, r, "accession")),
            "pazy_id": col(v2_header, r, "pazy_id"),
            "aa_sequence": seq,
            "source": source,
            "_component": component,
            "_manual": manual,
            "_idx": idx,
        })

    # 2. Append unmatched v1.1 rows.
    #    accession = `retrieved` only (may be blank); pazy_id blank.
    #    aa_sequence is left BLANK for Stage 2 to fill from the accession, except
    #    for manual assignments -- which includes every accession-less Erickson
    #    row, for which v1.1 is the only possible source.
    for idx, r in enumerate(v1_rows):
        if id(r) in consumed:
            continue
        name = col(v1_header, r, "enzyme_name")
        manual = is_manual(name)
        seq = col(v1_header, r, "aa_sequence").strip() if manual else ""
        out.append({
            "enzyme_name": name,
            "accession": norm_accession(col(v1_header, r, "retrieved")),
            "pazy_id": "",
            "aa_sequence": seq,
            "source": SRC_V1,
            "_component": col(v1_header, r, "component"),
            "_manual": manual,
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

    print(f"wrote   : {args.out}")
    counts = Counter(row["source"] for row in out)
    print(f"{SRC_V2:<8}: {counts[SRC_V2]}")
    print(f"{SRC_BOTH:<8}: {counts[SRC_BOTH]}")
    print(f"{SRC_V1:<8}: {counts[SRC_V1]}")
    print(f"total   : {len(out)}")

    manual_rows = [r for r in out if r["_manual"]]
    print(f"\nmanual-assignment rows (v1.1 sequence kept): {len(manual_rows)} "
          f"({sum(1 for r in manual_rows if r['aa_sequence'])} with a sequence)")
    if kept_manual:
        print(f"merged rows where the manual v1.1 sequence overrode v260701's: {len(kept_manual)}")
        for name, acc, l1, l2 in kept_manual:
            print(f"  {name:<20} {acc:<18} v1.1 len={l1:<5} (v260701 len={l2})")

    blank = [r for r in out if not r["aa_sequence"].strip()]
    no_acc = [r for r in blank if not r["accession"].strip()]
    print(f"\nrows needing Stage 2 (blank aa_sequence): {len(blank)}")
    print(f"  of which blank accession (Stage 2 CANNOT fill): {len(no_acc)}")
    for r in no_acc:
        print(f"    plasticome_id={r['plasticome_id']} {r['enzyme_name']} [{r['source']}]")


if __name__ == "__main__":
    main()
