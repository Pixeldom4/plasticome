#!/usr/bin/env python3
"""Build plasticome.v1-v2-union.tsv by folding v1 into v2.

Join key: v2 `accession` matched against v1 `retrieved` ONLY. No other v1
accession column (pazy_accession / pdb_accession / uniprot_accession) is used,
either for matching or for the output accession. v2 is the base table; unmatched
v1 rows are appended in the 5-column output schema.

Deterministic: no network, no external state. Re-running regenerates the output
identically.
"""
import csv
from collections import Counter

HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))
V1 = HERE + "/plasticome.v1.1.tsv"
V2 = HERE + "/cleaned_pazy-260701_retrieving_from_accession.tsv"
OUT = HERE + "/plasticome.v1-v2-union.tsv"

OUT_COLS = ["enzyme_name", "accession", "pazy_id", "aa_sequence", "source"]


def norm(x):
    x = (x or "").strip()
    return "" if x in ("", "-") else x


def read_rows(path):
    """Return list of rows as list-of-(header,value) preserving duplicate columns."""
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


def main():
    v1_header, v1_rows = read_rows(V1)
    v2_header, v2_rows = read_rows(V2)

    # Index v1 rows by their `retrieved` accession (the ONLY join field).
    v1_by_ret = {}
    for r in v1_rows:
        k = norm(col(v1_header, r, "retrieved"))
        if k:
            v1_by_ret.setdefault(k, []).append(r)

    consumed = set()  # id() of v1 rows folded into a v2 row
    out = []

    # 1. Emit every v2 row (base table).
    for r in v2_rows:
        a = norm(col(v2_header, r, "accession"))
        source = "v2"
        if a and a in v1_by_ret:
            for cand in v1_by_ret[a]:
                if id(cand) not in consumed:
                    consumed.add(id(cand))
                    source = "merged"
                    break
            else:
                source = "merged"  # accession present in v1 but row already consumed
        out.append({
            "enzyme_name": col(v2_header, r, "enzyme_name"),
            "accession": col(v2_header, r, "accession"),
            "pazy_id": col(v2_header, r, "pazy_id"),
            "aa_sequence": col(v2_header, r, "aa_sequence", occurrence=0),  # primary
            "source": source,
        })

    # 2. Append unmatched v1 rows, mapped into the 5-column schema.
    #    accession = `retrieved` only (may be blank); pazy_id blank.
    #    aa_sequence is carried forward from v1 (the ONLY sequence source for
    #    accession-less rows such as the Erickson "Enzyme 1xx like" enzymes).
    for r in v1_rows:
        if id(r) in consumed:
            continue
        out.append({
            "enzyme_name": col(v1_header, r, "enzyme_name"),
            "accession": norm(col(v1_header, r, "retrieved")),
            "pazy_id": "",
            "aa_sequence": col(v1_header, r, "aa_sequence").strip(),
            "source": "v1",
        })

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    counts = Counter(row["source"] for row in out)
    print(f"v2     : {counts['v2']}")
    print(f"merged : {counts['merged']}")
    print(f"v1     : {counts['v1']}")
    print(f"total  : {len(out)}")
    blank = sum(1 for row in out if row["source"] == "v1" and not row["accession"])
    print(f"v1-only rows with blank accession: {blank}")


if __name__ == "__main__":
    main()
