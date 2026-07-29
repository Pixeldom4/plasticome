#!/usr/bin/env python3
"""Collapse 'Enzyme N [like]' rows into their bare-number 'truth' rows.

Rule:
  - The row whose enzyme_name is (or contains) the bare number is the TRUTH.
  - Merge metadata (enzyme_name, pazy_id, accession, organism) with ';'
    (order-preserving, dedup, drop empties).
  - aa_sequence is NOT merged: keep the TRUTH row's sequence only.
  - The merged row keeps the TRUTH row's identifier.
"""
import csv, sys

SRC = "plasticomev1_join_retrievingfromaccession.tsv"
OUT = "plasticomev1_collapsed.tsv"

# (truth_identifier, enzyme_identifier, number)
GROUPS = [
    (35, 491, 101), (36, 486, 102), (37, 465, 202), (38, 482, 204),
    (39, 481, 211), (40, 466, 214), (41, 489, 301), (42, 496, 305),
    (44, 487, 307), (45, 49, 403),  (48, 463, 407), (50, 47, 412),
]

def merge(*vals):
    """Join ';'-delimited fields, order-preserving, dedup, drop empties."""
    seen, out = set(), []
    for v in vals:
        for tok in v.split(";"):
            tok = tok.strip()
            if tok and tok not in seen:
                seen.add(tok); out.append(tok)
    return ";".join(out)

with open(SRC, newline="") as fh:
    rows = list(csv.reader(fh, delimiter="\t"))
header, data = rows[0], rows[1:]
by_id = {int(r[0]): r for r in data}

log = []
enzyme_ids_removed = set()
for truth_id, enz_id, num in GROUPS:
    t, e = by_id[truth_id], by_id[enz_id]
    merged = [
        str(truth_id),                       # identifier: keep truth's
        merge(t[1], e[1]),                   # enzyme_name
        merge(t[2], e[2]),                   # pazy_id
        merge(t[3], e[3]),                   # accession
        merge(t[4], e[4]),                   # organism
        t[5],                                # aa_sequence: TRUTH ONLY
    ]
    by_id[truth_id] = merged
    enzyme_ids_removed.add(enz_id)
    log.append((num, truth_id, enz_id, merged, e[5]))

# emit all rows except the collapsed-away enzyme rows, sorted by identifier
out_rows = [by_id[i] for i in sorted(by_id) if i not in enzyme_ids_removed]
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(header)
    w.writerows(out_rows)

print(f"input records:  {len(data)}")
print(f"rows collapsed: {len(GROUPS)}")
print(f"output records: {len(out_rows)}")
print("\n# per-group result (identifier | enzyme_name | pazy | acc | org | seqlen | dropped_seqlen)")
for num, tid, eid, m, dropped in log:
    print(f"{num}: id={tid} name=[{m[1]}] pazy=[{m[2]}] acc=[{m[3]}] org=[{m[4]}] "
          f"kept_seqlen={len(m[5])} dropped(id{eid})_seqlen={len(dropped)}")
