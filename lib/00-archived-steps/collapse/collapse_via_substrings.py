#!/usr/bin/env python3
"""Collapse the erickson-collapsed dataset by substring groups.

Input:
  - erickson-collapsed TSV (the current dataset, 6 cols:
    identifier, enzyme_name, pazy_id, accession, organism, aa_sequence)
  - the curated substring-matches TSV (plasticome.v2 - substring-matches.tsv),
    whose `most-wildtype` column (TRUE/FALSE) marks, per substring group, the
    single member whose sequence is the true wildtype to keep.

What it does, per substring group:
  - KEEP the most-wildtype member's sequence.
  - MERGE the metadata (enzyme_name, pazy_id, accession, organism) of every
    group member into that kept row (';'-joined, order-preserving, dedup),
    keeping the kept row's `identifier`.
  - DROP the other members' rows.

Two robustness details, both real in this dataset:
  * The curated group numbers are unreliable — "group 5" bundles two distinct
    biological groups (UMG-SP-3 and Chath_Est1). So we ignore the curated group
    numbers and RE-CLUSTER the curated member sequences by exact-substring
    containment (union-find), which splits them correctly.
  * Members are matched to erickson rows by NORMALIZED sequence, not by id
    (the two files use different id spaces).

Special cases handled generically:
  * Identical members (UMG-SP-3: both curated variants normalize to the same
    439-aa sequence -> one erickson row) collapse to a no-op: one distinct
    sequence, already the wildtype, nothing dropped.
  * Sequence-swap (Enzyme 202): the wildtype (tag-stripped) sequence is NOT in
    the dataset; only the His-tagged form (id 110) is. Here we keep that row but
    REPLACE its sequence with the curated wildtype (strip the LEHHHHHH tag).

Output: substring-collapsed TSV in the same 6-col erickson schema.
"""
import csv
import os
import re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))          # lib/collapse/
REPO = os.path.dirname(os.path.dirname(HERE))              # repo root
# Data lives in a run, not next to this script. Point PLASTICOME_RUN at another
# runs/<date>_<name>/ to re-run this against a different run's collapse data.
RUN = os.environ.get("PLASTICOME_RUN",
                     os.path.join(REPO, "runs", "2026-07-20_collapse-v9-v10"))
ROOT = os.path.join(RUN, "collapse")
ERICK = os.path.join(ROOT, "erickson-collapsed",
                     "plasticomev1_join_retrievingfromaccession - plasticomev1_erickson_collapsed.tsv")
MATCH = os.path.join(ROOT, "substring-collapsed", "plasticome.v2 - substring-matches.tsv")
OUT = os.path.join(ROOT, "substring-collapsed", "plasticome.v2 - substring-collapsed.tsv")

META_COLS = [1, 2, 3, 4]  # enzyme_name, pazy_id, accession, organism


def norm(s):
    return re.sub(r"\s+", "", s or "").upper()


def merge(*vals):
    """Join ';'-delimited fields, order-preserving, dedup, drop empties."""
    seen, out = set(), []
    for v in vals:
        for tok in (v or "").split(";"):
            tok = tok.strip()
            if tok and tok not in seen:
                seen.add(tok)
                out.append(tok)
    return ";".join(out)


def read_tsv(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    # strip CR (the source files are CRLF) and surrounding whitespace
    return [[(c.replace("\r", "").strip() if c else c) for c in r] for r in rows]


def load_erickson():
    rows = read_tsv(ERICK)
    header = rows[0]
    data = [r for r in rows[1:] if any(c for c in r)]
    return header, data


def load_wildtypes():
    """Return list of (norm_seq, raw_seq, enzyme_name, is_wildtype) from curated file."""
    rows = read_tsv(MATCH)
    hdr = rows[0]
    gi = {c: i for i, c in enumerate(hdr)}
    out = []
    for r in rows[1:]:
        if not any(c for c in r):
            continue
        seq = r[gi["aa_sequence"]]
        if not seq:
            continue
        out.append({
            "seq": norm(seq),
            "raw": seq,
            "name": r[gi["enzyme_name"]],
            "wt": r[gi["most-wildtype"]].strip().upper() == "TRUE",
        })
    return out


def cluster_by_substring(members):
    """Union-find members whose normalized sequences are substring-related.

    Collapse identical sequences first (so UMG-SP-3's two identical variants
    become one node), then union pairs where one seq contains the other.
    Returns list of clusters; each cluster is a dict of
    {norm_seq: {'wt':bool, 'raw':str, 'names':[...]}}.
    """
    by_seq = {}
    for m in members:
        slot = by_seq.setdefault(m["seq"], {"wt": False, "raw": m["raw"], "names": []})
        slot["wt"] = slot["wt"] or m["wt"]
        slot["names"].append(m["name"])
    seqs = list(by_seq)
    parent = list(range(len(seqs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    order = sorted(range(len(seqs)), key=lambda i: len(seqs[i]))
    for a in range(len(order)):
        i = order[a]
        for b in range(a + 1, len(order)):
            j = order[b]
            if len(seqs[i]) == len(seqs[j]):
                continue
            if seqs[i] in seqs[j]:  # shorter contained in longer
                union(i, j)

    comp = defaultdict(dict)
    for idx, s in enumerate(seqs):
        comp[find(idx)][s] = by_seq[s]
    return list(comp.values())


def main():
    header, data = load_erickson()
    seq2rows = defaultdict(list)
    for r in data:
        seq2rows[norm(r[5])].append(r)

    clusters = cluster_by_substring(load_wildtypes())

    drop_ids = set()          # erickson identifiers to remove
    edits = {}                # identifier -> new row (merged / swapped)
    log = []

    for cl in clusters:
        wt_seqs = [s for s, v in cl.items() if v["wt"]]
        if len(wt_seqs) != 1:
            raise SystemExit(f"cluster has {len(wt_seqs)} wildtypes (expected 1): "
                             f"{[v['names'] for v in cl.values()]}")
        target = wt_seqs[0]

        # erickson rows for each member sequence
        member_rows = {s: seq2rows.get(s, []) for s in cl}
        present_ids = {r[0] for rows in member_rows.values() for r in rows}
        base_rows = member_rows.get(target, [])

        if not present_ids:
            log.append(("SKIP-absent", cl[target]["names"], "no member in dataset"))
            continue

        # distinct erickson rows across the whole cluster
        distinct = {r[0]: r for rows in member_rows.values() for r in rows}
        if len(distinct) == 1 and base_rows and norm(base_rows[0][5]) == target:
            log.append(("SKIP-noop", cl[target]["names"], f"single row id={next(iter(distinct))}"))
            continue

        # choose the base (kept) row: the wildtype row if present, else the
        # sole present row whose sequence we will swap to the wildtype.
        if base_rows:
            base = list(base_rows[0])
            swap = False
        else:
            # sequence-swap case (e.g. Enzyme 202): keep the present row, use
            # the curated wildtype sequence in place of its tagged sequence.
            base = list(distinct[sorted(distinct)[0]])
            base[5] = cl[target]["raw"]
            swap = True

        # merge metadata from base first, then the other members (id order)
        others = [distinct[i] for i in sorted(distinct) if i != base[0]]
        for c in META_COLS:
            base[c] = merge(base[c], *[o[c] for o in others])
        base[5] = base[5] if not swap else base[5]  # keep wildtype sequence

        edits[base[0]] = base
        drop_ids.update(i for i in distinct if i != base[0])
        log.append((("SWAP" if swap else "COLLAPSE"), cl[target]["names"],
                    f"keep id={base[0]} len={len(norm(base[5]))} "
                    f"drop={sorted(drop_ids & set(distinct))}"))

    # emit: apply edits, drop removed rows, preserve original order
    out_rows = []
    for r in data:
        if r[0] in drop_ids:
            continue
        out_rows.append(edits.get(r[0], r))

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(out_rows)

    print(f"input rows:   {len(data)}")
    print(f"dropped rows: {len(drop_ids)}")
    print(f"output rows:  {len(out_rows)}")
    print(f"wrote {OUT}\n")
    for action, names, detail in log:
        print(f"  {action:13} {';'.join(dict.fromkeys(names))[:34]:34} {detail}")


if __name__ == "__main__":
    main()
