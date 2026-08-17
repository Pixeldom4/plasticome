#!/usr/bin/env python3
"""Step 5 - the union annotated with its cluster and component facts, one row per sequence.

Steps 2 to 4 are cluster-level: 411 rows for 609 sequences, and the 198 rows that
are not a centroid appear only as a label inside `member_labels`. Step 5 turns
that back into a per-sequence table, so every union row states which cluster and
component it landed in, whether it is the centroid, and how close it is to one.

    01-union.tsv     609 rows, one per sequence, no cluster facts
    03-alignment.tsv 411 rows, one per cluster, membership as labels
    05-...tsv        609 rows, one per sequence, WITH cluster facts   <- this

Row order is the union's, unchanged, so `05` and `01` are row-aligned and can be
pasted side by side. The first five columns are `01-union.tsv` verbatim.

Columns
-------
    plasticome_id enzyme_name accession pazy_id source   from the union, verbatim
    seq_len                                              len of the normalized sequence
    cluster_id component_id                              from the step-3 table
    is_centroid                                          yes / no
    pct_id_to_centroid                                   as the aligner reported it
    centroid_identifier centroid_label centroid_accession  the centroid this row sits under
    aa_sequence                                          from the union, verbatim

`pct_id_to_centroid` is `100.0` for a centroid. The clusters table stores one
`member_pct_ids` entry per *non-representative* member and nothing for the
representative, since its identity to itself is not something the aligner
reports. `100.0` is this script's convention for that empty slot, not a measured
value, and it is exact rather than approximate: a centroid is its own sequence.

`aa_sequence` is passed through from the union **unnormalized**, unlike step 6,
which emits the normalized form because a hash is only meaningful over the bytes
that were hashed. Step 5 is the union with columns added, so the sequence column
stays the union's. On the current union the two are identical anyway, since step
1 already writes clean sequences; the distinction only matters if that changes.

Identifier resolution
---------------------
`centroid_identifier` is `PL` plus the `plasticome_id` of the union row the
centroid's `U####` label points at. That is the same identifier step 4 emits for
the same cluster, so `05.centroid_identifier` joins to `04.identifier` exactly.

It is resolved positionally rather than from `rep_plasticome_id`, which is blank
on all 411 rows of the current clusters table, and rather than from the
centroid's accession, which is the route that renamed union row 18 to `PL17`.
Membership comes from `lib/common/membership.py`, which verifies the positional
map is a partition of the union before returning it.

Example
-------
  python "lib/05 annotate/annotate_union.py" runs/<run>/03-alignment.tsv \\
      --union runs/<run>/01-union.tsv -o runs/<run>/05-union-with-components.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "03 alignment" / "scripts"))
from membership import LABEL_RE, membership  # noqa: E402
from step1_nodes import normalize  # noqa: E402

csv.field_size_limit(10 ** 7)

FROM_UNION = ["plasticome_id", "enzyme_name", "accession", "pazy_id", "source"]

COLS = FROM_UNION + [
    "seq_len", "cluster_id", "component_id", "is_centroid", "pct_id_to_centroid",
    "centroid_identifier", "centroid_label", "centroid_accession", "aa_sequence",
]

# A centroid's identity to itself. Not measured, and not the aligner's; see the
# module docstring.
SELF_PCT = "100.0"


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build(union: list[dict], clusters: list[dict], name: str, prefix: str) -> list[dict]:
    at = membership(clusters, len(union), name)

    # position -> the PL identifier of the union row there, which is what a
    # centroid label resolves to.
    ident_at = {i: f"{prefix}{(r.get('plasticome_id') or '').strip()}"
                for i, r in enumerate(union, 1)}

    out = []
    for pos, urow in enumerate(union, 1):
        m = at[pos]
        c = m["cluster"]
        rep_label = (c.get("rep_label") or "").strip()
        # membership() already proved every label parses and is in range.
        rep_pos = int(LABEL_RE.match(rep_label).group(1))
        row = {k: (urow.get(k) or "").strip() for k in FROM_UNION}
        row.update({
            "seq_len": len(normalize(urow.get("aa_sequence", ""))),
            "cluster_id": (c.get("cluster_id") or "").strip(),
            "component_id": (c.get("component_id") or "").strip(),
            "is_centroid": "yes" if m["is_centroid"] else "no",
            "pct_id_to_centroid": SELF_PCT if m["is_centroid"] else (m["pct_id"] or ""),
            "centroid_identifier": ident_at[rep_pos],
            "centroid_label": rep_label,
            "centroid_accession": (c.get("rep_accession") or "").strip(),
            "aa_sequence": urow.get("aa_sequence", ""),
        })
        out.append(row)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("clusters", type=Path,
                    help="Step-03 alignment TSV (Step-02 works, but leaves component_id blank)")
    ap.add_argument("--union", type=Path, default=None,
                    help="Step-01 union TSV (default: 01-union.tsv beside the input)")
    ap.add_argument("-o", "--out", type=Path, default=None,
                    help="output TSV (default: 05-union-with-components.tsv beside the input)")
    ap.add_argument("--prefix", default="PL", help="identifier prefix (default: PL)")
    args = ap.parse_args()

    if not args.clusters.exists():
        print(f"error: no such file: {args.clusters}", file=sys.stderr)
        return 1
    union_path = args.union or args.clusters.parent / "01-union.tsv"
    if not union_path.exists():
        print(f"error: no such union table: {union_path}", file=sys.stderr)
        return 1

    union = read_tsv(union_path)
    clusters = read_tsv(args.clusters)
    if "component_id" not in (clusters[0] if clusters else {}):
        print(f"note: {args.clusters.name} has no component_id; that column will be blank",
              file=sys.stderr)

    rows = build(union, clusters, args.clusters.name, args.prefix)

    out = args.out or args.clusters.parent / "05-union-with-components.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    n_cent = sum(1 for r in rows if r["is_centroid"] == "yes")
    print(f"{out}: {len(rows)} union rows over "
          f"{len({r['cluster_id'] for r in rows})} clusters / "
          f"{len({r['component_id'] for r in rows if r['component_id']})} components "
          f"({n_cent} centroid, {len(rows) - n_cent} not)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(0)
