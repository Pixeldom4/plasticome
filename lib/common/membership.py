#!/usr/bin/env python3
"""Reconstruct cluster membership from a step-2/3 table, positionally.

`02-clusters.tsv` and `03-alignment.tsv` are one row per *cluster*, not per
sequence. They carry `rep_aa_sequence` for the centroid only, so there is no
sequence on that side to hash for a non-centroid member. The only handle on the
other members is `member_labels`, whose `U####` prefix is the union row's 1-based
position:

    U0018|WP_015787089.1   ->  union row 18

The accession after the pipe is there to keep labels unique and is not reliable
as a key. Cluster 3 of `2026-08-06.final-usearch.2` has `U0024|CBY05530` as its
representative and `U0514|CBY05530.1` as a member, the same accession at two
versions in two positions. Read the number, ignore the rest.

`rep_label` is NOT repeated inside `member_labels`, so the representative has to
be prepended to get the full membership. `member_pct_ids` is parallel to
`member_labels` and therefore has one fewer entry than the membership does; the
representative's identity to itself is not stored.

Both step 5 and step 6 need exactly this, so it lives here rather than in either
of them. The repo already carries one duplicated definition (`normalize()`, in
`lib/03 alignment/scripts/step1_nodes.py` and `lib/04 fasta/tsv_to_fasta.py`); a
second one would be a second way for two steps to silently disagree about which
union row a cluster contains.

The map is verified to be a **partition of the union** before it is returned.
Anything built on top of it, the md5-containment check in `06`, the per-row
annotation in `05`, is meaningless if a union row is in two clusters or in none,
so that is a hard failure rather than something to discover downstream.
"""

from __future__ import annotations

import re

# `U0018|WP_015787089.1`; the number is the union row's 1-based position. Same
# pattern clusters_to_fasta.py reads off `rep_label` for its route-2 resolution.
LABEL_RE = re.compile(r"^U0*(\d+)\|")


def labels_of(row: dict) -> list[str]:
    """A cluster's full membership: the representative first, then the members."""
    members = [x.strip() for x in (row.get("member_labels") or "").split(";")]
    return [(row.get("rep_label") or "").strip()] + [m for m in members if m]


def pct_ids_of(row: dict) -> list[str]:
    """`member_pct_ids`, parallel to `member_labels` and so one shorter than membership."""
    return [x.strip() for x in (row.get("member_pct_ids") or "").split(";") if x.strip()]


def membership(rows: list[dict], n_union: int, name: str) -> dict[int, dict]:
    """union position (1-based) -> {cluster, is_centroid, pct_id}, checked to be a partition.

    `pct_id` is the member's identity to its centroid as the aligner reported it,
    and None for the centroid itself, which has no stored value. Callers that
    want a number for the centroid supply their own convention rather than having
    one invented here.

    Raises SystemExit with every fault it found, not just the first: a table that
    is wrong in one way is usually wrong in several, and fixing them one run at a
    time is worse than seeing them together.
    """
    at: dict[int, dict] = {}
    doubled, unparsed, out_of_range, wrong_size, short_pct = [], [], [], [], []

    for row in rows:
        cid = (row.get("cluster_id") or "?").strip()
        labels = labels_of(row)
        pcts = pct_ids_of(row)
        size = (row.get("size") or "").strip()
        if size.isdigit() and int(size) != len(labels):
            wrong_size.append(f"{cid} (size {size}, {len(labels)} labels)")
        # One pct per non-representative member. A short list means the table was
        # edited without keeping the two columns in step.
        if len(pcts) < len(labels) - 1:
            short_pct.append(f"{cid} ({len(labels) - 1} members, {len(pcts)} pct_ids)")

        for i, label in enumerate(labels):
            m = LABEL_RE.match(label)
            if not m:
                unparsed.append(f"{cid}: {label!r}")
                continue
            pos = int(m.group(1))
            if not 1 <= pos <= n_union:
                out_of_range.append(f"{cid}: {label} (union has {n_union} rows)")
                continue
            if pos in at:
                doubled.append(f"U{pos:04d} in {at[pos]['cluster'].get('cluster_id')} and {cid}")
            at[pos] = {
                "cluster": row,
                "is_centroid": i == 0,
                "pct_id": None if i == 0 else (pcts[i - 1] if i - 1 < len(pcts) else None),
            }

    errs = []
    for what, bad in (("label(s) with no U#### prefix", unparsed),
                      ("label(s) pointing outside the union", out_of_range),
                      ("union row(s) in more than one cluster", doubled),
                      ("cluster(s) whose label count differs from `size`", wrong_size),
                      ("cluster(s) with fewer member_pct_ids than members", short_pct)):
        if bad:
            errs.append(f"{len(bad)} {what}: " + ", ".join(bad[:8])
                        + (" ..." if len(bad) > 8 else ""))
    uncovered = [p for p in range(1, n_union + 1) if p not in at]
    if uncovered:
        errs.append(f"{len(uncovered)} union row(s) in no cluster: "
                    + ", ".join(f"U{p:04d}" for p in uncovered[:8])
                    + (" ..." if len(uncovered) > 8 else ""))

    if errs:
        raise SystemExit(
            f"error: the labels in {name} are not a partition of the union.\n"
            + "".join(f"       {e}\n" for e in errs)
            + "       membership is reconstructed positionally from U#### labels; if that "
              "map is not a partition,\n       nothing built on top of it means anything."
        )
    return at
