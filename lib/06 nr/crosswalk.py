#!/usr/bin/env python3
"""Step 6, optional half - the B x C crosswalk, sequence md5 to cluster and component.

`build_nr.py` deliberately puts no `cluster_id`, `component_id` or `is_centroid`
in `06-nr.tsv`. Those are branch-B facts, and branch B depends on an unpinned
`--engine` choice. C's single most valuable property is that it is the one set in
the pipeline whose row count does not move when the engine does, so the branch-B
facts live here instead, in a file that carries the provenance of the branch it
was built against.

This script is therefore optional by construction. It exits 0 with a note when
branch B is absent from the run directory, which is what makes `--only 6` work in
a directory holding nothing but `01-union.tsv`.

Output
------
`06-nr-to-clusters.tsv`, one row per distinct sequence, in the same order as
`06-nr.tsv`:

    seq_md5  cluster_id  component_id  is_centroid  centroid_identifier  engine

`component_id` is blank when only `02-clusters.tsv` is available, since the
component partition is a step-3 fact.

`engine` repeats on every row. It is constant within a run, and carrying it keeps
the file self-describing once its rows have been joined onto something else --
`cluster_id` and `component_id` mean nothing without knowing which engine drew
them, and diamond and usearch runs of the same union are both on disk here.

The clustering **identity** is not a column. It was one, and it was always `0.9`:
`--id` is a knob nobody has moved in seven runs, so a per-row copy of it bought
nothing that `engine` was not already buying. It lives in
`06-nr.intermediates/crosswalk.json` instead, together with the rest of the
provenance this file used to be the only record of.

Two joins, not one
------------------
**The centroid is matched on `seq_md5`.** Both sides carry the sequence, so both
get hashed and the hashes are compared. This bypasses the `rep_plasticome_id` /
`rep_label` / accession resolution cascade in `clusters_to_fasta.py` entirely,
including the accession route that renamed union row 18 to `PL17`. It is the
route that produces `centroid_identifier`, which is the field a caller is most
likely to trust.

**Membership is matched on label position.** `02-clusters.tsv` and
`03-alignment.tsv` carry `rep_aa_sequence` for the centroid only, never for the
other members, so there is no sequence on the B side to hash for a non-centroid
row. Membership has to come from `rep_label` plus `member_labels`, whose `U####`
prefix is the union row's 1-based position. That is the same reconstruction the
containment measurement of 2026-08-17 used. It lives in `lib/common/membership.py`
because step 5 needs exactly the same thing, and that module checks the positional
map is a partition of the union before returning it: every union row in exactly
one cluster, no row in two, label count equal to the cluster's stated `size`.

Checks
------
All hard failures. None of them are cosmetic.

  * **md5 containment.** Every `seq_md5` group falls inside exactly one
    `cluster_id` and one `component_id`. A violation means the clustering split
    identical sequences, which invalidates B rather than merely relabelling it.
    This properly belongs as a step-2 postcondition, so that it fails where the
    fault is introduced rather than three steps downstream in an optional file;
    it is re-asserted here because it is free. Measured 0 of 109 groups spanning
    more than one cluster or component on `2026-08-06.final-usearch.2`.
  * **Centroid md5 uniqueness.** The 411 centroids of B must already be
    md5-distinct. Step 3 reported 411 md5-unique nodes from 411 clusters, so this
    holds trivially, but it is what makes the B-to-C join one-to-one on the
    centroid side, so it is asserted rather than assumed.
  * **Coverage.** Every row of `06-nr.tsv` gets a cluster, and every centroid's
    md5 is found in `06-nr.tsv`.

Example
-------
  python "lib/06 nr/crosswalk.py" runs/2026-08-06.final-usearch.2/06-nr.tsv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "03 alignment" / "scripts"))
try:
    from build_nr import write_fasta  # noqa: E402
    from membership import membership  # noqa: E402
    from step1_nodes import normalize  # noqa: E402
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        f"error: cannot import normalize() from 'lib/03 alignment/scripts/step1_nodes.py' "
        f"or membership() from 'lib/common/membership.py': {exc}"
    )

csv.field_size_limit(10 ** 7)

COLS = ["seq_md5", "cluster_id", "component_id", "is_centroid",
        "centroid_identifier", "engine"]


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def md5_of(seq: str) -> str:
    return hashlib.md5(normalize(seq).encode()).hexdigest()


def find_b(run_dir: Path) -> tuple[Path | None, bool]:
    """(branch-B table, has_components). Step 3 preferred, since only it has them."""
    for name, comps in (("03-alignment.tsv", True), ("02-clusters.tsv", False)):
        p = run_dir / name
        if p.exists():
            return p, comps
    for pat, comps in (("03-*.tsv", True), ("02-*.tsv", False)):
        hits = sorted(run_dir.glob(pat))
        if len(hits) == 1:
            return hits[0], comps
    return None, False


def provenance(run_dir: Path) -> tuple[str, str]:
    """(engine, identity) as step 2 recorded them, else ('', '')."""
    p = run_dir / "02-clusters.intermediates" / "provenance.json"
    try:
        with p.open() as fh:
            d = json.load(fh)
    except Exception:
        return "", ""
    return str(d.get("engine_version") or d.get("engine") or ""), str(d.get("id") or "")


def write_provenance(path: Path, rows: list[dict], nr_path: Path, b_name: str,
                     engine: str, identity: str, has_comps: bool) -> dict:
    """What branch B this crosswalk was built against.

    A separate file from `build_nr.py`'s `provenance.json` rather than a block
    inside it: the two are written by different scripts at different times, and
    re-running `build_nr.py` would silently drop a block it does not own.
    """
    n_cent = sum(1 for r in rows if r["is_centroid"] == "yes")
    prov = {
        "step": 6,
        "part": "crosswalk",
        "built_against": b_name,
        "engine": engine,
        "identity": identity,
        "nr_input": nr_path.name,
        "n_sequences": len(rows),
        "n_clusters": len({r["cluster_id"] for r in rows}),
        "n_components": len({r["component_id"] for r in rows if r["component_id"]}) if has_comps else None,
        "n_centroid": n_cent,
        "n_non_centroid": len(rows) - n_cent,
        # The file only gets written when these passed, so recording them is a
        # statement about this build, not a claim re-checked at read time.
        "checks_passed": ["union partition", "md5 containment",
                          "centroid md5 uniqueness", "coverage"],
    }
    with path.open("w") as fh:
        json.dump(prov, fh, indent=2)
        fh.write("\n")
    return prov


def build(nr: list[dict], union: list[dict], b_rows: list[dict], b_name: str,
          has_comps: bool, prefix: str) -> list[dict]:
    at = membership(b_rows, len(union), b_name)

    # md5 join, centroid side. This is the route that yields centroid_identifier.
    centroid_md5: dict[str, dict] = {}
    collisions = []
    for row in b_rows:
        h = md5_of(row.get("rep_aa_sequence", ""))
        if h in centroid_md5:
            collisions.append(f"{centroid_md5[h].get('cluster_id')} and {row.get('cluster_id')}")
        centroid_md5[h] = row
    if collisions:
        raise SystemExit(
            f"error: {len(collisions)} pair(s) of clusters in {b_name} have the same "
            f"centroid sequence: " + ", ".join(collisions[:8])
            + (" ..." if len(collisions) > 8 else "")
            + "\n       the B-to-C join is one-to-one on the centroid side only if the "
              "centroids are md5-distinct."
        )

    by_cluster = {(row.get("cluster_id") or "").strip(): row for row in b_rows}
    # md5 -> the rep_plasticome_id C chose for that sequence, so the identifier
    # comes from the union row rather than from B's own labelling.
    nr_rep = {r["seq_md5"]: (r.get("rep_plasticome_id") or "").strip() for r in nr}

    missing_centroids = [(row.get("cluster_id") or "?") for h, row in centroid_md5.items()
                         if h not in nr_rep]
    if missing_centroids:
        raise SystemExit(
            f"error: {len(missing_centroids)} centroid sequence(s) in {b_name} are not in "
            f"06-nr.tsv: " + ", ".join(missing_centroids[:8])
            + (" ..." if len(missing_centroids) > 8 else "")
            + "\n       every centroid is a union row, so its sequence must appear in C; "
              "the two were built from different unions."
        )

    # Membership side: which cluster(s) does each seq_md5 land in.
    seen: dict[str, set] = defaultdict(set)
    for pos, m in at.items():
        h = md5_of(union[pos - 1].get("aa_sequence", ""))
        c = m["cluster"]
        seen[h].add((c.get("cluster_id") or "", (c.get("component_id") or "") if has_comps else ""))

    split = sorted(h for h, s in seen.items()
                   if len({c for c, _ in s}) > 1 or len({p for _, p in s}) > 1)
    if split:
        detail = []
        for h in split[:8]:
            where = ", ".join(f"{c}/{p or '-'}" for c, p in sorted(seen[h]))
            detail.append(f"{h[:12]} ({nr_rep.get(h, '?')}) -> {where}")
        raise SystemExit(
            f"error: {len(split)} md5 group(s) span more than one cluster or component:\n"
            + "".join(f"       {d}\n" for d in detail)
            + ("       ...\n" if len(split) > 8 else "")
            + "       identical sequences were split by the clustering, which invalidates B "
              "rather than relabelling it.\n"
              "       On usearch this is a -maxrejects early-exit, not a threshold effect; "
              "rerun step 2 with -maxrejects 0."
        )

    uncovered = [r["seq_md5"] for r in nr if r["seq_md5"] not in seen]
    if uncovered:
        raise SystemExit(
            f"error: {len(uncovered)} sequence(s) in 06-nr.tsv appear in no cluster of "
            f"{b_name}: " + ", ".join(h[:12] for h in uncovered[:8])
            + (" ..." if len(uncovered) > 8 else "")
        )

    out = []
    for r in nr:                      # 06-nr.tsv order, so the two files stay row-aligned
        h = r["seq_md5"]
        cid, comp = next(iter(seen[h]))          # containment: the set is a singleton
        crep = md5_of(by_cluster[cid].get("rep_aa_sequence", ""))
        out.append({
            "seq_md5": h,
            "cluster_id": cid,
            "component_id": comp,
            "is_centroid": "yes" if crep == h else "no",
            "centroid_identifier": f"{prefix}{nr_rep[crep]}",
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("nr", type=Path, help="06-nr.tsv from build_nr.py")
    ap.add_argument("-o", "--out", type=Path, default=None,
                    help="output TSV (default: 06-nr-to-clusters.tsv beside the input)")
    ap.add_argument("--union", type=Path, default=None,
                    help="Step-01 union TSV (default: 01-union.tsv beside the input)")
    ap.add_argument("--clusters", type=Path, default=None,
                    help="branch-B table (default: 03-alignment.tsv, else 02-clusters.tsv, "
                         "beside the input)")
    ap.add_argument("--engine", default=None,
                    help="override the engine recorded by step 2")
    ap.add_argument("--id", default=None,
                    help="override the clustering identity recorded by step 2")
    ap.add_argument("--intermediates", type=Path, default=None,
                    help="sidecar directory for crosswalk.json "
                         "(default: 06-nr.intermediates beside the input)")
    ap.add_argument("--fasta", type=Path, default=None,
                    help="06-nr.fasta to re-emit with the component filled into header "
                         "field 5 (default: the .fasta beside the input)")
    ap.add_argument("--no-fasta", action="store_true",
                    help="leave 06-nr.fasta as build_nr.py wrote it, field 5 empty")
    ap.add_argument("--width", type=int, default=0,
                    help="FASTA line-wrap width; 0 = one line per record (default: 0)")
    ap.add_argument("--prefix", default="PL",
                    help="identifier prefix for centroid_identifier (default: PL)")
    args = ap.parse_args()

    if not args.nr.exists():
        print(f"error: no such file: {args.nr}", file=sys.stderr)
        return 1
    run_dir = args.nr.parent

    union_path = args.union or run_dir / "01-union.tsv"
    if not union_path.exists():
        print(f"error: no such union table: {union_path}", file=sys.stderr)
        return 1

    if args.clusters:
        b_path, has_comps = args.clusters, "03" in args.clusters.name
        if not b_path.exists():
            print(f"error: no such clusters table: {b_path}", file=sys.stderr)
            return 1
    else:
        b_path, has_comps = find_b(run_dir)

    # Optional by construction: branch B may simply not exist in this run directory.
    if b_path is None:
        print(f"skip: no 02/03 table in {run_dir.name}; the crosswalk needs branch B "
              f"and 06 does not.", file=sys.stderr)
        return 0

    nr = read_tsv(args.nr)
    union = read_tsv(union_path)
    b_rows = read_tsv(b_path)

    rows = build(nr, union, b_rows, b_path.name, has_comps, args.prefix)

    engine, identity = provenance(run_dir)
    engine = args.engine if args.engine is not None else engine
    identity = args.id if args.id is not None else identity
    for r in rows:
        r["engine"] = engine

    out = args.out or run_dir / "06-nr-to-clusters.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    inter = args.intermediates or run_dir / "06-nr.intermediates"
    inter.mkdir(parents=True, exist_ok=True)
    write_provenance(inter / "crosswalk.json", rows, args.nr, b_path.name,
                     engine, identity, has_comps)

    # Field 5 of the FASTA header is the component, which is a branch-B fact and so
    # unknown to build_nr.py. It writes the file with field 5 empty; this fills it.
    # The file therefore has two byte-states under one name -- see --no-fasta.
    fasta = args.fasta or args.nr.with_suffix(".fasta")
    if not args.no_fasta and has_comps and fasta.exists():
        # rep_pazy_id is header field 4 and deliberately not a column of
        # 06-nr.tsv, so it is re-derived from the union rather than read back.
        pazy = {(r.get("plasticome_id") or "").strip(): (r.get("pazy_id") or "").strip()
                for r in union}
        groups = read_tsv(args.nr)
        for g in groups:
            g["rep_pazy_id"] = pazy.get(g["rep_plasticome_id"], "")
        comps = {r["seq_md5"]: r["component_id"] for r in rows}
        n = write_fasta(fasta, groups, args.width, comps)
        print(f"{fasta}: field 5 filled on {n} of {len(groups)} records", file=sys.stderr)

    n_cent = sum(1 for r in rows if r["is_centroid"] == "yes")
    print(f"{out}: {len(rows)} sequences -> {len({r['cluster_id'] for r in rows})} clusters"
          + (f", {len({r['component_id'] for r in rows})} components" if has_comps else "")
          + f" ({n_cent} centroid, {len(rows) - n_cent} not)", file=sys.stderr)
    print(f"built against {b_path.name}"
          + (f", engine {engine}" if engine else "")
          + (f", identity {identity}" if identity else ""), file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(0)
