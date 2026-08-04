#!/usr/bin/env python3
"""Write a concise summary.md for a pipeline run directory.

Reads only what is already on disk -- the four deliverables plus the sidecar
JSON that steps 1 and 3 leave behind -- so it is safe to re-run at any time and
correct for a partial run: a step whose deliverable is absent is reported as
"not built" rather than guessed at.

Usage
-----
  python3 summarize_run.py --run-dir runs/<run> [--s1 F] [--s2 F] [--s3 F] [--s4 F]
                           [--v1 F] [--v2 F] [--v2-seq MODE] [--seeds F]
                           [--id PCT] [--date YYYY-MM-DD] [-o summary.md]

Deliverable paths default to the canonical 01-union.tsv / 02-clusters.tsv /
03-alignment.tsv names under --run-dir; the driver passes them explicitly
because a run made before those names settled may use others.
"""
import argparse
import csv
import glob
import json
import os
import sys
from collections import Counter

csv.field_size_limit(10 ** 7)  # rep_aa_sequence makes for long fields


def read_tsv(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None


def median(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else 0


def fasta_count(path):
    if not path or not os.path.exists(path):
        return None
    with open(path) as f:
        return sum(1 for line in f if line.startswith(">"))


def rel(path, base):
    """Path relative to the repo root when it sits under it, else as given."""
    try:
        r = os.path.relpath(path, base)
    except ValueError:
        return path
    return path if r.startswith(os.pardir) else r


def counts_line(counter, keys):
    """'v260701 398 / both 75 / v1.1 137', keys first, then any extras seen."""
    order = list(keys) + [k for k in counter if k not in keys]
    return " / ".join(f"{k} {counter[k]}" for k in order if counter[k])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--run-dir", required=True)
    p.add_argument("--s1"); p.add_argument("--s2")
    p.add_argument("--s3"); p.add_argument("--s4")
    p.add_argument("--v1"); p.add_argument("--v2"); p.add_argument("--seeds")
    p.add_argument("--v2-seq", default="")
    p.add_argument("--id", default="")
    p.add_argument("--date", default="", help="run date for the header (the driver passes today's)")
    p.add_argument("-o", "--out", default="")
    return p.parse_args()


def main():
    a = parse_args()
    run_dir = os.path.abspath(a.run_dir)
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    out_path = a.out or os.path.join(run_dir, "summary.md")

    s1 = a.s1 or os.path.join(run_dir, "01-union.tsv")
    s2 = a.s2 or os.path.join(run_dir, "02-clusters.tsv")
    s3 = a.s3 or os.path.join(run_dir, "03-alignment.tsv")
    s4 = a.s4

    L = []  # output lines
    L.append(f"# {os.path.basename(run_dir)}")
    L.append("")
    if a.date:
        L.append(f"Summary written {a.date} by `lib/summary/summarize_run.py`.")
        L.append("")

    # -- inputs ---------------------------------------------------------------
    L.append("## Inputs")
    L.append("")
    for label, path in (("v1.1", a.v1), ("v260701", a.v2), ("seeds", a.seeds)):
        if path:
            L.append(f"- {label}: `{rel(path, repo)}`")
    if a.v2_seq:
        L.append(f"- v260701 sequence rule: `{a.v2_seq}`")
    if a.id:
        L.append(f"- clustering identity: {a.id}")
    L.append("")

    # -- deliverables ---------------------------------------------------------
    union = read_tsv(s1)
    clusters = read_tsv(s2)
    aligned = read_tsv(s3)
    n_fasta = fasta_count(s4)

    L.append("## Deliverables")
    L.append("")
    L.append("| step | file | rows |")
    L.append("| --- | --- | --- |")
    for n, path, rows in (
        (1, s1, len(union) if union is not None else None),
        (2, s2, len(clusters) if clusters is not None else None),
        (3, s3, len(aligned) if aligned is not None else None),
        (4, s4, n_fasta),
    ):
        name = f"`{os.path.basename(path)}`" if path else "—"
        L.append(f"| {n} | {name} | {rows if rows is not None else 'not built'} |")
    L.append("")

    # -- step 1 ---------------------------------------------------------------
    if union is not None:
        src = Counter(r.get("source", "") for r in union)
        with_seq = sum(1 for r in union if (r.get("aa_sequence") or "").strip())
        blank = len(union) - with_seq
        L.append("## 1. Union")
        L.append("")
        L.append(f"- {len(union)} rows — {counts_line(src, ['v260701', 'both', 'v1.1'])}")
        L.append(f"- sequences: {with_seq} present, {blank} blank")
        fetched = read_json(os.path.join(run_dir, "fetched_sequences.json"))
        if fetched:
            unresolved = fetched.get("unresolved") or []
            L.append(f"- accession lookups: {len(fetched.get('resolved') or {})} resolved, "
                     f"{len(unresolved)} unresolved")
            by_stage = fetched.get("by_stage") or {}
            if by_stage:
                L.append("  - resolved by: " + ", ".join(f"{k} {v}" for k, v in by_stage.items()))
            if unresolved:
                shown = ", ".join(f"`{x}`" for x in unresolved[:12])
                more = f" (+{len(unresolved) - 12} more)" if len(unresolved) > 12 else ""
                L.append(f"  - unresolved: {shown}{more}")
        L.append("")

    # -- step 2 ---------------------------------------------------------------
    if clusters is not None:
        sizes = [int(r["size"]) for r in clusters if (r.get("size") or "").isdigit()]
        lens = [int(r["rep_seq_len"]) for r in clusters if (r.get("rep_seq_len") or "").isdigit()]
        rep_src = Counter(r.get("rep_source", "") for r in clusters)
        L.append(f"## 2. Clusters{f' ({a.id} identity)' if a.id else ''}")
        L.append("")
        L.append(f"- {len(clusters)} clusters over {sum(sizes)} sequences")
        L.append(f"- {sum(1 for s in sizes if s == 1)} singleton clusters, "
                 f"largest {max(sizes) if sizes else 0}")
        L.append(f"- representatives — {counts_line(rep_src, ['v260701', 'both', 'v1.1', 'seed'])}")
        if lens:
            L.append(f"- representative length: min {min(lens)} / median {median(lens)} / max {max(lens)} aa")
        L.append("")

    # -- step 3 ---------------------------------------------------------------
    if aligned is not None:
        comp = Counter(r.get("component_id", "") for r in aligned)
        stats = None
        for p in sorted(glob.glob(os.path.join(run_dir, "03-alignment.intermediates", "stats_*.json"))):
            stats = read_json(p) or stats
        L.append("## 3. Components")
        L.append("")
        L.append(f"- {len(comp)} components over {len(aligned)} clusters")
        top = ", ".join(f"{cid} ({n})" for cid, n in comp.most_common(3))
        L.append(f"- {sum(1 for n in comp.values() if n == 1)} single-cluster components; largest: {top}")
        if stats:
            if stats.get("edge_criterion"):
                L.append(f"- edge criterion: {stats['edge_criterion']}")
            if stats.get("edges_passing") is not None:
                L.append(f"- {stats.get('aligned_pairs', '?')} aligned pairs, "
                         f"{stats['edges_passing']} edges passing")
            if stats.get("n_nodes_with_v1_overlay") is not None:
                L.append(f"- {stats['n_nodes_with_v1_overlay']} clusters carry a v1.1 overlay; "
                         f"{stats.get('n_v1_labels_split_across_components', 0)} v1 labels split across components")
        L.append("")

    # -- step 4 ---------------------------------------------------------------
    if n_fasta is not None:
        L.append("## 4. FASTA")
        L.append("")
        L.append(f"- {n_fasta} centroid records in `{os.path.basename(s4)}`")
        if clusters is not None and n_fasta != len(clusters):
            L.append(f"- NOTE: record count differs from the {len(clusters)} clusters in step 2")
        L.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(L).rstrip() + "\n")
    print(f"wrote   : {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
