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
    p.add_argument("--s5"); p.add_argument("--s6")
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
    s5 = a.s5 or os.path.join(run_dir, "05-union-with-components.tsv")
    s6 = a.s6 or os.path.join(run_dir, "06-nr.tsv")

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
    annot = read_tsv(s5)
    nr = read_tsv(s6)

    L.append("## Deliverables")
    L.append("")
    L.append("| step | file | rows |")
    L.append("| --- | --- | --- |")
    rows_by_step = [
        (1, s1, len(union) if union is not None else None),
        (2, s2, len(clusters) if clusters is not None else None),
        (3, s3, len(aligned) if aligned is not None else None),
        (4, s4, n_fasta),
    ]
    # Steps 1-4 are the driver's default range, so an absent one is a partial run and
    # worth naming. Step 6 is opt-in (`--only 6`), so a run that never asked for it is
    # not partial; listing it as "not built" would report an absence nobody chose.
    if annot is not None:
        rows_by_step.append((5, s5, len(annot)))
    if nr is not None:
        rows_by_step.append((6, s6, len(nr)))
    for n, path, rows in rows_by_step:
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
        prov = read_json(os.path.join(run_dir, "02-clusters.intermediates", "provenance.json"))
        if prov:
            cov = prov.get("member_cov")
            L.append(f"- engine: {prov.get('engine_version', prov.get('engine', '?'))}"
                     f" — centroids ordered by {prov.get('greedy_ordering', '?')}, member "
                     f"coverage {cov if cov else 'not enforced (identity only)'}")
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
            if stats.get("aligner_command"):
                L.append(f"- aligner: {stats['aligner_command']}")
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

    # -- step 5 ---------------------------------------------------------------
    # The cluster-level result put back at union-row level. Same partition as step 2,
    # so the interesting number is the non-centroid count, which steps 2-4 never show.
    if annot is not None:
        cent = sum(1 for r in annot if r.get("is_centroid") == "yes")
        comps = {r.get("component_id") for r in annot if r.get("component_id")}
        L.append("## 5. Union with components")
        L.append("")
        L.append(f"- {len(annot)} rows, one per union sequence — "
                 f"{cent} centroids, {len(annot) - cent} non-centroid members")
        L.append(f"- over {len({r.get('cluster_id') for r in annot})} clusters "
                 f"and {len(comps)} components")
        if union is not None and len(annot) != len(union):
            L.append(f"- NOTE: {len(annot)} rows against {len(union)} union rows; "
                     f"step 5 is meant to be row-aligned with step 1")
        if clusters is not None and cent != len(clusters):
            L.append(f"- NOTE: {cent} centroids against {len(clusters)} clusters in step 2")
        L.append("")

    # -- step 6 ---------------------------------------------------------------
    # Branch C. Parallel to steps 2-4, not downstream of them, so this section says
    # nothing about the engine; the crosswalk is where the branch-B provenance lives.
    if nr is not None:
        L.append("## 6. Non-redundant set")
        L.append("")
        prov = read_json(os.path.join(run_dir, "06-nr.intermediates", "provenance.json"))
        n_union = (prov or {}).get("n_union_rows") or (len(union) if union else None)
        L.append(f"- {len(nr)} distinct sequences"
                 + (f" from {n_union} union rows" if n_union else "")
                 + " — 100% identity by md5 of the normalized sequence")
        if prov:
            L.append(f"- {prov.get('n_duplicate_groups', '?')} duplicate groups covering "
                     f"{prov.get('n_rows_in_duplicate_groups', '?')} union rows, "
                     f"largest {prov.get('largest_group', '?')}")
            L.append(f"- built from `{prov.get('input', '?')}` "
                     f"(md5 `{str(prov.get('input_md5', ''))[:12]}`), "
                     f"accession versions {'kept' if prov.get('keeps_accession_versions') else 'stripped'}")
        # The universe must be conserved: every union row lands in exactly one group.
        members = sum(int(r["n_members"]) for r in nr if (r.get("n_members") or "").isdigit())
        if union is not None and members != len(union):
            L.append(f"- NOTE: members sum to {members}, not the {len(union)} union rows")

        # 06-nr.fasta has two byte-states under one name: build_nr.py writes field 5
        # of the header empty, crosswalk.py fills it. Report which one is on disk so
        # the two are distinguishable without re-reading the file.
        fa = os.path.join(run_dir, "06-nr.fasta")
        if os.path.exists(fa):
            with open(fa) as f:
                heads = [l[1:].rstrip("\n") for l in f if l.startswith(">")]
            filled = sum(1 for h in heads if h.count("|") == 4 and h.rsplit("|", 1)[1])
            L.append(f"- `06-nr.fasta`: {len(heads)} records, "
                     f"`>identifier|accession|alt_accessions|pazy_id|component`; "
                     f"component filled on {filled}"
                     + ("" if filled else " (crosswalk not run)"))

        xw = read_tsv(os.path.join(run_dir, "06-nr-to-clusters.tsv"))
        if xw is None:
            L.append("- crosswalk: not built (needs step 2/3 in this run directory)")
        else:
            # engine is still a column; identity moved to the crosswalk sidecar,
            # since it is a per-run constant and was 0.9 on every run to date.
            xwp = read_json(os.path.join(run_dir, "06-nr.intermediates", "crosswalk.json")) or {}
            eng = next((r.get("engine") for r in xw if r.get("engine")), "") or xwp.get("engine", "")
            ident = xwp.get("identity", "")
            cent = sum(1 for r in xw if r.get("is_centroid") == "yes")
            L.append(f"- crosswalk: {len(xw)} sequences over "
                     f"{len({r.get('cluster_id') for r in xw})} clusters / "
                     f"{len({r.get('component_id') for r in xw if r.get('component_id')})} components; "
                     f"{cent} are centroids of branch B, {len(xw) - cent} are not")
            if eng or ident or xwp.get("built_against"):
                L.append(f"  - built against {xwp.get('built_against', 'branch B')}, "
                         f"{eng or 'unrecorded engine'}"
                         + (f" at {ident} identity" if ident else ""))
            # crosswalk.py writes the file only after the containment check passes, so
            # the file existing IS the result. This is not a check re-run here.
            L.append("  - md5 containment held at build time: no identical sequence "
                     "split across two clusters or components")
        L.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(L).rstrip() + "\n")
    print(f"wrote   : {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
