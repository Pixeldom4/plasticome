#!/usr/bin/env python3
"""Compare this run's partition against the 2026-07-23 singleton-cleaned-union passes.

Keyed on **sequence_md5**, not node_id: `m####` ids are assigned in md5-sorted order
within a run, so they are meaningless across runs. Every metric below is computed on
the node set the two runs actually share.

Reports, for each (this run vs pass1/usearch, pass3/usearch, pass3/diamond, ...):
  - node-set overlap: shared / this-only / other-only sequences
  - edge-set Jaccard on the passing edges, remapped to md5 pairs
  - partition agreement: Rand index + Adjusted Rand Index (same estimator as the
    2026-07-23 compare.py, so the numbers are directly comparable to that REPORT)
  - per-node disagreements after greedy best-overlap component matching

Writes comparison_vs_2026-07-23.json + node_disagreements_vs_2026-07-23.tsv here.
"""
import csv
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE.parents[1]
OLD = RUNS / "2026-07-23_singleton-cleaned-union.v1.v260701"

THIS = {"assign": HERE / "component_assignment_clusters_2026-07-30.csv",
        "edges": HERE / "component_edges_clusters_2026-07-30.csv",
        "nodes": HERE / "combined_nodes.tsv"}

OTHERS = {
    "pass1/usearch": {
        "assign": OLD / "pass1-alignment/usearch/results/component_assignment_singleton_2026-07-23.csv",
        "edges": OLD / "pass1-alignment/usearch/results/component_edges_singleton_2026-07-23.csv",
        "nodes": OLD / "pass1-alignment/usearch/results/combined_nodes.tsv"},
    "pass1/diamond": {
        "assign": OLD / "pass1-alignment/diamond/results/component_assignment.csv",
        "edges": OLD / "pass1-alignment/diamond/results/component_edges.csv",
        "nodes": OLD / "pass1-alignment/usearch/results/combined_nodes.tsv"},
    "pass3/usearch": {
        "assign": OLD / "pass3-alignment/usearch/results/component_assignment_singleton_2026-07-23.csv",
        "edges": OLD / "pass3-alignment/usearch/results/component_edges_singleton_2026-07-23.csv",
        "nodes": OLD / "pass3-alignment/usearch/results/combined_nodes.tsv"},
    "pass3/diamond": {
        "assign": OLD / "pass3-alignment/diamond/results/component_assignment.csv",
        "edges": OLD / "pass3-alignment/diamond/results/component_edges.csv",
        "nodes": OLD / "pass3-alignment/usearch/results/combined_nodes.tsv"},
}


def node_md5(nodes_tsv):
    """node_id -> sequence_md5, from the run's own node roster."""
    with Path(nodes_tsv).open(newline="") as fh:
        return {r["node_id"]: r["sequence_md5"] for r in csv.DictReader(fh, delimiter="\t")}


def load_assign(spec):
    """sequence_md5 -> component label (labels are opaque; only the grouping matters)."""
    by_node = node_md5(spec["nodes"])
    out = {}
    with Path(spec["assign"]).open(newline="") as fh:
        rdr = csv.DictReader(fh)
        col = "component" if "component" in rdr.fieldnames else "component_id"
        for r in rdr:
            md5 = r.get("sequence_md5") or by_node.get(r["node_id"], "")
            if md5:
                out[md5] = r[col]
    return out, by_node


def load_edges(spec, by_node):
    """set of frozenset({md5_a, md5_b}) over the passing edges."""
    s = set()
    with Path(spec["edges"]).open(newline="") as fh:
        for r in csv.DictReader(fh):
            a, b = by_node.get(r["node_a"]), by_node.get(r["node_b"])
            if a and b:
                s.add(frozenset((a, b)))
    return s


def clusters(assign):
    c = defaultdict(set)
    for n, k in assign.items():
        c[k].add(n)
    return c


def rand_indices(a, b):
    """Rand index and Adjusted Rand Index over the shared node set."""
    nodes = sorted(set(a) & set(b))
    n = len(nodes)
    ca, cb = clusters({k: a[k] for k in nodes}), clusters({k: b[k] for k in nodes})
    idx_b = {x: k for k, s in cb.items() for x in s}
    cont = defaultdict(int)
    for k, s in ca.items():
        for x in s:
            cont[(k, idx_b[x])] += 1
    comb2 = lambda x: x * (x - 1) // 2
    sum_ij = sum(comb2(v) for v in cont.values())
    sum_ai = sum(comb2(len(s)) for s in ca.values())
    sum_bj = sum(comb2(len(s)) for s in cb.values())
    total = comb2(n)
    ri = (total + 2 * sum_ij - sum_ai - sum_bj) / total
    expected = sum_ai * sum_bj / total
    max_index = (sum_ai + sum_bj) / 2
    ari = (sum_ij - expected) / (max_index - expected) if max_index != expected else 1.0
    return n, ri, ari


def match_components(a, b):
    """Greedy best-overlap matching a-comp -> b-comp; per-node disagreements on shared nodes."""
    nodes = sorted(set(a) & set(b))
    ca = clusters({k: a[k] for k in nodes})
    cb = clusters({k: b[k] for k in nodes})
    match = {}
    for k, s in ca.items():
        overlap = defaultdict(int)
        for x in s:
            overlap[b[x]] += 1
        match[k] = max(overlap.items(), key=lambda kv: kv[1])[0]
    return [(x, a[x], len(ca[a[x]]), b[x], len(cb[b[x]])) for x in nodes if match[a[x]] != b[x]]


def summarize(assign):
    c = clusters(assign)
    sizes = sorted((len(s) for s in c.values()), reverse=True)
    return {"n_nodes": len(assign), "n_components": len(c),
            "singletons": sum(1 for x in sizes if x == 1),
            "largest": sizes[0] if sizes else 0, "top10_sizes": sizes[:10]}


def main():
    this_a, this_by_node = load_assign(THIS)
    this_e = load_edges(THIS, this_by_node)
    # label lookup for the disagreement table
    with Path(THIS["nodes"]).open(newline="") as fh:
        label = {r["sequence_md5"]: (r.get("identifier") or r.get("accession") or r["node_id"])
                 for r in csv.DictReader(fh, delimiter="\t")}

    out = {"this_run": {"path": str(HERE), **summarize(this_a), "edges": len(this_e)},
           "comparisons": {}}
    disagreements = []

    for name, spec in OTHERS.items():
        other_a, other_by_node = load_assign(spec)
        other_e = load_edges(spec, other_by_node)
        shared = set(this_a) & set(other_a)
        se, oe = {e for e in this_e if e <= shared}, {e for e in other_e if e <= shared}
        n, ri, ari = rand_indices(this_a, other_a)
        dis = match_components(this_a, other_a)
        for md5, ca, sa, cb, sb in dis:
            disagreements.append([name, label.get(md5, ""), md5, ca, sa, cb, sb])
        out["comparisons"][name] = {
            "other": {**summarize(other_a), "edges": len(other_e)},
            "node_overlap": {"shared": len(shared),
                             "this_only": len(set(this_a) - set(other_a)),
                             "other_only": len(set(other_a) - set(this_a))},
            "edge_overlap_on_shared_nodes": {
                "shared": len(se & oe), "this_only": len(se - oe), "other_only": len(oe - se),
                "jaccard": round(len(se & oe) / len(se | oe), 5) if (se | oe) else 1.0},
            "partition_agreement": {"n_nodes_compared": n,
                                    "rand_index": round(ri, 6),
                                    "adjusted_rand_index": round(ari, 6),
                                    "nodes_in_differently_matched_component": len(dis)},
        }

    (HERE / "comparison_vs_2026-07-23.json").write_text(json.dumps(out, indent=2))
    with (HERE / "node_disagreements_vs_2026-07-23.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["vs", "label", "sequence_md5", "this_comp", "this_comp_size",
                    "other_comp", "other_comp_size"])
        w.writerows(sorted(disagreements, key=lambda r: (r[0], r[3], r[1])))

    print(json.dumps(out, indent=2))
    print(f"\n{len(disagreements)} node/comparison disagreement(s) "
          f"-> node_disagreements_vs_2026-07-23.tsv")


if __name__ == "__main__":
    main()
