#!/usr/bin/env python3
"""Compare the usearch vs diamond partitions on the shared 412-node set.

Both tools ran the identical method (allpairs local, >=30% aaid AND e<1e-5 post-filter,
single-linkage) on the identical md5-unique node set (m#### ids). This quantifies how
much the choice of aligner moves the final component assignment:
  - edge-set overlap (Jaccard) on the passing edge lists
  - partition agreement: pair-counting Rand index + Adjusted Rand Index
  - per-node disagreements after greedy best-overlap component matching
Writes comparison.json and node_disagreements.tsv into this folder.
"""
import csv, json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
U_ASSIGN = HERE / "usearch/results/component_assignment_singleton_2026-07-23.csv"
D_ASSIGN = HERE / "diamond/results/component_assignment.csv"
U_EDGES  = HERE / "usearch/results/component_edges_singleton_2026-07-23.csv"
D_EDGES  = HERE / "diamond/results/component_edges.csv"


def load_assign(path):
    """node_id -> component label. Handles usearch (component_id='C001') and
    diamond (component='1') column shapes; labels are opaque, only grouping matters."""
    m = {}
    with path.open() as f:
        r = csv.DictReader(f)
        col = "component" if "component" in r.fieldnames else "component_id"
        for row in r:
            m[row["node_id"]] = row[col]
    return m


def load_edges(path):
    """set of frozenset({a,b}) passing edges."""
    s = set()
    with path.open() as f:
        r = csv.DictReader(f)
        for row in r:
            s.add(frozenset((row["node_a"], row["node_b"])))
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
    # contingency
    idx_b = {}
    for k, s in cb.items():
        for x in s:
            idx_b[x] = k
    cont = defaultdict(int)
    for k, s in ca.items():
        for x in s:
            cont[(k, idx_b[x])] += 1
    comb2 = lambda x: x * (x - 1) // 2
    sum_ij = sum(comb2(v) for v in cont.values())
    sum_ai = sum(comb2(len(s)) for s in ca.values())
    sum_bj = sum(comb2(len(s)) for s in cb.values())
    total = comb2(n)
    # Rand index
    agree = total + 2 * sum_ij - sum_ai - sum_bj
    ri = agree / total
    # Adjusted Rand
    expected = sum_ai * sum_bj / total
    max_index = (sum_ai + sum_bj) / 2
    ari = (sum_ij - expected) / (max_index - expected) if max_index != expected else 1.0
    return n, ri, ari


def match_components(a, b):
    """Greedy best-overlap matching of usearch comps -> diamond comps; return per-node disagreements."""
    nodes = sorted(set(a) & set(b))
    ca, cb = clusters(a), clusters(b)
    # For each usearch comp, find diamond comp with max overlap = its "match"
    match = {}
    for k, s in ca.items():
        best, bestn = None, -1
        overlap = defaultdict(int)
        for x in s:
            overlap[b[x]] += 1
        for dk, cnt in overlap.items():
            if cnt > bestn:
                best, bestn = dk, cnt
        match[k] = best
    disagree = []
    for x in nodes:
        if match[a[x]] != b[x]:
            disagree.append((x, a[x], len(ca[a[x]]), b[x], len(cb[b[x]])))
    return match, disagree


def main():
    ua, da = load_assign(U_ASSIGN), load_assign(D_ASSIGN)
    ue, de = load_edges(U_EDGES), load_edges(D_EDGES)

    inter = len(ue & de)
    union = len(ue | de)
    jac = inter / union if union else 1.0

    n, ri, ari = rand_indices(ua, da)
    match, disagree = match_components(ua, da)

    out = {
        "n_nodes_shared": n,
        "usearch": {"n_components": len(clusters(ua)), "edges": len(ue),
                    "singletons": sum(1 for s in clusters(ua).values() if len(s) == 1),
                    "largest": max(len(s) for s in clusters(ua).values())},
        "diamond": {"n_components": len(clusters(da)), "edges": len(de),
                    "singletons": sum(1 for s in clusters(da).values() if len(s) == 1),
                    "largest": max(len(s) for s in clusters(da).values())},
        "edge_overlap": {"shared": inter, "usearch_only": len(ue - de),
                         "diamond_only": len(de - ue), "union": union, "jaccard": round(jac, 5)},
        "partition_agreement": {"rand_index": round(ri, 6),
                                "adjusted_rand_index": round(ari, 6),
                                "nodes_in_differently_matched_component": len(disagree)},
    }
    (HERE / "comparison.json").write_text(json.dumps(out, indent=2))

    with (HERE / "node_disagreements.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["node_id", "usearch_comp", "usearch_comp_size", "diamond_comp", "diamond_comp_size"])
        for row in sorted(disagree, key=lambda r: (r[1], r[0])):
            w.writerow(row)

    print(json.dumps(out, indent=2))
    print(f"\n{len(disagree)} node(s) fall in a non-corresponding component -> node_disagreements.tsv")


if __name__ == "__main__":
    main()
