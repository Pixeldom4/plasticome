#!/usr/bin/env python3
"""Steps 2-3 - recompute the graph over the combined node set and classify the delta.

Edge criterion (applied in post-processing, never as a search parameter):
    >=30% amino-acid identity AND e-value < 1e-5
Partition: undirected single-linkage connected components.

Component IDs are merge-aware and canonical (D4):
    inherited  -> the v1 component number
    split      -> v1 number + a/b/... in canonical order
    merged     -> "M<smallest parent v1 number>", merged_from lists all parents
    novel      -> "N<k>", k by canonical order (smallest member accession)

The v1-only baseline partition is the induced subgraph on v1 nodes *of this same
edge set*, so merges attributable to delta nodes are separated from merges that
Step 0 already saw at the 30%/1e-5 boundary.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

ID_MIN = 30.0
EVALUE_MAX = 1e-5

# usearch11 reports E = K * m * N * 2^-bits, with N = total db letters (verified:
# median log10(e_492/e_213) = 0.3537 vs predicted log10(N2/N0) = 0.3538, identical
# alignments). The e-value therefore tightens as the node set grows. To keep the
# edge criterion a property of the *pair* rather than of the database size, rescale
# every e-value back to the paper's 213-node search space before thresholding.
N_REF = 74254  # db letters of step0/v1_nodes_213.fasta (the paper's seed set)


class DSU:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def components(node_ids, edges):
    d = DSU()
    for n in node_ids:
        d.find(n)
    for q, t, *_ in edges:
        d.union(q, t)
    groups = defaultdict(set)
    for n in node_ids:
        groups[d.find(n)].add(n)
    return groups


def db_letters(fasta):
    return sum(len(l.strip()) for l in Path(fasta).open() if not l.startswith(">"))


def load_pairs(paths, scale=1.0):
    """Deduplicate to the best-scoring HSP per pair, then apply the dual threshold.

    Two usearch quirks are corrected here.

    (1) Multiple HSPs per pair. Selecting by max *identity* picks short, high-id,
    statistically insignificant fragments and silently drops strong edges - this is
    what step0 did: it kept n157<->n160's 69.6%-id / e=2.7e-5 / 40-bit HSP over the
    63.1%-id / e=4e-78 / 282-bit one, then failed it on e-value. The best alignment
    for a pair is the highest-scoring one; ties break on lower e-value.

    (2) Orientation dependence. `-allpairs_local` always makes the earlier FASTA
    sequence the query, and the HSP it returns depends on that role (m0275/m0438:
    29.8% id / 195 bits one way, 35.8% id / 179 bits the other). Pass both the
    forward and reverse-order pair files so every pair is scored in both
    orientations; the graph then depends on the node set, not on file order.

    `scale` multiplies observed e-values to renormalize them onto a fixed search
    space (N_REF / N_run).
    """
    if isinstance(paths, (str, Path)):
        paths = [paths]
    best = {}
    total = 0
    for path in paths:
        for line in Path(path).open():
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            q, t, pid, ev, bits = f[0], f[1], float(f[2]), float(f[3]), float(f[4])
            total += 1
            if q == t:
                continue
            k = (q, t) if q < t else (t, q)
            cur = best.get(k)
            if cur is None or bits > cur[2] or (bits == cur[2] and ev < cur[1]):
                best[k] = (pid, ev * scale, bits)
    return ({k: v for k, v in best.items()
             if v[0] >= ID_MIN and v[1] < EVALUE_MAX}, total, len(best))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--date", default="2026-07-09")
    ap.add_argument("--tag", default="v2")
    ap.add_argument("--evalue-mode", choices=["paper_fixed", "raw_db"],
                    default="paper_fixed",
                    help="paper_fixed: renormalize e-values to the 213-node search "
                         "space (criterion independent of node-set size). "
                         "raw_db: threshold the e-value as reported for this db.")
    args = ap.parse_args()

    nodes = {r["node_id"]: r for r in
             csv.DictReader(Path(f"{args.prefix}_nodes.tsv").open(), delimiter="\t")}
    n_run = db_letters(f"{args.prefix}.fasta")
    scale = (N_REF / n_run) if args.evalue_mode == "paper_fixed" else 1.0
    pair_files = [f"{args.prefix}_pairs.tsv", f"{args.prefix}_rev_pairs.tsv"]
    edges, n_raw, n_pairs = load_pairs(pair_files, scale=scale)
    edge_list = [(q, t, *v) for (q, t), v in sorted(edges.items())]

    v1_nodes = {n for n, r in nodes.items() if r["source"] == "v1"}
    delta_nodes = {n for n, r in nodes.items() if r["source"] == "delta"}

    # --- baseline: v1-only induced subgraph of THIS edge set
    v1_edges = [e for e in edge_list if e[0] in v1_nodes and e[1] in v1_nodes]
    base_groups = components(v1_nodes, v1_edges)
    base_of = {n: root for root, mem in base_groups.items() for n in mem}

    # --- full partition over all nodes
    full_groups = components(set(nodes), edge_list)

    # --- canonical component labelling (D4)
    def acc_key(n):
        return (nodes[n]["accession"], n)

    def canon(mem):
        return min(mem, key=acc_key)

    comps = []
    for root, mem in full_groups.items():
        v1_mem = mem & v1_nodes
        v1_labels = sorted({nodes[n]["v1_component"] for n in v1_mem if nodes[n]["v1_component"]},
                           key=lambda s: int(s))
        base_ids = {base_of[n] for n in v1_mem}
        comps.append({"members": mem, "v1_labels": v1_labels, "n_base": len(base_ids),
                      "canon": canon(mem)})

    # detect v1 labels split across >1 final component
    label_spread = defaultdict(list)
    for c in comps:
        for lab in c["v1_labels"]:
            label_spread[lab].append(c)
    for lab, cs in label_spread.items():
        if len(cs) > 1:
            for i, c in enumerate(sorted(cs, key=lambda c: acc_key(c["canon"]))):
                c.setdefault("split_suffix", {})[lab] = chr(ord("a") + i)

    novel_sorted = sorted([c for c in comps if not c["v1_labels"]],
                          key=lambda c: acc_key(c["canon"]))
    for k, c in enumerate(novel_sorted, start=1):
        c["component_id"] = f"N{k}"
        c["kind"] = "novel"
    for c in comps:
        if c["v1_labels"]:
            if len(c["v1_labels"]) == 1:
                lab = c["v1_labels"][0]
                sfx = c.get("split_suffix", {}).get(lab, "")
                c["component_id"] = f"{lab}{sfx}"
                c["kind"] = "split" if sfx else "inherited"
            else:
                c["component_id"] = "M" + c["v1_labels"][0]
                c["kind"] = "merged"

    comps.sort(key=lambda c: acc_key(c["canon"]))

    # --- Step 3: per-delta outcome, relative to the v1 baseline partition
    outcome_of, comp_of = {}, {}
    for c in comps:
        for n in c["members"]:
            comp_of[n] = c
        if not (c["members"] & v1_nodes):
            oc = "novel"
        elif c["n_base"] >= 2:
            oc = "merge"
        else:
            oc = "join"
        for n in c["members"] & delta_nodes:
            outcome_of[n] = oc

    # merges that already existed in the v1-only baseline (not caused by delta nodes)
    baseline_merges = []
    for root, mem in base_groups.items():
        labs = sorted({nodes[n]["v1_component"] for n in mem if nodes[n]["v1_component"]},
                      key=lambda s: int(s))
        if len(labs) > 1:
            baseline_merges.append(labs)

    # ---- outputs
    stem = f"{args.tag}_{args.date}"
    with open(f"component_edges_{stem}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["node_a", "node_b", "pct_id", "evalue", "bits",
                    "acc_a", "acc_b", "a_is_new", "b_is_new"])
        for q, t, pid, ev, bits in edge_list:
            w.writerow([q, t, f"{pid:.1f}", f"{ev:.3e}", f"{bits:.1f}",
                        nodes[q]["accession"], nodes[t]["accession"],
                        nodes[q]["is_new"], nodes[t]["is_new"]])

    with open(f"component_assignment_{stem}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["node_id", "sequence_md5", "accession", "sequence_length",
                    "component_id", "component_kind", "is_new", "outcome",
                    "v1_component", "merged_from", "cath", "pazy_id", "activity"])
        for n in sorted(nodes, key=lambda n: (comps.index(comp_of[n]), acc_key(n))):
            r, c = nodes[n], comp_of[n]
            w.writerow([n, r["sequence_md5"], r["accession"], r["sequence_length"],
                        c["component_id"], c["kind"], r["is_new"],
                        outcome_of.get(n, ""), r["v1_component"],
                        ";".join(c["v1_labels"]) if c["kind"] == "merged" else "",
                        r["cath"], r["pazy_id"], r["activity"]])

    with open(f"delta_component_summary_{stem}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["node_id", "accession", "pazy_id", "sequence_length",
                    "outcome", "component_id", "merged_from", "activity",
                    "n_edges_to_v1", "n_edges_to_delta"])
        for n in sorted(delta_nodes, key=acc_key):
            c = comp_of[n]
            e_v1 = sum(1 for q, t, *_ in edge_list
                       if (q == n and t in v1_nodes) or (t == n and q in v1_nodes))
            e_d = sum(1 for q, t, *_ in edge_list
                      if (q == n and t in delta_nodes) or (t == n and q in delta_nodes))
            w.writerow([n, nodes[n]["accession"], nodes[n]["pazy_id"],
                        nodes[n]["sequence_length"], outcome_of[n], c["component_id"],
                        ";".join(c["v1_labels"]) if c["kind"] == "merged" else "",
                        nodes[n]["activity"], e_v1, e_d])

    with open(f"component_members_{stem}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["component_id", "kind", "n_members", "n_v1", "n_new",
                    "v1_labels", "merged_from_baseline_groups", "cath_folds"])
        for c in comps:
            cf = sorted({nodes[n]["cath"] for n in c["members"] if nodes[n]["cath"]})
            w.writerow([c["component_id"], c["kind"], len(c["members"]),
                        len(c["members"] & v1_nodes), len(c["members"] & delta_nodes),
                        ";".join(c["v1_labels"]), c["n_base"], ";".join(cf)])

    oc_counts = defaultdict(int)
    for n in delta_nodes:
        oc_counts[outcome_of[n]] += 1
    stats = {
        "prefix": args.prefix, "evalue_mode": args.evalue_mode,
        "db_letters": n_run, "evalue_scale": round(scale, 6),
        "n_nodes": len(nodes), "n_v1": len(v1_nodes),
        "n_delta": len(delta_nodes), "raw_alignments": n_raw,
        "aligned_pairs": n_pairs, "edges_passing": len(edge_list),
        "v1_induced_edges": len(v1_edges),
        "v1_induced_components": len(base_groups),
        "baseline_merges": baseline_merges,
        "n_components": len(comps),
        "kinds": dict(sorted(defaultdict(int, {
            k: sum(1 for c in comps if c["kind"] == k)
            for k in {c["kind"] for c in comps}}).items())),
        "delta_outcomes": dict(oc_counts),
        "singletons": sum(1 for c in comps if len(c["members"]) == 1),
        "largest_component": max(len(c["members"]) for c in comps),
    }
    Path(f"stats_{stem}.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
