#!/usr/bin/env python3
"""Steps 2-3 (v6) - flat similarity graph + connected-component partition.

Edge criterion (post-processing only, never a search parameter):
    >=30% amino-acid identity AND e-value < 1e-5
Partition: undirected single-linkage connected components over ALL nodes.

FLAT de-novo run (decisions 1-4), identical method to v4/v5:
  * No v1 anchor, no delta, no join/merge/novel typing.
  * E-values thresholded AS REPORTED by usearch at this run's search space
    (no down-scaling to the paper's 213-seq space; decision 2). `--scale 1.0` = identity.
  * HSP selection: best HSP per unordered pair by MAX BITS (ties -> lower e-value),
    over both FASTA orientations (adopted v3 s7 convention; accepted B4 cost).
  * Canonical component labels/IDs resolve on the PL identifier (v6: every node has a
    stable, unique PL<n>, so no accession/enzyme-name fallback is ever needed).

v1 component/cath columns are an OVERLAY only (Step 3 interpretation), never used to
seed or gate the partition.
"""
import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

ID_MIN = 30.0
EVALUE_MAX = 1e-5


def pl_num(ident):
    m = re.search(r"(\d+)", ident or "")
    return int(m.group(1)) if m else 10**9


def node_label(node):
    """Canonical label: PL identifier, else accession, else enzyme name, else node id."""
    return (node.get("identifier") or node.get("accession")
            or node.get("enzyme_name") or node["node_id"])


def canon_key(nodes, n):
    """Canonical ordering: by PL number (stable, unique), node_id as final tiebreak."""
    node = nodes[n]
    return (pl_num(node.get("identifier", "")), node_label(node), n)


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


def load_pairs(paths, scale=1.0):
    """Best HSP per unordered pair by max bits (ties -> lower e), then the dual
    threshold. `scale` multiplies e-values to renormalize the search space; v6 uses
    1.0 (threshold the e-value as reported at this run's db)."""
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
    passing = {k: v for k, v in best.items() if v[0] >= ID_MIN and v[1] < EVALUE_MAX}
    return passing, total, len(best)


def partition(prefix, scale=1.0):
    nodes = {r["node_id"]: r for r in
             csv.DictReader(Path(f"{prefix}_nodes.tsv").open(), delimiter="\t")}
    pair_files = [f"{prefix}_pairs.tsv", f"{prefix}_rev_pairs.tsv"]
    edges, n_raw, n_pairs = load_pairs(pair_files, scale=scale)
    edge_list = [(q, t, *v) for (q, t), v in sorted(edges.items())]

    groups = components(set(nodes), edge_list)

    comps = []
    for root, mem in groups.items():
        comps.append({"members": mem, "canon": min(mem, key=lambda n: canon_key(nodes, n))})
    comps.sort(key=lambda c: canon_key(nodes, c["canon"]))
    for rank, c in enumerate(comps, start=1):
        c["component_id"] = f"C{rank:03d}"
        c["canonical_label"] = node_label(nodes[c["canon"]])
    return nodes, edge_list, comps, {"raw_alignments": n_raw, "aligned_pairs": n_pairs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="stem for _nodes.tsv/_pairs.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--date", default="2026-07-17")
    ap.add_argument("--tag", default="v6")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="e-value renorm multiplier (v6=1.0: threshold at the run db)")
    args = ap.parse_args()

    nodes, edge_list, comps, frag = partition(args.prefix, scale=args.scale)
    outdir = Path(args.outdir)
    stem = f"{args.tag}_{args.date}"

    key = lambda n: canon_key(nodes, n)

    comp_of = {n: c for c in comps for n in c["members"]}
    for sr, c in enumerate(sorted(comps, key=lambda c: (-len(c["members"]),
                                                        key(c["canon"]))), start=1):
        c["size_rank"] = sr

    N = sum(int(nodes[n]["sequence_length"]) for n in nodes)

    # ---- edges
    with (outdir / f"component_edges_{stem}.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["node_a", "node_b", "pl_a", "pl_b", "pct_id", "evalue", "bits",
                    "acc_a", "acc_b", "component_id"])
        for q, t, pid, ev, bits in edge_list:
            w.writerow([q, t, nodes[q]["identifier"], nodes[t]["identifier"],
                        f"{pid:.1f}", f"{ev:.3e}", f"{bits:.1f}",
                        nodes[q]["accession"], nodes[t]["accession"],
                        comp_of[q]["component_id"]])

    # ---- per-node assignment (one row per node)
    with (outdir / f"component_assignment_{stem}.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["node_id", "sequence_md5", "identifier", "identifier_all",
                    "accession", "accession_all", "pazy_id", "component_id",
                    "size_rank", "sequence_length", "enzyme_name", "organism",
                    "v1_component", "cath"])
        for n in sorted(nodes, key=lambda n: (comp_of[n]["size_rank"], key(n))):
            r, c = nodes[n], comp_of[n]
            w.writerow([n, r["sequence_md5"], r["identifier"], r["identifier_all"],
                        r["accession"], r["accession_all"], r["pazy_id"],
                        c["component_id"], c["size_rank"], r["sequence_length"],
                        r["enzyme_name"], r["organism"], r["v1_component"], r["cath"]])

    # ---- per-component members
    with (outdir / f"component_members_{stem}.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["component_id", "size_rank", "canonical_label", "n_members",
                    "n_with_v1_overlay", "v1_component_spread", "cath_folds",
                    "member_identifiers"])
        for c in sorted(comps, key=lambda c: c["size_rank"]):
            mem = sorted(c["members"], key=key)
            v1labs = sorted({nodes[n]["v1_component"] for n in mem
                             if nodes[n]["v1_component"]}, key=lambda s: int(s))
            cath = sorted({nodes[n]["cath"] for n in mem if nodes[n]["cath"]})
            w.writerow([c["component_id"], c["size_rank"], c["canonical_label"],
                        len(mem), sum(1 for n in mem if nodes[n]["v1_component"]),
                        ";".join(v1labs), ";".join(cath),
                        ";".join(node_label(nodes[n]) for n in mem)])

    sizes = sorted((len(c["members"]) for c in comps), reverse=True)
    v1_spread = defaultdict(set)
    for c in comps:
        for n in c["members"]:
            if nodes[n]["v1_component"]:
                v1_spread[nodes[n]["v1_component"]].add(c["component_id"])
    v1_split = {lab: sorted(cs) for lab, cs in v1_spread.items() if len(cs) > 1}

    stats = {
        "tag": args.tag, "date": args.date, "prefix": str(args.prefix),
        "partition_mode": "flat de novo (no v1 anchor)",
        "evalue_scale": args.scale, "db_letters_N": N,
        "usearch": "v11.0.667 -allpairs_local -acceptall, both orientations",
        "docker_platform": "linux/amd64",
        "edge_criterion": f"pct_id>={ID_MIN} AND evalue<{EVALUE_MAX:g} @ {N}-residue db",
        "n_nodes": len(nodes),
        "raw_alignments": frag["raw_alignments"], "aligned_pairs": frag["aligned_pairs"],
        "edges_passing": len(edge_list),
        "n_components": len(comps),
        "n_singletons": sum(1 for c in comps if len(c["members"]) == 1),
        "largest_component": sizes[0] if sizes else 0,
        "top10_component_sizes": sizes[:10],
        "n_nodes_with_v1_overlay": sum(1 for n in nodes if nodes[n]["v1_component"]),
        "n_v1_labels_split_across_components": len(v1_split),
        "v1_labels_split": v1_split,
    }
    (outdir / f"stats_{stem}.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
