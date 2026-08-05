#!/usr/bin/env python3
"""Steps 2-3 - flat similarity graph + connected-component partition.

Edge criterion (post-processing only, NEVER a search parameter):
    >= 30% amino-acid identity AND e-value < 1e-5
Partition: undirected single-linkage connected components over ALL nodes.

Flat de-novo run (decisions 1-4), identical method to v4-v8:
  * No v1 anchor, no delta, no join/merge/novel typing.
  * E-values thresholded AS REPORTED by the aligner at this run's search space (no
    down-scaling to the paper's 213-seq space; decision 2). --scale 1.0 = identity.
  * HSP selection: best HSP per unordered pair by MAX BITS (ties -> lower e-value),
    over every pair file present (adopted v3 s7 convention; accepted B4 cost). That
    is one DIAMOND file, which already holds both directions of each pair, or the
    two usearch orientation files, whose asymmetry is what the reduction repairs.
  * Canonical labels/IDs resolve on the PL identifier when present, else fall back
    to accession, enzyme name, node id. The 213-node Step-0 roster has no PL
    identifiers, so it labels by accession -- the component/edge COUNTS are label-
    independent, so the 42 / 3,178 sanity target is unaffected.

The v1_component/cath columns are a descriptive OVERLAY only (Step 3 reporting),
never used to seed or gate the partition.
"""
import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402

ID_MIN = config.ID_MIN
EVALUE_MAX = config.EVALUE_MAX

# What produced the pair files, recorded in stats_*.json so a partition can always be
# traced back to its aligner (the edge sets differ slightly between the two).
ALIGNER_PROVENANCE = {
    "diamond": "diamond blastp all-vs-all, " + " ".join(config.DIAMOND_FLAGS),
    "usearch": f"usearch {config.USEARCH_VERSION} -allpairs_local -acceptall, both orientations",
}


def pl_num(ident):
    m = re.search(r"(\d+)", ident or "")
    return int(m.group(1)) if m else 10**9


def node_label(node):
    """Canonical label: PL identifier, else accession, else enzyme name, else node id."""
    return (node.get("identifier") or node.get("accession")
            or node.get("enzyme_name") or node["node_id"])


def canon_key(nodes, n):
    """Canonical ordering: by PL number (stable, unique when present), then label,
    node_id as final tiebreak. Rosters without PL fall to (10**9, label, node_id)
    == accession order."""
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
    threshold. `scale` multiplies e-values to renormalize the search space (1.0 =
    threshold the e-value as reported at this run's db).

    A path that does not exist is skipped: the reversed-orientation file is a usearch
    device (-allpairs_local is asymmetric in query/target), and DIAMOND reports both
    directions of every pair inside its single all-vs-all run, so it writes only
    `_pairs.tsv`. The max-bits reduction below is what makes either input symmetric.
    """
    best = {}
    total = 0
    for path in paths:
        if not Path(path).exists():
            continue
        for line in Path(path).open():
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                q, t, pid, ev, bits = f[0], f[1], float(f[2]), float(f[3]), float(f[4])
            except ValueError:
                continue  # skip a labelled header row (query/target/pct_id/...)
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
    pair_files = [p for p in (f"{prefix}_pairs.tsv", f"{prefix}_rev_pairs.tsv")
                  if Path(p).exists()]
    if not pair_files:
        sys.exit(f"error: no alignment output at {prefix}_pairs.tsv (run step2_align.py first)")
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
    return nodes, edge_list, comps, {"raw_alignments": n_raw, "aligned_pairs": n_pairs,
                                     "pair_files": [Path(p).name for p in pair_files]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prefix", required=True, help="stem for _nodes.tsv/_pairs.tsv")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--date", default=config.DATE)
    ap.add_argument("--tag", default=config.TAG)
    ap.add_argument("--scale", type=float, default=config.EVALUE_SCALE,
                    help="e-value renorm multiplier (1.0: threshold at the run db)")
    ap.add_argument("--engine", default=config.ALIGNER, choices=["diamond", "usearch"],
                    help="aligner that produced the pair files; recorded in stats only")
    args = ap.parse_args()

    nodes, edge_list, comps, frag = partition(args.prefix, scale=args.scale)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.tag}_{args.date}"
    key = lambda n: canon_key(nodes, n)

    comp_of = {n: c for c in comps for n in c["members"]}
    for sr, c in enumerate(sorted(comps, key=lambda c: (-len(c["members"]), key(c["canon"]))),
                           start=1):
        c["size_rank"] = sr

    N = sum(int(nodes[n]["sequence_length"]) for n in nodes)

    # ---- edges
    with (outdir / f"component_edges_{stem}.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["node_a", "node_b", "pl_a", "pl_b", "pct_id", "evalue", "bits",
                    "acc_a", "acc_b", "component_id"])
        for q, t, pid, ev, bits in edge_list:
            w.writerow([q, t, nodes[q].get("identifier", ""), nodes[t].get("identifier", ""),
                        f"{pid:.1f}", f"{ev:.3e}", f"{bits:.1f}",
                        nodes[q].get("accession", ""), nodes[t].get("accession", ""),
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
            w.writerow([n, r["sequence_md5"], r.get("identifier", ""),
                        r.get("identifier_all", ""), r.get("accession", ""),
                        r.get("accession_all", ""), r.get("pazy_id", ""),
                        c["component_id"], c["size_rank"], r["sequence_length"],
                        r.get("enzyme_name", ""), r.get("organism", ""),
                        r.get("v1_component", ""), r.get("cath", "")])

    # ---- per-component members
    with (outdir / f"component_members_{stem}.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["component_id", "size_rank", "canonical_label", "n_members",
                    "n_with_v1_overlay", "v1_component_spread", "cath_folds",
                    "member_identifiers"])
        for c in sorted(comps, key=lambda c: c["size_rank"]):
            mem = sorted(c["members"], key=key)
            v1labs = sorted({nodes[n]["v1_component"] for n in mem if nodes[n].get("v1_component")},
                            key=lambda s: int(s))
            cath = sorted({nodes[n]["cath"] for n in mem if nodes[n].get("cath")})
            w.writerow([c["component_id"], c["size_rank"], c["canonical_label"], len(mem),
                        sum(1 for n in mem if nodes[n].get("v1_component")),
                        ";".join(v1labs), ";".join(cath),
                        ";".join(node_label(nodes[n]) for n in mem)])

    sizes = sorted((len(c["members"]) for c in comps), reverse=True)
    v1_spread = defaultdict(set)
    for c in comps:
        for n in c["members"]:
            if nodes[n].get("v1_component"):
                v1_spread[nodes[n]["v1_component"]].add(c["component_id"])
    v1_split = {lab: sorted(cs) for lab, cs in v1_spread.items() if len(cs) > 1}

    stats = {
        "tag": args.tag, "date": args.date, "prefix": str(args.prefix),
        "partition_mode": "flat de novo (no v1 anchor)",
        "evalue_scale": args.scale, "db_letters_N": N,
        "aligner": args.engine,
        "aligner_command": ALIGNER_PROVENANCE.get(args.engine, args.engine),
        "pair_files": frag["pair_files"],
        "docker_platform": config.DOCKER_PLATFORM,
        "edge_criterion": f"pct_id>={ID_MIN} AND evalue<{EVALUE_MAX:g} @ {N}-residue db",
        "n_nodes": len(nodes),
        "raw_alignments": frag["raw_alignments"], "aligned_pairs": frag["aligned_pairs"],
        "edges_passing": len(edge_list),
        "n_components": len(comps),
        "n_singletons": sum(1 for c in comps if len(c["members"]) == 1),
        "largest_component": sizes[0] if sizes else 0,
        "top10_component_sizes": sizes[:10],
        "n_nodes_with_v1_overlay": sum(1 for n in nodes if nodes[n].get("v1_component")),
        "n_v1_labels_split_across_components": len(v1_split),
        "v1_labels_split": v1_split,
    }
    (outdir / f"stats_{stem}.json").write_text(json.dumps(stats, indent=2))
    print(f"[step23:{args.tag}] {len(nodes)} nodes -> {len(comps)} components, "
          f"{len(edge_list)} edges (largest {sizes[0] if sizes else 0})")


if __name__ == "__main__":
    main()
