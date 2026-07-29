#!/usr/bin/env python3
"""Post-filter DIAMOND all-vs-all hits and build single-linkage components.

Paper-faithful criterion carried from the usearch pipeline (v4-v8):
    keep an edge iff  pct_id >= 30.0  AND  evalue < 1e-5
Self-hits (q == t) dropped. Filtering happens HERE, never as a search parameter,
exactly like usearch's `-acceptall` + post-filter. Nodes are the m#### ids from the
shared usearch node set, so output is directly comparable to the usearch partition.
Emits both a per-node assignment CSV and the passing edge list.
"""
from collections import defaultdict
from pathlib import Path
import csv, json

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
PAIRS = OUT / "allpairs.tsv"
FASTA = OUT / "seqs.faa"
ID_MIN, EVALUE_MAX = 30.0, 1e-5

# best HSP per unordered pair by max bits (ties -> lower e), mirroring usearch load_pairs
best = {}
n_raw = 0
with PAIRS.open() as f:
    for line in f:
        q, s, pid, ev, bit = line.rstrip("\n").split("\t")
        n_raw += 1
        if q == s:
            continue
        pid, ev, bit = float(pid), float(ev), float(bit)
        k = (q, s) if q < s else (s, q)
        cur = best.get(k)
        if cur is None or bit > cur[2] or (bit == cur[2] and ev < cur[1]):
            best[k] = (pid, ev, bit)

edges = {k: v for k, v in best.items() if v[0] >= ID_MIN and v[1] < EVALUE_MAX}

ids = [l[1:].strip() for l in FASTA.read_text().splitlines() if l.startswith(">")]
parent = {i: i for i in ids}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
for (a, b) in edges:
    parent[find(a)] = find(b)

comp = defaultdict(list)
for i in ids:
    comp[find(i)].append(i)

def mnum(x):
    return int(x[1:]) if x[1:].isdigit() else 10**9
ordered = sorted(comp.values(), key=lambda v: (-len(v), min(mnum(x) for x in v)))

# per-node assignment
with (OUT / "component_assignment.csv").open("w", newline="") as w:
    wr = csv.writer(w)
    wr.writerow(["node_id", "component", "component_size"])
    for cid, members in enumerate(ordered, 1):
        for m in sorted(members, key=mnum):
            wr.writerow([m, cid, len(members)])

# passing edge list
with (OUT / "component_edges.csv").open("w", newline="") as w:
    wr = csv.writer(w)
    wr.writerow(["node_a", "node_b", "pct_id", "evalue", "bits"])
    for (a, b), (pid, ev, bit) in sorted(edges.items()):
        wr.writerow([a, b, f"{pid:.1f}", f"{ev:.3e}", f"{bit:.1f}"])

sizes = [len(v) for v in ordered]
stats = {
    "tool": "diamond blastp --ultra-sensitive (Docker linux/amd64)",
    "edge_criterion": f"pct_id>={ID_MIN} AND evalue<{EVALUE_MAX:g}",
    "n_nodes": len(ids),
    "raw_alignments": n_raw,
    "aligned_pairs": len(best),
    "edges_passing": len(edges),
    "n_components": len(ordered),
    "n_singletons": sum(1 for v in ordered if len(v) == 1),
    "largest_component": sizes[0] if sizes else 0,
    "top10_component_sizes": sizes[:10],
}
(OUT / "stats_diamond.json").write_text(json.dumps(stats, indent=2))

print(f"raw alignments   : {n_raw}")
print(f"aligned pairs    : {len(best)}")
print(f"edges kept       : {len(edges)}")
print(f"sequences        : {len(ids)}")
print(f"components       : {len(ordered)}")
print(f"singletons       : {stats['n_singletons']}")
print(f"top sizes        : {sizes[:10]}")
