#!/usr/bin/env python3
"""Post-filter DIAMOND all-vs-all hits and build single-linkage components.

Paper-faithful criterion carried from the usearch pipeline (v4-v8):
    keep an edge iff  pct_id >= 30.0  AND  evalue < 1e-5
Self-hits (q == t) are dropped (they are self-loops). Filtering is done HERE,
never as a search parameter, exactly like usearch's `-acceptall` + post-filter.
Writes component_assignment.csv next to this script.
"""
from collections import defaultdict
from pathlib import Path
import csv

HERE = Path(__file__).resolve().parent
PAIRS = HERE / "allpairs.tsv"
FASTA = HERE / "seqs.faa"
OUT = HERE / "component_assignment.csv"
ID_MIN, EVALUE_MAX = 30.0, 1e-5

edges = set()
with PAIRS.open() as f:
    for q, s, pid, ev, _bit in (l.rstrip("\n").split("\t") for l in f):
        if q != s and float(pid) >= ID_MIN and float(ev) < EVALUE_MAX:
            edges.add(frozenset((q, s)))

ids = [l[1:].strip() for l in FASTA.read_text().splitlines() if l.startswith(">")]
parent = {i: i for i in ids}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
for e in edges:
    a, b = tuple(e)
    parent[find(a)] = find(b)

comp = defaultdict(list)
for i in ids:
    comp[find(i)].append(i)
# stable component ids: largest first
ordered = sorted(comp.values(), key=lambda v: (-len(v), min(int(x) for x in v)))
with OUT.open("w", newline="") as w:
    wr = csv.writer(w)
    wr.writerow(["identifier", "component", "component_size"])
    for cid, members in enumerate(ordered, 1):
        for m in members:
            wr.writerow([m, cid, len(members)])

print(f"edges kept       : {len(edges)}")
print(f"sequences        : {len(ids)}")
print(f"components       : {len(ordered)}")
print(f"singletons       : {sum(1 for v in ordered if len(v) == 1)}")
print(f"top sizes        : {[len(v) for v in ordered[:8]]}")
print(f"wrote            : {OUT}")
