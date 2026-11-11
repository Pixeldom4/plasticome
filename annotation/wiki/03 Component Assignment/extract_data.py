"""Pull the figure and demo data for the Component Assignment draft out of the
final plasticome run. Writes data.json next to this script.

  python3 extract_data.py [run dir] [2x2 archive dir]

Defaults: runs/2026-08-06.final-usearch.3 and archive/diamond-usearch-comparison.
Per-component CATH comes from ../02 CATH Assignment/data.json (catalytic CATH per
enzyme). The identity sweep re-partitions the run's raw HSPs with the same rule as
lib/03 alignment/scripts/step23_graph.py (best HSP per unordered pair by max bits,
ties to lower E; edge if pct_id >= cutoff and E < 1e-5; single linkage).
"""
import csv, json, sys, pathlib, collections as C

here = pathlib.Path(__file__).parent
repo = here.parents[2]          # annotation/wiki/<page>/ -> repo root
run = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else repo / "runs" / "2026-08-06.final-usearch.3"
arc = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else repo / "archive" / "diamond-usearch-comparison"
I = run / "03-alignment.intermediates"
rd = lambda p, d="\t": list(csv.DictReader(open(p, newline=""), delimiter=d))
blank = lambda v: v in ("", "-", "NA", None)
out = {}

# ---- run stats ----
st = json.loads(next(I.glob("stats_*.json")).read_text())
out["stats"] = {k: st[k] for k in ("n_nodes", "raw_alignments", "aligned_pairs", "edges_passing", "n_components",
                                   "n_singletons", "largest_component", "db_letters_N", "aligner_command",
                                   "n_nodes_with_v1_overlay", "n_v1_labels_split_across_components")}

# ---- per enzyme and per component ----
cath = json.loads((here.parent / "02 CATH Assignment" / "data.json").read_text())
cat = {str(e["id"]): e for e in cath["enz"]}
U = rd(run / "05-union-with-components.tsv")
out["union"] = {"n": len(U), "centroids": sum(r["is_centroid"] == "yes" for r in U)}
comp = C.defaultdict(lambda: {"cent": 0, "enz": 0, "cath": C.Counter(), "dom": C.Counter(), "names": []})
for r in U:
    c = comp[r["component_id"]]; e = cat[r["plasticome_id"]]
    c["enz"] += 1; c["cent"] += r["is_centroid"] == "yes"
    c["cath"][e["cat"]] += 1; c["dom"][e["dom"].split(";")[0]] += 1
    if r["is_centroid"] == "yes": c["names"].append(r["enzyme_name"])
out["comps"] = [{"c": k, "cent": v["cent"], "enz": v["enz"], "cath": v["cath"].most_common(1)[0][0],
                 "dom": v["dom"].most_common(1)[0][0], "pure": len(v["cath"]) == 1, "ex": v["names"][:4]}
                for k, v in sorted(comp.items(), key=lambda kv: (-kv[1]["cent"], -kv[1]["enz"], kv[0]))]
out["size_hist"] = C.Counter(c["cent"] for c in out["comps"])
out["names"] = {c["cath"]: cath["names"].get(c["cath"], "") for c in out["comps"] if not blank(c["cath"])}

# ---- identity sweep over the run's raw HSPs ----
nodes = rd(I / "combined_nodes.tsv")
asg = {r["node_id"]: r["component_id"] for r in rd(next(I.glob("component_assignment_*.csv")), ",")}
by_pl = {}
for r in nodes:
    pl = r["identifier"].split("|")[0]  # CL0001: cluster id
    by_pl[r["node_id"]] = pl
cl_rep = {r["cluster_id"]: r["rep_plasticome_id"] for r in rd(run / "03-alignment.tsv")}
node_cat = {n: cat[cl_rep[str(int(by_pl[n][2:]))]]["cat"] for n in by_pl}
best = {}
for f in ("combined_pairs.tsv", "combined_rev_pairs.tsv"):
    for line in open(I / f):
        q, t, pid, ev, bits = line.rstrip("\n").split("\t")[:5]
        if q == t: continue
        k = (q, t) if q < t else (t, q); pid, ev, bits = float(pid), float(ev), float(bits)
        cur = best.get(k)
        if cur is None or bits > cur[2] or (bits == cur[2] and ev < cur[1]): best[k] = (pid, ev, bits)

def partition(cut):
    par = {n: n for n in by_pl}
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for (a, b), (pid, ev, _) in best.items():
        if pid >= cut and ev < 1e-5: par[find(a)] = find(b)
    g = C.defaultdict(list)
    for n in by_pl: g[find(n)].append(n)
    return list(g.values())

sweep = []
for cut in (20, 25, 30, 35, 40, 50, 60, 70):
    P = partition(cut)
    # purity: centroids whose catalytic CATH equals their component's most common one
    pure = sum(C.Counter(node_cat[n] for n in p).most_common(1)[0][1] for p in P)
    mixed = sum(len(set(node_cat[n] for n in p if not blank(node_cat[n]))) > 1 for p in P)
    sweep.append({"cut": cut, "comps": len(P), "single": sum(len(p) == 1 for p in P), "largest": max(map(len, P)),
                  "purity": round(100 * pure / len(by_pl), 1), "mixed": mixed})
out["sweep"] = sweep
assert next(s for s in sweep if s["cut"] == 30)["comps"] == st["n_components"], "sweep does not reproduce the run"

# identity of best HSP per pair, among pairs that pass the E-value half of the rule
hist = [0] * 20  # 5% bins
for pid, ev, _ in best.values():
    if ev < 1e-5: hist[min(int(pid // 5), 19)] += 1
out["pid_hist"] = {"w": 5, "counts": hist, "n": sum(hist), "pairs": len(best)}

# ---- USEARCH / DIAMOND 2x2 ----
fac = []
for d in sorted(arc.glob("clust-*.align-*")):
    s = json.loads(next((d / "03-alignment.intermediates").glob("stats_*.json")).read_text())
    cl, al = d.name.split(".")
    fac.append({"clust": cl.split("-")[1], "align": al.split("-")[1], "nodes": s["n_nodes"], "pairs": s["aligned_pairs"],
                "edges": s["edges_passing"], "comps": s["n_components"], "single": s["n_singletons"],
                "largest": s["largest_component"]})
out["factorial"] = fac

# ---- component map (Figure 1), from figures/ssn/out ----
# Node positions are build_ssn.py's packed layout. Halos are redrawn in the page
# the way draw_blob() does it: a disc round every node plus a narrower bridge
# along a Euclidean minimum spanning tree of each map component.
S = repo / "figures" / "ssn" / "out"
sn = rd(S / "ssn_nodes_c30_d40_k3.tsv")
ix = {r["plasticome_id"]: i for i, r in enumerate(sn)}
upl = {r["plasticome_id"]: r for r in U}
seeds = {}
by_seq = {r["aa_sequence"].strip().upper(): r["plasticome_id"] for r in U}
for s in rd(repo / "source-data" / "plasticome-curated-seeds" / "plasticome-curated-seeds.tsv"):
    pid = next((r["plasticome_id"] for r in sn if r["accession"].strip() == s["accession"].strip()), None) \
        or by_seq.get(s["aa_sequence"].strip().upper())
    if pid in ix: seeds[pid] = s["enzyme_name"].strip()
nodes = [{"id": int(r["plasticome_id"]), "name": r["enzyme_name"], "x": round(float(r["x"]), 3), "y": round(float(r["y"]), 3),
          "sub": r["substrate_group"], "raw": r["substrate_raw"], "shape": r["shape"], "fam": r["component_family"],
          "m": int(r["component"]), "c": upl[r["plasticome_id"]]["component_id"], "dom": r["domain"],
          "seed": seeds.get(r["plasticome_id"], "")} for r in sn]
bridges = []
for m in {n["m"] for n in nodes}:
    mem = [i for i, n in enumerate(nodes) if n["m"] == m]
    if len(mem) < 2: continue
    d2 = lambda a, b: (nodes[a]["x"] - nodes[b]["x"]) ** 2 + (nodes[a]["y"] - nodes[b]["y"]) ** 2
    inside, rest = {mem[0]}, set(mem[1:])
    best_to = {j: (d2(mem[0], j), mem[0]) for j in rest}
    while rest:  # Prim
        j = min(rest, key=lambda k: best_to[k][0]); bridges.append([best_to[j][1], j])
        rest.discard(j); inside.add(j)
        for k in rest:
            dk = d2(j, k)
            if dk < best_to[k][0]: best_to[k] = (dk, j)
edges = [[ix[r["source"]], ix[r["target"]], round(float(r["pident"]), 1)] for r in rd(S / "ssn_edges_c30_d40_k3.tsv") if r["drawn"] == "1"]
out["map"] = {"nodes": nodes, "edges": edges, "bridges": bridges, "pad": 0.95}

p = here / "data.json"
p.write_text(json.dumps(out, separators=(",", ":")))
print(p, p.stat().st_size, "bytes")
print(json.dumps({k: v for k, v in out.items() if k not in ("comps", "names")}, indent=1))
