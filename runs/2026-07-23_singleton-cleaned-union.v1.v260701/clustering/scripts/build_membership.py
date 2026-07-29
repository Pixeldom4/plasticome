import csv, collections, os

# resolve paths relative to the run dir (parent of this scripts/ folder), CWD-independent
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def p(*parts): return os.path.join(ROOT, *parts)

# --- read seed labels (reference centroids), preserve file order for stable ids ---
seeds = []
for line in open(p("inputs", "curated_seed_set.fasta")):
    if line.startswith(">"):
        seeds.append(line[1:].strip())
seed_idx = {s: i for i, s in enumerate(seeds)}

rows = []  # (cluster_id, origin, role, label, centroid, pct_id)

# --- reference clusters: seed = centroid ---
ref_members = collections.defaultdict(list)  # seed -> [(member,pctid)]
for line in open(p("work", "02-closedref.uc")):
    f = line.rstrip("\n").split("\t")
    if f[0] == "H":
        query, target, pid = f[8], f[9], f[3]
        ref_members[target].append((query, pid))

for s in seeds:
    cid = f"R{seed_idx[s]:03d}"
    rows.append((cid, "reference", "centroid", s, s, "*"))
    for m, pid in ref_members.get(s, []):
        rows.append((cid, "reference", "member", m, s, pid))

# --- de novo clusters ---
# parse S/H; S establishes centroid for cluster number, H are members
denovo_centroid = {}   # clusternum -> centroid label
denovo_members = collections.defaultdict(list)
for line in open(p("work", "03-denovo.uc")):
    f = line.rstrip("\n").split("\t")
    if f[0] == "S":
        denovo_centroid[int(f[1])] = f[8]
    elif f[0] == "H":
        denovo_members[int(f[1])].append((f[8], f[3], f[9]))

for cnum in sorted(denovo_centroid):
    cid = f"D{cnum:03d}"
    cen = denovo_centroid[cnum]
    rows.append((cid, "de_novo", "centroid", cen, cen, "*"))
    for m, pid, tgt in denovo_members.get(cnum, []):
        rows.append((cid, "de_novo", "member", m, cen, pid))

# --- compute sizes and write, sorted by size desc then cluster id ---
size = collections.Counter(r[0] for r in rows)
rows.sort(key=lambda r: (-size[r[0]], r[0], r[2] != "centroid"))

os.makedirs(p("results"), exist_ok=True)
with open(p("results", "cluster_membership.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["cluster_id", "origin", "size", "role", "label", "centroid", "pct_id_to_centroid"])
    for r in rows:
        w.writerow([r[0], r[1], size[r[0]], r[2], r[3], r[4], r[5]])

# --- summary ---
clusters = sorted(size)
ref_clusters = [c for c in clusters if c.startswith("R")]
den_clusters = [c for c in clusters if c.startswith("D")]
print(f"total records written : {len(rows)}")
print(f"total clusters        : {len(clusters)}  (reference {len(ref_clusters)} + de_novo {len(den_clusters)})")
print(f"reference cluster sizes:")
for c in sorted(ref_clusters, key=lambda c:-size[c]):
    cen = next(r[3] for r in rows if r[0]==c and r[3] in seed_idx)
    print(f"   {c}  size={size[c]:>2}  {cen}")
sd = collections.Counter(size[c] for c in clusters)
print("cluster size distribution (size:count):", dict(sorted(sd.items())))
singl = sum(1 for c in clusters if size[c]==1)
print(f"singletons            : {singl} ({100*singl/len(clusters):.1f}% of clusters)")
print(f"largest cluster       : {max(size.values())} members -> {max(size, key=size.get)}")
