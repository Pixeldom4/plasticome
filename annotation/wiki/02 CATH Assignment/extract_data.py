"""Pull the figure and demo data for the CATH Assignment draft out of the
plasticome-cath-assignment results. Writes data.json next to this script.

  python3 extract_data.py ~/plasticome-cath-assignment

Reads only committed result tables, plus results/route_a_ted50/foldseek/aln.tsv
(gitignored, 42 MB) for the E-value calibration figure; that block is skipped if
the file is absent. Reference of record is CATH v4.3.0 throughout; the v4.4.0
audit is summarised from results/cath_release_comparison.tsv only.
"""
import csv, json, sys, math, pathlib, collections as C

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "~/plasticome-cath-assignment").expanduser()
R = root / "results"
rd = lambda p: list(csv.DictReader(open(p, newline=""), delimiter="\t"))
blank = lambda v: v in ("", "-", "NA", None)
out = {}

# CATH v4.3.0 superfamily names
names = {}
for line in open(root / "v1-reference" / "cath-names.txt"):
    if line.startswith("#"): continue
    p = line.rstrip("\n").split(None, 2)
    if len(p) >= 3: names[p[0]] = p[2].lstrip(":").strip()

# ---- per enzyme ----
prov = {r["plasticome_id"]: r for r in rd(R / "union_cath_provenance.tsv")}
enz = []
for r in rd(R / "plasticome_three_route_annotation.tsv"):
    p = prov[r["plasticome_id"]]
    enz.append({"id": int(r["plasticome_id"]), "name": r["enzyme_name"], "acc": r["accession"], "src": r["source"],
                "len": int(r["seq_len"]), "comp": p["component"], "cat": p["catalytic-cath"], "acc_cath": p["accessory-cath"],
                "dom": p["domain"], "tier": p["tier"], "cent": p["source_centroid_PL"], "pid": p["pident"],
                "A": r["A_s40_cath"], "T": r["A_ted_cath"], "Tl": r["A_ted_label"], "B": r["B_gene3d_cath"],
                "Cn": r["C_pfam_name"], "Ca": r["C_pfam_acc"], "Cc": r["C_pfam_clan_id"],
                "nr": int(r["n_routes_with_hit"]), "agree": r["cath_route_agreement"]})
out["enz"] = enz

# ---- per domain (Merizo layout; Pfam/Gene3D reconciled onto it) ----
dom = C.defaultdict(list)
for r in rd(R / "plasticome_three_route_domains.tsv"):
    if r["domain_id"] in ("", "-"): continue
    dom[r["plasticome_id"]].append({
        "d": r["domain_id"].split("_", 1)[1], "seg": r["segments"],
        "A": r["A_s40_cath"], "Atm": r["A_s40_tm"], "Aev": r["A_s40_evalue"], "As": r["A_s40_status"],
        "T": r["A_ted_cath"], "Tl": r["A_ted_label"], "Ttm": r["A_ted_tm"], "Ts": r["A_ted_status"],
        "B": r["B_gene3d_cath"], "Bsc": r["B_gene3d_score"], "Bseg": r["B_gene3d_segments"],
        "C": r["C_pfam_name"], "Ca": r["C_pfam_acc"], "Cc": r["C_pfam_clan_id"], "Cenv": r["C_pfam_env"],
        "agr": r["domain_agreement"]})
out["dom"] = dom
out["dom_agree"] = C.Counter(d["agr"] for v in dom.values() for d in v)
out["n_merizo_domains"] = sum(len(v) for v in dom.values())
out["enz_agree"] = C.Counter(e["agree"] for e in enz)
out["n_routes"] = C.Counter(e["nr"] for e in enz)

# route coverage, enzymes with >= 1 call
cov = {k: sum(not blank(e[k]) for e in enz) for k in ("A", "T", "B", "Ca")}
cov["union"] = sum(any(not blank(e[k]) for k in ("A", "T", "B", "Ca")) for e in enz)
cov["cath_union"] = sum(any(not blank(e[k]) for k in ("A", "T", "B")) for e in enz)
out["cov"] = cov

# route-level domain counts
ra = rd(R / "route_a" / "route_a_domains.tsv")
rt = rd(R / "route_a_ted50" / "route_a_domains.tsv")
rb = rd(R / "route_b" / "route_b_domains.tsv")
rc = rd(R / "route_c_pfam" / "route_c_domains.tsv")
out["routes"] = {
    "A": {"domains": len(ra), "hit": sum(r["status"] == "hit" for r in ra)},
    "T": {"domains": len(rt), "status": C.Counter(r["status"] for r in rt)},
    "B": {"domains": len(rb), "seqs": len({r["structure_id"] for r in rb})},
    "C": {"domains": len(rc), "seqs": len({r[next(k for k in rc[0] if "struct" in k or "id" in k)] for r in rc})},
}
ndom = C.Counter(int(r["structure_id"]) for r in ra)
out["merizo_ndom"] = C.Counter(ndom.get(i, 0) for i in range(1, 610))
out["merizo_discontinuous"] = sum(";" in r["segments"] or "_" in r["segments"] for r in ra)

# ---- A vs B cross-check ----
out["xc"] = {r["tier"]: int(r["count"]) for r in rd(R / "crosscheck" / "summary.tsv")}
out["ted_cmp"] = C.Counter(r["verdict"] for r in rd(R / "route_a_ted50" / "structure_comparison.tsv"))

# ---- fold distribution: final catalytic CATH per union row ----
folds = C.Counter()
for e in enz:
    for s in e["cat"].split(";"):
        folds[s.strip()] += 1
out["folds"] = [[s, names.get(s, ""), n] for s, n in folds.most_common()]
acc = C.Counter(s.strip() for e in enz for s in e["acc_cath"].split(";") if not blank(s.strip()))
out["accessory"] = [[s, names.get(s, ""), n] for s, n in acc.most_common()]
used = {s for e in enz for k in ("A", "T", "B", "cat", "acc_cath") for s in e[k].split(";")}
out["names"] = {s: names[s] for s in used if s in names}

# ---- centroid level and transfer ----
cen = rd(R / "centroids_cath_assignment.tsv")
out["centroids"] = {"n": len(cen), "assigned": sum(r["cath_assigned"] == "Y" for r in cen),
                    "sources": C.Counter(r["sources"] if not blank(r["sources"]) else "none" for r in cen),
                    "components": len({r["component"] for r in cen})}
out["union_tier"] = C.Counter(p["tier"] for p in prov.values())
comp = rd(R / "component_representative_cath.tsv")
out["components"] = {"n": len(comp), "homogeneous": sum(r["homogeneous"] == "Y" for r in comp)}

# ---- validation and audit ----
out["v1"] = C.Counter(r["category"] for r in rd(R / "crosscheck" / "v1_validation.tsv"))
out["release"] = C.Counter(r["verdict"] for r in rd(R / "cath_release_comparison.tsv"))

# ---- TED50 E-value calibration (raw HSPs) ----
aln = R / "route_a_ted50" / "foldseek" / "aln.tsv"
if aln.exists():
    LO, HI, W = -60, 2, 1
    hist = [0] * ((HI - LO) // W)
    best = {}  # domain -> best E among TM>=0.5 HSPs with a 4-level superfamily label
    n = 0
    for line in open(aln):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12: continue
        tm, ev = float(f[10]), max(float(f[4]), 1e-300)
        if tm < 0.5: continue
        n += 1
        le = math.log10(ev)
        hist[min(max(int((le - LO) // W), 0), len(hist) - 1)] += 1
        lab = f[11].split(" CATH", 1)[1].split("_RES")[0] if " CATH" in f[11] else ""
        if lab.count(".") == 3:
            q = f[0].replace(".pdb", "")
            best[q] = min(best.get(q, 1e9), ev)
    out["ted_hist"] = {"lo": LO, "w": W, "counts": hist, "n": n}
    out["ted_yield"] = [[c, sum(v <= c for v in best.values())] for c in (1e-20, 1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1)]

p = pathlib.Path(__file__).with_name("data.json")
p.write_text(json.dumps(out, separators=(",", ":")))
print(p, p.stat().st_size, "bytes")
print(json.dumps({k: v for k, v in out.items() if k not in ("enz", "dom", "names")}, indent=1)[:4000])
