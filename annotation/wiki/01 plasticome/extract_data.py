"""Pull the figure and demo data for the Plasticome Reference Curation draft out of
plasticome.v2.xlsx. Writes data.json next to this script.

Run with the plasticome conda env (needs openpyxl):
  /opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python extract_data.py ../../union/plasticome.v2.xlsx ../../../runs/2026-08-06.final-usearch.3/fetched_sequences.json

Facts in the page text follow petadex-manuscript/petadex-manuscript.tex.
"""
import json, sys, collections as C, pathlib, warnings
import openpyxl

warnings.filterwarnings("ignore")
xlsx = sys.argv[1]
fetched = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None
wb = openpyxl.load_workbook(xlsx, data_only=True)

def sheet(name, hdr=0):
    R = list(wb[name].iter_rows())
    h = [str(c.value).strip() if c.value is not None else "" for c in R[hdr]]
    body = [r for r in R[hdr + 1:] if any(c.value not in (None, "") for c in r)]
    return h, body

def s(v):
    if v is None: return ""
    if isinstance(v, float) and v.is_integer(): return str(int(v))
    return str(v).strip()

def fill(c):
    return c.fill.fgColor.rgb if c.fill is not None and c.fill.fill_type else None

out = {}

# PAZy pull
h, R = sheet("pazy-260701_pull"); ix = {k: i for i, k in enumerate(h)}
out["pull"] = {"n": len(R), "status": C.Counter(s(r[ix["status"]].value) for r in R)}

# manual curation: colour codes on the sequence column
h, R = sheet("pazy-260701_manual-curation")
SEQ = {"FFB7E1CD": "from_accession", "FFF4C7C3": "no_sequence", "FFFFE599": "fragment"}
seq_flag, name_changed, acc_added, activity, curation = C.Counter(), 0, 0, C.Counter(), {}
for r in R:
    f = SEQ.get(fill(r[9]), "unmarked"); seq_flag[f] += 1
    name_changed += fill(r[0]) == "FFB7E1CD"; acc_added += fill(r[1]) == "FFB7E1CD"
    act = s(r[2].value); activity[act] += 1
    pid = s(r[4].value)
    if pid:
        curation[pid] = {"seq": f, "activity": act, "note": s(r[3].value), "seq_note": s(r[8].value),
                         "renamed": fill(r[0]) == "FFB7E1CD", "acc_added": fill(r[1]) == "FFB7E1CD"}
out["manual"] = {"n": len(R), "seq_flag": seq_flag, "name_changed": name_changed,
                 "acc_added": acc_added, "activity": activity,
                 "activity_comments": sum(1 for r in R if r[3].value),
                 "seq_comments": sum(1 for r in R if r[8].value)}
_, R = sheet("cleaned_pazy-260701"); out["cleaned"] = {"n": len(R), "activity": C.Counter(s(r[2].value) for r in R)}

# v1 -> v1.1
h1, R1 = sheet("plasticome.v1"); h2, R2 = sheet("plasticome.v1.1")
ret1 = {s(r[0].value): s(r[1].value) for r in R1}
_, E = sheet("erickson-changes")
out["v1"] = {"n_v1": len(R1), "n_v11": len(R2),
             "retrieved_changed": sum(1 for r in R2 if s(r[0].value) in ret1 and ret1[s(r[0].value)] != s(r[1].value)),
             "doi_filled": sum(1 for r in R2 if s(r[h2.index("doi")].value) not in ("", "-")),
             "kingdom_labels": C.Counter(s(r[h2.index("kingdom")].value) for r in R2),
             "erickson": [{"name": s(r[8].value), "new": len(s(r[20].value)), "old": len(s(r[21].value))} for r in E]}

# union, annotated: the demo table
h, R = sheet("plasticome-v1-v260701-union-ann", 1); ix = {k: i for i, k in enumerate(h)}
rows = []
for r in R:
    g = lambda k: s(r[ix[k]].value)
    acc, src, pid, name = g("accession"), g("source"), g("pazy_id"), g("enzyme_name")
    seq = g("aa_sequence")
    if src == "v1.1":
        if name.startswith("Enzyme ") or name == "jmPE13": seqsrc = "Manual (paper supplement)"
        elif fetched and acc in fetched["resolved_by"]: seqsrc = fetched["resolved_by"][acc]
        else: seqsrc = "Retrieved from accession"
    else:
        seqsrc = "Retrieved from accession (PAZy pull)"
    cur = curation.get(pid, {})
    rows.append({"id": int(g("plasticome_id")), "name": name, "acc": acc, "pazy": pid, "src": src,
                 "doi": g("doi"), "sub": g("substrate"), "org": g("organism"), "sk": g("superkingdom"),
                 "meta": g("metagenomic") == "yes", "pdb": g("pdb_accession"), "uni": g("uniprot_accession"),
                 "cath": g("catalytic-cath"), "dom": g("domain"), "cent": g("centroid") == "yes",
                 "cid": g("centroid_id"), "len": len(seq), "seqsrc": seqsrc,
                 "act": cur.get("activity", ""), "note": cur.get("note", ""), "seqnote": cur.get("seq_note", ""),
                 "renamed": cur.get("renamed", False), "accfix": cur.get("acc_added", False)})
rows.sort(key=lambda x: x["id"])
out["union"] = rows

_, R = sheet("plasticome-union.non-redundant", 1)
out["nr"] = {"n": len(R), "members": C.Counter(int(float(r[2].value)) for r in R)}
_, R = sheet("plasticome_centroids_final", 1)
out["centroids"] = {"n": len(R), "components": len({s(r[4].value) for r in R})}
if fetched: out["fetch_by_stage"] = fetched["by_stage"]

p = pathlib.Path(__file__).with_name("data.json")
p.write_text(json.dumps(out, separators=(",", ":")))
print(p, p.stat().st_size, "bytes;", len(rows), "union rows")
print(json.dumps({k: v for k, v in out.items() if k != "union"}, indent=1)[:3000])
