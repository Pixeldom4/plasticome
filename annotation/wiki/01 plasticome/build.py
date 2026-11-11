"""Inline the page data into template.html -> plasticome-curation.html.

  /*DATA*/  data.json from extract_data.py (curation workbook)
  /*CATH*/  per-enzyme CATH calls and domains, from ../02 CATH Assignment/data.json
  /*TREE*/  fold groups of the forest figure (Supplementary Fig. 1), from
            figures/tree/runs/centroids/results/
"""
import csv, json, pathlib

here = pathlib.Path(__file__).parent
repo = here.parents[2]          # annotation/wiki/<page>/ -> repo root
cath = json.loads((here.parent / "02 CATH Assignment" / "data.json").read_text())
cath = {"enz": {e["id"]: e for e in cath["enz"]}, "dom": cath["dom"], "names": cath["names"]}

T = repo / "figures" / "tree" / "runs" / "centroids" / "results"
rd = lambda p: list(csv.DictReader(open(p, newline=""), delimiter="\t"))
tree = {"groups": [{"g": r["fold_group"], "fold": r["fold"], "cath": r["cath"], "n": int(r["sequences"]),
                    "comp": int(r["components"]), "homology": r["homology"]} for r in rd(T / "fold_table.tsv")],
        "tip": {r["tip"]: r["fold_group"] for r in rd(T / "tree_tips.tsv")}}

dump = lambda o: json.dumps(o, separators=(",", ":"))
html = ((here / "template.html").read_text()
        .replace("/*DATA*/null", (here / "data.json").read_text())
        .replace("/*CATH*/null", dump(cath))
        .replace("/*TREE*/null", dump(tree)))
(here / "plasticome-curation.html").write_text(html)
print(here / "plasticome-curation.html")
