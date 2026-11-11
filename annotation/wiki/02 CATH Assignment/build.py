"""Inline the page data into template.html -> cath-assignment.html.

  /*DATA*/  data.json from extract_data.py
  /*TREE*/  fold groups of the forest figure (Supplementary Fig. 1), from
            figures/tree/runs/centroids/results/fold_table.tsv
"""
import csv, json, pathlib

here = pathlib.Path(__file__).parent
T = here.parents[2] / "figures" / "tree" / "runs" / "centroids" / "results"
tree = {"groups": [{"g": r["fold_group"], "fold": r["fold"], "cath": r["cath"], "n": int(r["sequences"])}
                   for r in csv.DictReader(open(T / "fold_table.tsv", newline=""), delimiter="\t")]}
html = ((here / "template.html").read_text()
        .replace("/*DATA*/null", (here / "data.json").read_text())
        .replace("/*TREE*/null", json.dumps(tree, separators=(",", ":"))))
(here / "cath-assignment.html").write_text(html)
print(here / "cath-assignment.html")
