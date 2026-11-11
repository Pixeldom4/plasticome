#!/usr/bin/env python3
"""Check that the centroid run finished and looks sane. Exits non-zero if not."""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(HERE, "logan-plasticome-figure", "tree")
WORK, OUT = os.path.join(TREE, "work_centroids"), os.path.join(TREE, "out_centroids")
bad = []


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check(ok, msg):
    print(("ok    " if ok else "FAIL  ") + msg)
    if not ok:
        bad.append(msg)


need = ["figure_pie.svg", "figure_pie.pdf", "figure_pie.png", "figure_pie_table.svg",
        "figure_tree.svg", "forest.nwk", "tree_tips.tsv", "fold_table.tsv",
        "backbone_branches.tsv", "backbone_summary.tsv", "homology_check.tsv",
        "stems.tsv", "cath_domains.xlsx"]
for f in need:
    check(os.path.exists(os.path.join(OUT, f)), f"out_centroids/{f} exists")

try:
    comps = {int(r["component"]): r for r in tsv(os.path.join(WORK, "component_trees.tsv"))}
    check(len(comps) == 56, f"56 components in component_trees.tsv (found {len(comps)})")
    c0 = comps.get(0, {})
    check(c0.get("n_tips") == "225", f"component 0 has 225 tips (found {c0.get('n_tips')})")
    check(c0.get("alignments") == "16", f"component 0 used 16 alignments (found {c0.get('alignments')})")
    check(c0.get("model") not in (None, "", "NA"), f"component 0 has a model ({c0.get('model')})")
    c0dir = os.path.join(WORK, "components", "c00")
    for f in ("tree.nwk", "iq.iqtree", "boot.trees", "ensemble.efa", "reference.afa"):
        check(os.path.exists(os.path.join(c0dir, f)), f"work_centroids/components/c00/{f} exists")
    if os.path.exists(os.path.join(c0dir, "boot.trees")):
        n = sum(1 for _ in open(os.path.join(c0dir, "boot.trees")))
        check(n == 100, f"100 ensemble bootstrap trees for component 0 (found {n})")
    tips = tsv(os.path.join(OUT, "tree_tips.tsv"))
    check(len(tips) == 411, f"411 tips in tree_tips.tsv (found {len(tips)})")
    check(len({t["tip"] for t in tips}) == len(tips), "no tip appears twice")
    n_trees = sum(1 for _ in open(os.path.join(OUT, "forest.nwk")))
    folds = len({t["fold_group"] for t in tips})
    check(n_trees == folds, f"one tree per fold group ({n_trees} trees, {folds} fold groups)")
    print(f"\ncomponent 0: columns kept {c0.get('columns_kept')} of {c0.get('columns')}, "
          f"median UFBoot {c0.get('median_ufboot')}, median ensemble {c0.get('median_ensemble')}, "
          f"{c0.get('branches_ensemble_ge70')} of {c0.get('branches')} branches at ensemble >= 70")
    for r in tsv(os.path.join(OUT, "backbone_summary.tsv")):
        print(f"backbone {r['fold_group']}: {r['tips']} tips, {r['ufboot_ge95_and_cross_checked']} "
              f"of {r['branches']} branches kept")
except FileNotFoundError as e:
    check(False, f"missing file: {e.filename}")

print("\n" + ("ALL CHECKS PASSED" if not bad else f"{len(bad)} CHECK(S) FAILED"))
sys.exit(1 if bad else 0)
