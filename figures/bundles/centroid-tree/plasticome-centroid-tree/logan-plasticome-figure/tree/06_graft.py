#!/usr/bin/env python3
"""
Step 6: graft the component trees onto their backbone.

One rooted tree per fold group. A fold group with a single component is just
that component's tree. A fold group with several is built like this:

  1. Root the backbone at its midpoint.
  2. Collapse every backbone branch that is not trustworthy into a polytomy.
     A branch is kept when its UFBoot is at least --min-boot (default 95)
     and at least one of the two cross-check trees from step 5 recovers it.
  3. Root each component tree at its midpoint.
  4. Replace each backbone tip with the rooted tree of its component.

The stem that carries a component is the backbone's terminal branch for that
component's representative. That branch runs to the representative itself, and
the component tree then runs from its own root back out to the representative,
so the root to representative distance is taken off the stem. Otherwise that
stretch would be counted twice. The stem is never shortened below
--min-stem (default 0.02 substitutions per site).

Internal node names carry what kind of branch sits above the node and its
support, so nothing downstream has to guess:

  B_<ufboot>_<aa>_<tm>      backbone branch; aa and tm are 1 when the amino
                            acid only tree and the TM-score tree recover it
  S_c<NN>                   stem of component NN (root of its subtree)
  C_<ufboot>_<ensemble>     branch inside a component
  X                         no support defined (root, or a 2 or 3 tip tree)

Midpoint rooting is a convention, not a result. Neither the backbone root nor
the root of a component tree is supported by an outgroup.

Writes out/trees/gNN.nwk, out/forest.nwk and out/tree_tips.tsv.
"""

import argparse
import collections
import os
from io import StringIO

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from common import (OUT, WORK, cdir, final_groups, read_tsv, splits, write_tree,
                    write_tsv)


def label_by_split(tree, ref):
    """Internal labels of an unrooted tree, keyed by the split they describe."""
    tips = frozenset(t.name for t in tree.get_terminals())
    return {sp: cl.name for sp, cl in splits(tree, tips, ref).items() if cl.name}


def root_midpoint(tree, prefix):
    """Midpoint root, then put every support label back on the branch it
    belongs to. Rerooting moves nodes relative to their labels, so the labels
    are looked up again by the tips on each side of the branch."""
    tips = frozenset(t.name for t in tree.get_terminals())
    if len(tips) < 3:
        for cl in tree.find_clades(terminal=False):
            cl.name = "X"
        return tree
    ref = sorted(tips, key=lambda t: (len(t), t))[0]
    labels = label_by_split(tree, ref)
    tree.root_at_midpoint()
    unary = True
    while unary:
        unary = False
        for cl in list(tree.find_clades(terminal=False)):
            for i, ch in enumerate(cl.clades):
                if not ch.is_terminal() and len(ch.clades) == 1:
                    only = ch.clades[0]
                    only.branch_length = (only.branch_length or 0) + (ch.branch_length or 0)
                    cl.clades[i] = only
                    unary = True
    for cl in tree.find_clades(terminal=False):
        below = frozenset(t.name for t in cl.get_terminals())
        sp = tips - below if ref in below else below
        lab = labels.get(sp) if 1 < len(sp) < len(tips) - 1 else None
        cl.name = f"{prefix}_{lab.replace('/', '_')}" if lab else "X"
    tree.root.name = "X"
    for cl in tree.find_clades():
        cl.confidence = None
        if cl.branch_length is not None and cl.branch_length < 0:
            cl.branch_length = 0.0
    tree.root.branch_length = None
    return tree


def component_tree(c):
    """Rooted tree of one component, or a single tip."""
    path = os.path.join(cdir(c), "tree.nwk")
    if os.path.exists(path):
        return root_midpoint(Phylo.read(path, "newick"), "C")
    tip = next(iter(open(os.path.join(cdir(c), "seqs.faa")))).strip()[1:]
    return Tree(root=Clade(name=tip))


def keep_backbone(name, min_boot, cross_check=True):
    _, boot, aa, tm = name.split("_")
    return int(boot) >= min_boot and (not cross_check or aa == "1" or tm == "1")


def collapse_backbone(tree, min_boot, cross_check):
    """Turn untrusted backbone branches into polytomies."""
    kept = dropped = 0
    changed = True
    while changed:
        changed = False
        for parent in tree.find_clades(terminal=False):
            for ch in list(parent.clades):
                if ch.is_terminal() or not ch.name.startswith("B_"):
                    continue
                if keep_backbone(ch.name, min_boot, cross_check):
                    continue
                i = parent.clades.index(ch)
                for g in ch.clades:
                    g.branch_length = (g.branch_length or 0) + (ch.branch_length or 0)
                parent.clades[i:i + 1] = ch.clades
                dropped += 1
                changed = True
    # an unlabelled internal node below the root is one half of the root
    # branch whose other half was collapsed; it carries no claim of its own
    for parent in tree.find_clades(terminal=False):
        for ch in list(parent.clades):
            if not ch.is_terminal() and ch.name == "X" and parent is tree.root:
                i = parent.clades.index(ch)
                for g in ch.clades:
                    g.branch_length = (g.branch_length or 0) + (ch.branch_length or 0)
                parent.clades[i:i + 1] = ch.clades
    kept = sum(1 for cl in tree.find_clades(terminal=False) if cl.name.startswith("B_"))
    return kept, dropped


def graft(g, comps, reps, min_boot, min_stem, cross_check):
    if len(comps) == 1:
        return component_tree(comps[0]), None
    bb = Phylo.read(os.path.join(WORK, "backbone", f"g{g:02d}", "tree.nwk"), "newick")
    n_branches = sum(1 for cl in bb.find_clades(terminal=False)
                     if cl.name and "/" in cl.name)
    bb = root_midpoint(bb, "B")
    # the two root children describe one unrooted branch; count it once
    kept, dropped = collapse_backbone(bb, min_boot, cross_check)
    root_pair = [ch for ch in bb.root.clades if ch.name.startswith("B_")]
    kept_unique = kept - (1 if len(bb.root.clades) == 2 and len(root_pair) == 2 else 0)

    stems = []
    for tip in bb.get_terminals():
        c = int(tip.name[1:])
        sub = component_tree(c)
        if sub.root.is_terminal():
            tip.name = sub.root.name
            continue
        d_rep = sub.distance(sub.root, reps[c])
        stem = max((tip.branch_length or 0) - d_rep, min_stem)
        stems.append([f"c{c:02d}", f"{tip.branch_length:.4f}", f"{d_rep:.4f}",
                      f"{stem:.4f}", int((tip.branch_length or 0) - d_rep < min_stem)])
        tip.name = f"S_c{c:02d}"
        tip.branch_length = stem
        tip.clades = sub.root.clades
    return bb, (n_branches, kept_unique, stems)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-boot", type=int, default=95)
    ap.add_argument("--min-stem", type=float, default=0.02)
    ap.add_argument("--no-cross-check", action="store_true",
                    help="keep backbone branches on UFBoot alone")
    args = ap.parse_args()

    if os.path.exists(os.path.join(cdir(0), "PROVISIONAL")):
        print("WARNING: component 0 is using a provisional tree")

    comps = read_tsv(os.path.join(WORK, "components.tsv"))
    groups = collections.defaultdict(list)
    final, verdict = final_groups()
    for c, g in final.items():
        groups[g].append(c)
    info = {int(r["component"]): r for r in comps}
    reps = {int(r["component"]): r["representative"] for r in
            read_tsv(os.path.join(WORK, "structures.tsv"))}
    tips = {r["tip"]: r for r in read_tsv(os.path.join(WORK, "tips.tsv"))}

    os.makedirs(os.path.join(OUT, "trees"), exist_ok=True)
    rows, forest, stem_rows = [], [], []
    for g in sorted(groups):
        tree, bb = graft(g, sorted(groups[g]), reps, args.min_boot, args.min_stem,
                         not args.no_cross_check)
        tree.ladderize()
        path = os.path.join(OUT, "trees", f"g{g:02d}.nwk")
        write_tree(tree, path)
        forest.append(open(path).read().strip())
        names = [t.name for t in tree.get_terminals()]
        expected = sorted(t for t, r in tips.items() if r["status"] == "kept"
                          and final[int(r["component"])] == g)
        if sorted(names) != expected:
            raise SystemExit(f"g{g:02d}: grafted tree does not hold exactly the "
                             "kept tips of its fold group")
        for i, t in enumerate(names):
            r = tips[t]
            c = int(r["component"])
            rows.append([t, f"g{g:02d}", info[c]["fold_key"], f"c{c:02d}", i,
                         int(reps[c] == t), r["enzyme_name"], r["accession"],
                         r["records"], r["substrate_group"], r["shape"],
                         r["family_group"], verdict[c]])
        fam = info[sorted(groups[g])[0]]["family"]
        if bb:
            n, kept, stems = bb
            stem_rows += [[f"g{g:02d}"] + s for s in stems]
            print(f"g{g:02d} {fam}: {len(names)} tips from {len(groups[g])} components; "
                  f"{kept} of {n} backbone branches kept, {n - kept} collapsed; "
                  f"{sum(s[4] for s in stems)} stems at the minimum length")
        else:
            print(f"g{g:02d} {fam}: {len(names)} tips, one component")

    with open(os.path.join(OUT, "forest.nwk"), "w") as fh:
        fh.write("\n".join(forest) + "\n")
    write_tsv(os.path.join(OUT, "tree_tips.tsv"),
              ["tip", "fold_group", "fold_key", "component", "order_in_tree",
               "is_representative", "enzyme_name", "accession", "records",
               "substrate_group", "shape", "family_group", "homology"], rows)
    write_tsv(os.path.join(OUT, "stems.tsv"),
              ["fold_group", "component", "backbone_terminal_branch",
               "root_to_representative", "stem", "at_minimum"], stem_rows)
    left_out = [r for r in tips.values() if r["status"] != "kept"]
    print(f"{len(rows)} tips in {len(groups)} trees; "
          f"{len(left_out)} sequences left out as fragments")


if __name__ == "__main__":
    main()
