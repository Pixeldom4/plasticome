#!/usr/bin/env python3
"""
Step 5: the backbone that joins the components of one fold group.

One tip per component, the representative structure from step 3. Between
components sequence identity is around 10 percent on the alignable core, too
low for a sequence alignment to mean much, so the alignment comes from
structure and the tree uses structure as data too:

  1. FoldMason structural alignment of the representatives. It writes the
     same alignment twice, as amino acids and as 3Di structural states.
  2. Columns present in fewer than half the representatives are dropped.
  3. IQ-TREE on both alphabets as two partitions with shared topology and
     proportional branch lengths. Amino acid model chosen by ModelFinder from
     LG, WAG, JTT and Q.pfam; 3Di under the 3DI matrix (models/3DI).
     1000 ultrafast bootstraps with the NNI correction (-bnni), which is the
     variant meant for data where the model fits poorly. Standard bootstrap
     was tried first and took over a minute per replicate here.
  4. Branch lengths are then re-estimated on the amino acid partition alone,
     topology fixed, so the backbone is in amino acid substitutions per site
     like the component trees it will carry.

A deep branch can have a good bootstrap and still be an artefact of the
method, which is the criticism made of the RdRp global tree. So every backbone
branch is also checked against two trees that use the data differently:

  aa_only   ML on the amino acid partition alone
  tm_nj     neighbour joining on 1 - TM-score from step 4, no alignment
            columns and no substitution model

Step 6 keeps a backbone branch only when its UFBoot is 95 or more (the
threshold UFBoot is calibrated for, not the 70 used for standard bootstrap)
and at least one of the two also recovers it. The rest become polytomies.

A component with no structure is added to the amino acid alignment by sequence
(mafft --add --keeplength) and is missing data in the 3Di partition.

Groups of 3 have one unrooted topology, so no support is computed. Groups of
2 have one branch, and its length is the FastTree LG+gamma distance.

Writes work/backbone/gNN/tree.nwk and out/backbone_branches.tsv.
"""

import collections
import os
import re
import shutil
import subprocess
from io import StringIO

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from common import (HERE, OUT, THREADS, WORK, cdir, final_groups, read_fasta,
                    read_tsv, run, splits, tool, write_fasta, write_tsv, write_tree)

BDIR = os.path.join(WORK, "backbone")
MIN_OCCUPANCY = 0.5
UFBOOT_N = 1000
UFBOOT_KEEP = 95
SEED = 1


def best_model(iqtree_report):
    with open(iqtree_report) as fh:
        m = re.search(r"Best-fit model according to BIC: (\S+)", fh.read())
    return m.group(1) if m else "LG+G4"


def tm_tree(names, tm):
    """Neighbour joining on 1 - TM, tips named like the backbone tips."""
    order = sorted(names)
    rows = [[max(0.0, 1.0 - tm.get((min(a, b), max(a, b)), 0.0)) if a != b else 0.0
             for b in order[:i + 1]] for i, a in enumerate(order)]
    dm = DistanceMatrix([f"c{c:02d}" for c in order], rows)
    return DistanceTreeConstructor().nj(dm)


def project(tree, keep, ref):
    """Splits of `tree` restricted to the tips in `keep`."""
    out = set()
    for cl in tree.find_clades(terminal=False):
        tips = frozenset(t.name for t in cl.get_terminals()) & keep
        if ref in tips:
            tips = keep - tips
        if 1 < len(tips) < len(keep) - 1:
            out.add(tips)
    return out


def build(g, comps, has_struct, tm, log):
    d = os.path.join(BDIR, f"g{g:02d}")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "structures"))
    st = sorted(c for c in comps if c in has_struct)
    for c in st:
        shutil.copy(os.path.join(WORK, "structures", f"c{c:02d}.pdb"),
                    os.path.join(d, "structures"))
    run([tool("foldmason"), "easy-msa"] +
        [os.path.join(d, "structures", f"c{c:02d}.pdb") for c in st] +
        [os.path.join(d, "fm"), os.path.join(d, "tmp"), "--threads", str(THREADS),
         "-v", "1"], log)
    shutil.rmtree(os.path.join(d, "tmp"), ignore_errors=True)
    aa = read_fasta(os.path.join(d, "fm_aa.fa"))
    di = read_fasta(os.path.join(d, "fm_3di.fa"))
    width = len(next(iter(aa.values())))

    # components without a structure join the amino acid alignment by sequence
    extra = sorted(c for c in comps if c not in has_struct)
    if extra:
        reps = {int(r["component"]): r["representative"] for r in
                read_tsv(os.path.join(WORK, "structures.tsv"))}
        add = {f"c{c:02d}": read_fasta(os.path.join(cdir(c), "seqs.faa"))[reps[c]]
               for c in extra}
        write_fasta(os.path.join(d, "add.faa"), add)
        r = subprocess.run([tool("mafft"), "--add", os.path.join(d, "add.faa"),
                            "--keeplength", "--quiet", os.path.join(d, "fm_aa.fa")],
                           capture_output=True, text=True, check=True)
        with open(os.path.join(d, "fm_aa_added.fa"), "w") as fh:
            fh.write(r.stdout)
        aa = {k: v.upper() for k, v in
              read_fasta(os.path.join(d, "fm_aa_added.fa")).items()}
        for k in add:
            di[k] = "-" * width

    order = sorted(aa)
    A = np.array([list(aa[t]) for t in order])
    D = np.array([list(di[t]) for t in order])
    st_rows = [i for i, t in enumerate(order) if int(t[1:]) in has_struct]
    keep = (A[st_rows] != "-").mean(axis=0) >= MIN_OCCUPANCY
    write_fasta(os.path.join(d, "aa.afa"), {t: "".join(A[i, keep]) for i, t in enumerate(order)})
    write_fasta(os.path.join(d, "3di.afa"), {t: "".join(D[i, keep]) for i, t in enumerate(order)})
    K = A[:, keep]
    ident = [np.mean(K[i][m] == K[j][m]) for i in range(len(order))
             for j in range(i) if (m := (K[i] != "-") & (K[j] != "-")).sum()]
    stats = {"tips": len(order), "columns": width, "columns_kept": int(keep.sum()),
             "median_identity": float(np.median(ident)), "added_by_sequence": len(extra)}

    tips = frozenset(order)
    if len(order) == 2:
        r = subprocess.run([tool("FastTree"), "-lg", "-gamma", "-quiet",
                            os.path.join(d, "aa.afa")],
                           capture_output=True, text=True, check=True)
        tree = Phylo.read(StringIO(r.stdout), "newick")
        for cl in tree.find_clades(terminal=False):
            cl.confidence = None
        write_tree(tree, os.path.join(d, "tree.nwk"))
        stats["model"] = "LG+G (FastTree)"
        return stats, []

    # amino acid model, and the amino-acid-only tree used as a cross-check
    boot = ["-B", str(UFBOOT_N), "-bnni"] if len(order) >= 4 else []
    run([tool("iqtree"), "-s", "aa.afa", "-m", "MFP", "-mset", "LG,WAG,JTT,Q.pfam",
         "-seed", str(SEED), "-T", "1", "-pre", "aa", "-redo", "-quiet"] + boot,
        log, cwd=d)
    model = best_model(os.path.join(d, "aa.iqtree"))
    stats["model"] = f"{model} | 3DI+G4"

    shutil.copy(os.path.join(HERE, "models", "3DI"), os.path.join(d, "3DI"))
    with open(os.path.join(d, "partitions.nex"), "w") as fh:
        fh.write("#nexus\nbegin sets;\n"
                 "  charset aa = aa.afa: AA, *;\n"
                 "  charset di = 3di.afa: AA, *;\n"
                 f"  charpartition both = {model}: aa, 3DI+G4: di;\n"
                 "end;\n")
    run([tool("iqtree"), "-p", "partitions.nex", "-seed", str(SEED),
         "-T", "1", "-pre", "both", "-redo", "-quiet"] + boot, log, cwd=d)

    # same topology, amino acid branch lengths
    run([tool("iqtree"), "-s", "aa.afa", "-m", model, "-te", "both.treefile",
         "-seed", str(SEED), "-T", "1", "-pre", "aa_lengths", "-redo", "-quiet"],
        log, cwd=d)

    both = Phylo.read(os.path.join(d, "both.treefile"), "newick")
    tree = Phylo.read(os.path.join(d, "aa_lengths.treefile"), "newick")
    if len(order) < 4:
        write_tree(tree, os.path.join(d, "tree.nwk"))
        return stats, []

    aa_tree = Phylo.read(os.path.join(d, "aa.treefile"), "newick")
    ref = order[0]
    sp_both = splits(both, tips, ref)
    sp_aa = splits(aa_tree, tips, ref)
    st_names = frozenset(t for t in order if int(t[1:]) in has_struct)
    nj = tm_tree([int(t[1:]) for t in st_names], tm)
    Phylo.write(nj, os.path.join(d, "tm_nj.nwk"), "newick")
    sref = sorted(st_names)[0]
    sp_nj = project(nj, st_names, sref)

    branches = []
    for sp, cl in splits(tree, tips, ref).items():
        b = int(round(sp_both[sp].confidence or 0))
        in_aa = sp in sp_aa
        aa_boot = int(round(sp_aa[sp].confidence or 0)) if in_aa else 0
        proj = sp & st_names
        if sref in proj:
            proj = st_names - proj
        in_nj = "NA" if not 1 < len(proj) < len(st_names) - 1 else int(proj in sp_nj)
        cl.confidence = None
        cl.name = f"{b}/{int(in_aa)}/{in_nj}"
        small = sp if len(sp) <= len(tips) - len(sp) else tips - sp
        branches.append([f"g{g:02d}", ";".join(sorted(small)), len(small), b,
                         int(in_aa), aa_boot, in_nj, f"{cl.branch_length:.4f}"])
    write_tree(tree, os.path.join(d, "tree.nwk"))
    return stats, branches


def main():
    os.makedirs(BDIR, exist_ok=True)
    groups = collections.defaultdict(list)
    final, verdict = final_groups()
    for c, g in final.items():
        groups[g].append(c)
    for c in sorted(c for c, v in verdict.items() if v == "unconfirmed"):
        print(f"c{c:02d} not confirmed as a member of its annotated fold group, "
              "left out of the backbone")
    has_struct = {int(r["component"]) for r in
                  read_tsv(os.path.join(WORK, "structures.tsv")) if r["source"] != "none"}
    tm = {(int(r["a"]), int(r["b"])): float(r["tm"]) for r in
          read_tsv(os.path.join(WORK, "homology", "pairs.tsv")) if r["tm"] != "NA"}
    all_branches, summary = [], []
    for g in sorted(groups):
        if len(groups[g]) < 2:
            continue
        log = os.path.join(BDIR, f"g{g:02d}.log")
        open(log, "w").close()
        stats, branches = build(g, groups[g], has_struct, tm, log)
        all_branches += branches
        ok = [b for b in branches if b[3] >= UFBOOT_KEEP]
        summary.append([f"g{g:02d}", stats["tips"], stats["added_by_sequence"],
                        stats["columns"], stats["columns_kept"],
                        f"{stats['median_identity']:.3f}", stats["model"],
                        len(branches), len(ok),
                        sum(1 for b in ok if b[4] == 1 or b[6] == 1)])
        print(f"g{g:02d}: {stats['tips']} tips, {stats['columns_kept']} of "
              f"{stats['columns']} columns, median identity "
              f"{stats['median_identity']:.3f}, {stats['model']}; "
              f"{len(ok)} of {len(branches)} branches at UFBoot >= {UFBOOT_KEEP}", flush=True)
    write_tsv(os.path.join(OUT, "backbone_branches.tsv"),
              ["fold_group", "smaller_side", "size", "ufboot", "in_aa_only_tree",
               "aa_only_ufboot", "in_tm_nj_tree", "length"], all_branches)
    write_tsv(os.path.join(OUT, "backbone_summary.tsv"),
              ["fold_group", "tips", "added_by_sequence", "columns", "columns_kept",
               "median_identity", "model", "branches", "ufboot_ge95",
               "ufboot_ge95_and_cross_checked"], summary)


if __name__ == "__main__":
    main()
