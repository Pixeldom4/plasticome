#!/usr/bin/env python3
"""
Step 2: one alignment ensemble and one ML tree per component.

For a component with 4 or more tips:
  1. Muscle5 stratified ensemble, 16 alignments (4 guide tree permutations by
     4 HMM perturbations).
  2. Column confidence (CC) of a column is the fraction of the 16 alignments
     that contain exactly that column. The alignment with the highest mean CC
     is the reference. Columns with more than 50 percent gaps or CC below 0.5
     are dropped.
  3. IQ-TREE on the filtered reference, ModelFinder over LG, WAG, JTT and
     Q.pfam, 1000 ultrafast bootstraps.
  4. Ensemble bootstrap: 100 alignments made by drawing columns with
     replacement from all 16 alignments pooled (same two filters), a FastTree
     LG+gamma tree on each, and for every branch of the ML tree the share of
     the 100 trees that contain it.

UFBoot resamples columns of one fixed alignment, so it cannot see alignment
error. The ensemble bootstrap draws from alignments that disagree, so a branch
that depends on one way of aligning the family scores low. Both are written to
the tree as the internal node label "UFBoot/ensemble".

Muscle5 has a -resample command that does step 4. It is redone here in numpy
because it takes no seed, and the figure should come out the same twice.

3 tips: one alignment, IQ-TREE LG+G4, no support (there is only one topology).
2 tips: one alignment, FastTree LG+gamma distance.
1 tip:  nothing to infer.

Usage: python 02_component_trees.py [--only 0 3 4] [--force]
"""

import argparse
import collections
import os
import re
import shutil
from concurrent.futures import ThreadPoolExecutor
from io import StringIO

import numpy as np
from Bio import Phylo

from common import (IQ_THREADS, THREADS, WORK, WORK_ALL, cdir, read_fasta, read_tsv, run,
                    splits, tool, write_fasta, write_tree, write_tsv)

ENSEMBLE_N = 4          # per guide tree permutation, so 16 alignments
BOOT_N = 100
MAX_GAP = 0.5
MIN_CC = 0.5
SEED = 1


def read_efa(path):
    """Ensemble FASTA: '<name' starts a replicate, then ordinary FASTA."""
    reps, cur, name = collections.OrderedDict(), None, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("<"):
                cur = reps.setdefault(line[1:], collections.OrderedDict())
            elif line.startswith(">"):
                name = line[1:].split()[0]
                cur[name] = []
            elif line:
                cur[name].append(line)
    return {r: {k: "".join(v) for k, v in m.items()} for r, m in reps.items()}


def columns(msa, order):
    """Character matrix and one signature per column.

    A signature is the residue index of every sequence in that column (0 for a
    gap), so two columns from different alignments are the same column exactly
    when they align the same residues.
    """
    arr = np.array([list(msa[t]) for t in order])
    nongap = arr != "-"
    idx = np.where(nongap, np.cumsum(nongap, axis=1), 0).astype(np.int32)
    sigs = [idx[:, j].tobytes() for j in range(arr.shape[1])]
    return arr, nongap.mean(axis=0), sigs


def big_component(c, d, order, log):
    efa = os.path.join(d, "ensemble.efa")
    if not os.path.exists(efa):
        # Muscle writes the ensemble one alignment at a time. Write to a
        # scratch name and rename at the end, so an interrupted run never
        # leaves a short ensemble that the next run would take as complete.
        run([tool("muscle"), "-align", os.path.join(d, "seqs.faa"), "-stratified",
             "-replicates", str(ENSEMBLE_N), "-output", efa + ".part",
             "-threads", str(THREADS)], log)
        os.replace(efa + ".part", efa)
    reps = read_efa(efa)

    mats = {r: columns(m, order) for r, m in reps.items()}
    seen = collections.Counter()
    for arr, gap, sigs in mats.values():
        seen.update(set(sigs))
    cc = {r: np.array([seen[s] / len(reps) for s in sigs])
          for r, (arr, gap, sigs) in mats.items()}
    ok = {r: (mats[r][1] >= 1 - MAX_GAP) & (cc[r] >= MIN_CC) for r in mats}

    # reference: highest mean column confidence, ties by replicate name
    ref = sorted(reps, key=lambda r: (-cc[r].mean(), r))[0]
    arr = mats[ref][0]
    write_fasta(os.path.join(d, "reference.afa"),
                {t: "".join(arr[i]) for i, t in enumerate(order)})
    ml = os.path.join(d, "ml.afa")
    write_fasta(ml, {t: "".join(arr[i, ok[ref]]) for i, t in enumerate(order)})

    pre = os.path.join(d, "iq")
    # An interrupted search leaves a checkpoint and no report. Resume it then;
    # otherwise start clean so a stale checkpoint is never reused.
    resume = os.path.exists(pre + ".ckp.gz") and not os.path.exists(pre + ".iqtree")
    run([tool("iqtree"), "-s", ml, "-m", "MFP", "-mset", "LG,WAG,JTT,Q.pfam",
         "-B", "1000", "-seed", str(SEED), "-T", str(IQ_THREADS if len(order) >= 50 else 1),
         "-pre", pre, "-quiet"] + ([] if resume else ["-redo"]), log)
    model = "NA"
    with open(pre + ".iqtree") as fh:
        m = re.search(r"Best-fit model according to BIC: (\S+)", fh.read())
        if m:
            model = m.group(1)

    # ensemble bootstrap
    pool = np.concatenate([mats[r][0][:, ok[r]] for r in sorted(reps)], axis=1)
    n_cols = int(round(np.mean([ok[r].sum() for r in reps])))
    rng = np.random.default_rng(SEED)
    bdir = os.path.join(d, "boot")
    shutil.rmtree(bdir, ignore_errors=True)
    os.makedirs(bdir)
    jobs = []
    for b in range(BOOT_N):
        pick = rng.integers(0, pool.shape[1], n_cols)
        afa = os.path.join(bdir, f"{b:03d}.afa")
        write_fasta(afa, {t: "".join(pool[i, pick]) for i, t in enumerate(order)})
        jobs.append((b, afa))

    def fast(job):
        b, afa = job
        with open(afa[:-4] + ".nwk", "w") as out, open(os.devnull, "w") as err:
            import subprocess
            subprocess.run([tool("FastTree"), "-lg", "-gamma", "-quiet",
                            "-seed", str(SEED + b), afa],
                           stdout=out, stderr=err, check=True)
    with ThreadPoolExecutor(THREADS) as ex:
        list(ex.map(fast, jobs))

    tips = frozenset(order)
    counts = collections.Counter()
    for b, afa in jobs:
        counts.update(splits(Phylo.read(afa[:-4] + ".nwk", "newick"),
                             tips, order[0]).keys())
    with open(os.path.join(d, "boot.trees"), "w") as fh:
        for b, afa in jobs:
            fh.write(open(afa[:-4] + ".nwk").read())
    shutil.rmtree(bdir)

    tree = Phylo.read(pre + ".treefile", "newick")
    ufb, ens = [], []
    for sp, cl in splits(tree, tips, order[0]).items():
        e = round(100 * counts[sp] / BOOT_N)
        u = int(round(cl.confidence)) if cl.confidence is not None else 0
        ufb.append(u)
        ens.append(e)
        cl.confidence = None
        cl.name = f"{u}/{e}"
    write_tree(tree, os.path.join(d, "tree.nwk"))

    return [len(reps), arr.shape[1], int(ok[ref].sum()), f"{cc[ref].mean():.3f}",
            model, int(np.median(ufb)) if ufb else "NA",
            int(np.median(ens)) if ens else "NA",
            sum(1 for e in ens if e >= 70), len(ens)]


def small_component(c, d, order, log):
    afa = os.path.join(d, "reference.afa")
    run([tool("muscle"), "-align", os.path.join(d, "seqs.faa"), "-output", afa], log)
    msa = read_fasta(afa)
    arr = np.array([list(msa[t]) for t in order])
    shutil.copy(afa, os.path.join(d, "ml.afa"))
    if len(order) == 3:
        pre = os.path.join(d, "iq")
        run([tool("iqtree"), "-s", afa, "-m", "LG+G4", "-seed", str(SEED), "-T", "1",
             "-pre", pre, "-redo", "-quiet"], log)
        tree, model = Phylo.read(pre + ".treefile", "newick"), "LG+G4"
    else:
        import subprocess
        r = subprocess.run([tool("FastTree"), "-lg", "-gamma", "-quiet", afa],
                           capture_output=True, text=True, check=True)
        tree, model = Phylo.read(StringIO(r.stdout), "newick"), "LG+G (FastTree)"
    for cl in tree.find_clades(terminal=False):
        cl.confidence = None
    write_tree(tree, os.path.join(d, "tree.nwk"))
    return [1, arr.shape[1], arr.shape[1], "NA", model, "NA", "NA", 0, 0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=int, nargs="*")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    table = os.path.join(WORK, "component_trees.tsv")
    header = ["component", "n_tips", "alignments", "columns", "columns_kept",
              "mean_cc", "model", "median_ufboot", "median_ensemble",
              "branches_ensemble_ge70", "branches"]

    def on_disk():
        if not os.path.exists(table):
            return {}
        return {int(r["component"]): [r[h] for h in header] for r in read_tsv(table)}

    def save(mine):
        # merge with what is on disk, so two runs on different components do
        # not drop each other's rows
        rows = on_disk()
        rows.update(mine)
        write_tsv(table, header, [rows[k] for k in sorted(rows)])

    mine = {}
    for comp in read_tsv(os.path.join(WORK, "components.tsv")):
        c, n = int(comp["component"]), int(comp["n_tips"])
        if args.only and c not in args.only:
            continue
        d = cdir(c)
        order = sorted(read_fasta(os.path.join(d, "seqs.faa")), key=int)
        if n < 2:
            mine[c] = [c, n, 0, "NA", "NA", "NA", "NA", "NA", "NA", 0, 0]
            continue
        if (os.path.exists(os.path.join(d, "tree.nwk")) and c in on_disk()
                and not args.force):
            continue
        # Same sequences as in the all sequences run: the result is the same
        # too, since every step is seeded. Copy it instead of redoing it.
        twin = os.path.join(WORK_ALL, "components", f"c{c:02d}")
        twin_row = None
        if WORK != WORK_ALL and os.path.exists(os.path.join(twin, "tree.nwk")):
            same = (open(os.path.join(twin, "seqs.faa")).read()
                    == open(os.path.join(d, "seqs.faa")).read())
            done = {int(r["component"]): [r[h] for h in header] for r in
                    read_tsv(os.path.join(WORK_ALL, "component_trees.tsv"))}
            if same and c in done:
                twin_row = done[c]
        if twin_row and not args.force:
            shutil.copytree(twin, d, dirs_exist_ok=True)
            mine[c] = twin_row
            print(f"c{c:02d} n={n:<3} same sequences as the all sequences run, copied",
                  flush=True)
            save(mine)
            continue
        log = os.path.join(d, "log.txt")
        open(log, "a").close()
        stats = (big_component if n >= 4 else small_component)(c, d, order, log)
        mine[c] = [c, n] + stats
        stale = os.path.join(d, "PROVISIONAL")
        if os.path.exists(stale):
            os.remove(stale)
        print(f"c{c:02d} n={n:<3} " + " ".join(f"{h}={v}" for h, v in
                                              zip(header[2:], stats)), flush=True)
        save(mine)
    save(mine)


if __name__ == "__main__":
    main()
