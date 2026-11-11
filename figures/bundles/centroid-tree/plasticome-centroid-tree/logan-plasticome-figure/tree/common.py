"""Shared paths and helpers for the plasticome grafted tree."""

import csv
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SSN = os.path.join(HERE, "..", "ssn")
CSV = os.path.join(HERE, "..", "plasticome.v2 - plasticome-v1-v260701-union-annotated (1).csv")
NODES = os.path.join(SSN, "out", "ssn_nodes_c30_d40_k3.tsv")
HITS = os.path.join(SSN, "work", "allvall.tsv")

# Which sequences the trees are built from. Each set has its own work and out
# directory, so the runs sit side by side and neither overwrites the other.
#   all        every record of the union table (work/, out/)
#   centroids  only the 90 percent centroids, centroid = yes
#              (work_centroids/, out_centroids/)
# Components are the SSN components of the full set in both cases, so a
# component means the same thing here as in panel a-i.
SET = os.environ.get("PLASTICOME_TREE_SET", "all")
if SET not in ("all", "centroids"):
    raise SystemExit(f"PLASTICOME_TREE_SET must be all or centroids, not {SET}")
SUFFIX = "" if SET == "all" else f"_{SET}"
WORK = os.path.join(HERE, "work" + SUFFIX)
OUT = os.path.join(HERE, "out" + SUFFIX)
WORK_ALL = os.path.join(HERE, "work")

# Binaries (muscle, iqtree, FastTree, foldseek, foldmason, hhsuite) live in
# their own conda env so the plasticome env is left alone.
BIN = os.environ.get("PLASTICOME_TREE_BIN",
                     "/opt/homebrew/Caskroom/miniconda/base/envs/plasticome-tree/bin")
THREADS = int(os.environ.get("PLASTICOME_TREE_THREADS", os.cpu_count() or 4))
# IQ-TREE stops gaining from threads early on alignments this short (a few
# hundred columns), so it gets its own, smaller count. Muscle and the
# FastTree bootstraps use THREADS.
IQ_THREADS = int(os.environ.get("PLASTICOME_TREE_IQ_THREADS", min(THREADS, 8)))

# The SSN script owns the substrate, family and shape encodings and the
# palette, so the tree and the network cannot drift apart.
sys.path.insert(0, SSN)
from build_ssn import (SUBSTRATE_COLOR, SUBSTRATE_ORDER, SHAPE_LABEL,  # noqa: E402,F401
                       INK, MUTED, SURFACE, NODE_STROKE,
                       substrate_of, family_of, shape_of, seed_nodes, slug)


def tool(name):
    """A binary from BIN if it is there, otherwise from PATH, so the same
    scripts run inside any activated environment."""
    path = os.path.join(BIN, name)
    if os.path.exists(path):
        return path
    found = shutil.which(name)
    if not found:
        raise SystemExit(f"{name} not found in {BIN} or on PATH. Activate the "
                         "environment or set PLASTICOME_TREE_BIN.")
    return found


def run(cmd, log=None, cwd=None):
    """Run a command, fail loudly, keep stdout and stderr in `log`."""
    if log:
        with open(log, "a") as fh:
            fh.write("$ " + " ".join(cmd) + "\n")
            fh.flush()
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=cwd)
    else:
        r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        raise SystemExit(f"failed ({r.returncode}): {' '.join(cmd)}"
                         + (f"\n  see {log}" if log else ""))


def in_set(row):
    """Is this record part of the set being built?"""
    return SET == "all" or (row["centroid"] or "").strip().lower() == "yes"


def load_rows():
    fh = open(CSV, newline="")
    next(fh)  # provenance line above the header
    rows = {r["plasticome_id"]: r for r in csv.DictReader(fh)}
    fh.close()
    return rows


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")


def cdir(c):
    return os.path.join(WORK, "components", f"c{int(c):02d}")


def splits(tree, all_tips, ref):
    """Non-trivial bipartitions of a tree, each named by the side without
    `ref`, mapped to the clade below that branch."""
    out = {}
    for cl in tree.find_clades(terminal=False):
        tips = frozenset(t.name for t in cl.get_terminals())
        if ref in tips:
            tips = all_tips - tips
        if 1 < len(tips) < len(all_tips) - 1:
            out[tips] = cl
    return out


def write_tree(tree, path):
    from Bio import Phylo
    with open(path, "w") as fh:
        Phylo.write(tree, fh, "newick", format_branch_length="%1.8f")


def final_groups():
    """Component -> fold group after the homology check.

    Step 1 proposes fold groups from the CATH annotation. A component that
    step 4 could not confirm as a member is taken out and stands alone under a
    new group id, so nothing is grafted on annotation only.
    """
    comps = read_tsv(os.path.join(WORK, "components.tsv"))
    group = {int(r["component"]): int(r["fold_group"]) for r in comps}
    verdict = {int(r["component"]): r["verdict"] for r in
               read_tsv(os.path.join(WORK, "homology", "groups.tsv"))}
    nxt = max(group.values()) + 1
    for c in sorted(group):
        if verdict[c] == "unconfirmed":
            group[c] = nxt
            nxt += 1
    return group, verdict
