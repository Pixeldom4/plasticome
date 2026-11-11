#!/usr/bin/env bash
# Finish the 90 percent centroid tree. Run inside the plasticome-tree
# environment, on a compute node, from anywhere.
#
#   PLASTICOME_TREE_THREADS      cores for Muscle and the bootstrap trees
#                                (default: all cores the node reports)
#   PLASTICOME_TREE_IQ_THREADS   cores for IQ-TREE (default: min(threads, 8))
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/logan-plasticome-figure/tree"

export PLASTICOME_TREE_SET=centroids
export PLASTICOME_TREE_THREADS="${PLASTICOME_TREE_THREADS:-$(nproc)}"
export MPLBACKEND=Agg
# no binaries at the Mac path, so the scripts take them from PATH
export PLASTICOME_TREE_BIN="${PLASTICOME_TREE_BIN:-/nonexistent}"

step() { echo "== $1 ($(date '+%F %T'))"; python "$1.py" "${@:2}"; }

# Steps 1 and 3 are already done and shipped (work_centroids/). Step 3 needs
# the internet unless every lookup is in work/structures/cache, which it is.
step 02_component_trees     # only component 0 is left; the other 55 are skipped
step 04_homology
step 05_backbone
step 06_graft
step 07_draw
step 08_draw_pie
step 09_cath_table

python "$ROOT/verify.py"
