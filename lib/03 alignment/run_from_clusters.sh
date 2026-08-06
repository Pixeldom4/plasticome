#!/usr/bin/env bash
# End-to-end driver: cluster-centroids TSV -> components.
#
# Same method as run.sh / the central notebook -- permissive all-vs-all search,
# >=30% aaid AND e-value < 1e-5 applied as a POST-filter, single-linkage
# connected components -- with two differences: the node set comes from one row
# per cluster (the representative) instead of the curated PAZy TSV, and there is
# no Step-0 engine sanity check (run.sh keeps that; it is a fixed 213-node
# usearch fixture, independent of this input).
#
# The aligner is DIAMOND (`diamond blastp`, all-vs-all, settings mapped one to one
# from the usearch method -- see scripts/step2_align.py and config.py). It runs
# natively when a `diamond` is on PATH, else bin/diamond through Docker
# linux/amd64. ENGINE=usearch switches back to usearch11 -allpairs_local, which
# always needs Docker on this host.
#
# Usage:
#   ./run_from_clusters.sh "<clusters.tsv>" ["<outdir>"]
#   ./run_from_clusters.sh "<union.tsv>"    ["<outdir>"]
#
# The UNION table (step 1's "01-union.tsv") is accepted in place of the clusters
# table, which is the align-BEFORE-cluster arrangement: every union row becomes a
# node and clustering never runs. The shape is detected from the header by
# scripts/clusters_to_tsv.py and scripts/annotate_clusters.py; the aligner, the
# filter and the partition are identical either way, so the two runs are directly
# comparable. With a union input, UNION defaults to that same file.
#
# The run deliverable is "<dir of clusters.tsv>/03 alignment.tsv" -- the input
# table with one column added, `component_id` -- and every intermediate lands in
# "<...>/03 alignment.intermediates/", mirroring how the clustering step pairs
# "02-clusters.tsv" with "02-clusters.intermediates/". The sidecar takes the
# deliverable's stem with its extension REPLACED (the older scheme appended
# ".work" instead). Passing <outdir> moves only the intermediates. Both paths
# must live under the repo root, which is the single Docker mount.
#
# Env overrides:
#   PY=...        python 3.9+           (default: the plasticome conda env)
#   ENGINE=...    diamond | usearch     (default: diamond)
#   OUT_TSV=...   the deliverable path  (default: "<dir>/03 alignment.tsv")
#   TAG=...       output file stem      (default: clusters)
#   DATE=...      output file stem      (default: today)
#   ID_MODE=...   cluster | label       (default: cluster; see clusters_to_tsv.py)
#   V1=...        v1 overlay table      (default: source-data/plasticome.v1.1/plasticome.v1.1.tsv)
#   UNION=...     union TSV supplying sequences for v1 rows that carry none, so the
#                 overlay can still key on sequence md5 (optional; see step1_nodes.py)
set -euo pipefail
# Absolute + physical path; only the parent has to exist. `pwd -P` matters: the
# paths are later matched against $REPO to build the Docker mount-relative form.
abspath() {
  local d b
  d="$(dirname "$1")"; b="$(basename "$1")"
  d="$(cd "$d" && pwd -P)" && printf '%s/%s\n' "$d" "$b"
}

# Drop the extension from the BASENAME only. Run directories carry dots
# ("runs/2026-07-30.plasticome.v1.212.union-spec/"), so a bare "${p%.*}" would
# chew into the parent path instead of the filename.
stem() {
  local d b
  d="$(dirname "$1")"; b="$(basename "$1")"
  printf '%s/%s\n' "$d" "${b%.*}"
}

# Resolve the arguments against the CALLER's cwd, before moving to the script dir.
CLUSTERS="$(abspath "${1:?usage: $0 \"<clusters.tsv>\" [\"<outdir>\"]}")"
OUT_TSV="$(abspath "${OUT_TSV:-$(dirname "$CLUSTERS")/03 alignment.tsv}")"
OUT="$(abspath "${2:-$(stem "$OUT_TSV").intermediates}")"
cd "$(dirname "$0")"                       # -> lib/03 alignment/
REPO="$(cd ../.. && pwd -P)"               # -> repo root == the Docker /d mount

PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python}"
ENGINE="${ENGINE:-diamond}"
TAG="${TAG:-clusters}"
DATE="${DATE:-$(date +%F)}"
ID_MODE="${ID_MODE:-cluster}"
V1="${V1:-$REPO/source-data/plasticome.v1.1/plasticome.v1.1.tsv}"

[ -f "$CLUSTERS" ]           || { echo "no such clusters TSV: $CLUSTERS" >&2; exit 1; }
[ -x "$PY" ] || command -v "$PY" >/dev/null || { echo "no python at: $PY" >&2; exit 1; }
case "$ENGINE" in
  diamond)
    # A native diamond skips Docker entirely; the bundled copy is a Linux ELF and does not.
    if ! command -v diamond >/dev/null; then
      [ -f "$REPO/bin/diamond" ] || { echo "no diamond on PATH and missing $REPO/bin/diamond" >&2; exit 1; }
      docker info >/dev/null 2>&1 || { echo "Docker is not running (bin/diamond is a Linux ELF binary)" >&2; exit 1; }
    fi ;;
  usearch)
    [ -f "$REPO/bin/usearch11" ] || { echo "missing $REPO/bin/usearch11" >&2; exit 1; }
    docker info >/dev/null 2>&1 || { echo "Docker is not running (usearch is a Linux ELF binary)" >&2; exit 1; } ;;
  *) echo "ENGINE must be diamond or usearch (got: $ENGINE)" >&2; exit 1 ;;
esac
case "$CLUSTERS" in "$REPO"/*) ;; *) echo "clusters TSV must be under $REPO" >&2; exit 1;; esac
case "$OUT"      in "$REPO"/*) ;; *) echo "outdir must be under $REPO" >&2; exit 1;; esac
mkdir -p "$OUT"

# A union input IS the sequence source the v1 overlay wants to key on, so it
# supplies its own --v1-seqs unless the caller named another one.
if [ -z "${UNION:-}" ] && head -1 "$CLUSTERS" | grep -q $'^plasticome_id\t'; then
  UNION="$CLUSTERS"
fi

echo "== Step A: input TSV -> curated node shape =="
"$PY" scripts/clusters_to_tsv.py --clusters "$CLUSTERS" --out "$OUT/nodes_input.tsv" \
      --id-mode "$ID_MODE"

echo "== Step 1: md5-unique node set + both-orientation FASTA =="
"$PY" scripts/step1_nodes.py --tsv "$OUT/nodes_input.tsv" --v1 "$V1" \
      ${UNION:+--v1-seqs "$UNION"} --outprefix "$OUT/combined"

echo "== Step 2: all-vs-all alignment ($ENGINE) =="
"$PY" scripts/step2_align.py --prefix "$OUT/combined" --engine "$ENGINE"

echo "== Step 3: filter + single-linkage partition =="
"$PY" scripts/step23_graph.py --prefix "$OUT/combined" --outdir "$OUT" --tag "$TAG" \
      --date "$DATE" --engine "$ENGINE"

echo "== Step 4: join components back onto every input row =="
"$PY" scripts/annotate_clusters.py --clusters "$CLUSTERS" \
      --assignment "$OUT/component_assignment_${TAG}_${DATE}.csv" \
      --out "$OUT/clusters_components_${DATE}.tsv" \
      --slim-out "$OUT_TSV"

echo "== done: $OUT_TSV (intermediates in $OUT/) =="
