#!/usr/bin/env bash
# All-vs-all alignment for the singleton-cleaned-union v1.v2 all-centroids set, using
# the frozen PETadex component pipeline (petadex-alignment-usearch/analysis/, v4-v8
# paper-faithful method):
#   usearch11 -allpairs_local via Docker linux/amd64, both orientations,
#   best-HSP by max-bits, post-filter >=30% aaid AND e-value < 1e-5, single-linkage.
# Consolidated layout: this lives at petadex-alignment/<dataset>/usearch/; the shared
# analysis/ tooling stays in petadex-alignment-usearch/. Results land in ./results/.
# The resulting results/combined.fasta is the shared node set the sibling diamond/ run
# reuses as seqs.faa, so the two partitions are node-for-node comparable.
set -euo pipefail
cd "$(dirname "$0")"
REPO="$PWD"; while [[ ! -d "$REPO/lib" || ! -d "$REPO/bin" ]]; do \
  REPO="$(dirname "$REPO")"; [[ "$REPO" == "/" ]] && { echo "repo root not found" >&2; exit 1; }; done
REL="${PWD#$REPO/}"                                    # this dir, repo-root-relative
HERE_ABS="$(pwd -P)"

REL="${HERE_ABS#"$ROOT"/}"                             # -> <dataset>/usearch
ANALYSIS="$REPO/lib/alignment"
PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python}"
TAG="${TAG:-singleton}"
DATE="${DATE:-2026-07-23}"
FASTA=singleton-cleaned-union.v1.v2.all-centroids.fasta
OUT=results
mkdir -p "$OUT"

usearch() {  # $1 = fasta, $2 = out pairs  (both relative to this usearch/ dir)
  docker run --rm --platform linux/amd64 -v "$ROOT":/d -w /d debian:bookworm-slim \
    /d/bin/usearch11 -allpairs_local "/d/$REL/$1" \
    -userout "/d/$REL/$2" \
    -userfields query+target+id+evalue+bits -acceptall 2>>"$OUT/usearch.log"
}

echo "== Step 1: centroids FASTA -> curated-shape TSV =="
$PY "$ANALYSIS/scripts/centroids_to_tsv.py" --fasta "$FASTA" --out "$OUT/input.tsv"

echo "== Step 1b: assemble md5-unique node set (+ both-orientation FASTAs) =="
$PY "$ANALYSIS/scripts/step1_nodes.py" \
    --tsv "$OUT/input.tsv" \
    --v1  "$REPO/sources/plasticome.v1.csv" \
    --outprefix "$OUT/combined"

echo "== Step 2: all-vs-all usearch (both orientations) =="
usearch "$OUT/combined.fasta"     "$OUT/combined_pairs.tsv"
usearch "$OUT/combined_rev.fasta" "$OUT/combined_rev_pairs.tsv"

echo "== Step 3: partition + deliverables =="
$PY "$ANALYSIS/scripts/step23_graph.py" \
    --prefix "$OUT/combined" --outdir "$OUT" --tag "$TAG" --date "$DATE"
echo "== done: results in $OUT/ =="
