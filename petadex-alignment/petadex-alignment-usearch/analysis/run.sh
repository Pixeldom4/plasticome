#!/usr/bin/env bash
# End-to-end driver for the consolidated PETadex-alignment pipeline.
# Parity with the central notebook (PETadex_alignment.ipynb) -- same steps, same
# order, same outputs. Paper-faithful method (v4-v8): usearch v11.0.667
# -allpairs_local via Docker linux/amd64, permissive search, >=30% aaid AND e<1e-5
# post-filter, single-linkage components, both FASTA orientations.
#
# Usage:  PY=/path/to/python ./run.sh        (defaults to the plasticome conda env)
set -euo pipefail
cd "$(dirname "$0")"                # -> analysis/
PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python}"
DATE="${DATE:-2026-07-21}"
OUT=outputs
S0="$OUT/step0_check"

usearch() {  # $1 = fasta (analysis-relative), $2 = out pairs (analysis-relative)
  docker run --rm --platform linux/amd64 -v "$PWD":/d -w /d debian:bookworm-slim \
    /d/bin/usearch11 -allpairs_local "/d/$1" -userout "/d/$2" \
    -userfields query+target+id+evalue+bits -acceptall 2>>"$OUT/usearch.log"
}

echo "== Step 0: engine sanity (213 v1 nodes -> expect 42 components / 3178 edges) =="
$PY scripts/step0_sanity.py
usearch "$S0/v1_213.fasta"     "$S0/v1_213_pairs.tsv"
usearch "$S0/v1_213_rev.fasta" "$S0/v1_213_rev_pairs.tsv"
$PY scripts/step23_graph.py --prefix "$S0/v1_213" --outdir "$S0" --tag step0 --date "$DATE"

echo "== Step 1: assemble md5-unique node set (expect 484 nodes) =="
$PY scripts/step1_nodes.py

echo "== Step 2: all-vs-all usearch (both orientations) =="
usearch "$OUT/combined.fasta"     "$OUT/combined_pairs.tsv"
usearch "$OUT/combined_rev.fasta" "$OUT/combined_rev_pairs.tsv"

echo "== Step 3: partition + deliverables (expect 46 components) =="
$PY scripts/step23_graph.py --prefix "$OUT/combined" --outdir "$OUT" --tag petadex --date "$DATE"
$PY scripts/annotate_source.py \
    --assignment "$OUT/component_assignment_petadex_${DATE}.csv" \
    --out "$OUT/cleaned_pazy_final_components_${DATE}.csv"
echo "== done: outputs in $OUT/ =="
