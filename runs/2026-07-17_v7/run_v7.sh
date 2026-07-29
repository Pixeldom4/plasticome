#!/usr/bin/env bash
# v7 pipeline driver. Same paper-faithful method as v4/v5/v6:
#   usearch v11.0.667 -allpairs_local (Docker linux/amd64), permissive search,
#   >=30% aaid AND e<1e-5 post-filter, single-linkage components, both orientations.
# Input: v7/cleaned_pazy_final.tsv (snapshot of the collapsed ~490-seq PAZy set).
set -euo pipefail
cd "$(dirname "$0")/.."            # -> petadex-alignment/
ROOT="$PWD"
DATE=2026-07-17
IMG=debian:bookworm-slim

usearch() {  # $1 = fasta (repo-relative), $2 = out pairs (repo-relative)
  docker run --rm --platform linux/amd64 -v "$ROOT":/d -w /d "$IMG" \
    /d/bin/usearch11 -allpairs_local "/d/$1" -userout "/d/$2" \
    -userfields query+target+id+evalue+bits -acceptall 2>>"$ROOT/v7/usearch.log"
}

echo "== Step 0: engine sanity (213 v1 nodes -> expect 42 components) =="
python3 v7/step0_v7.py
usearch v7/step0_check/v1_213.fasta      v7/step0_check/v1_213_pairs.tsv
usearch v7/step0_check/v1_213_rev.fasta  v7/step0_check/v1_213_rev_pairs.tsv
python3 v7/step23_graph_v5.py --prefix v7/step0_check/v1_213 \
    --outdir v7/step0_check --tag step0 --date "$DATE"

echo "== Step 1: assemble md5-unique node set =="
python3 v7/step1_nodes_v7.py

echo "== Step 2: all-vs-all usearch (both orientations) =="
usearch v7/combined.fasta      v7/combined_pairs.tsv
usearch v7/combined_rev.fasta  v7/combined_rev_pairs.tsv

echo "== Step 3: partition + deliverables =="
python3 v7/step23_graph_v7.py --prefix v7/combined --outdir v7 --tag v7 --date "$DATE"
python3 v7/annotate_source_v7.py
echo "== done =="
