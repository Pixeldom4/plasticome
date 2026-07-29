#!/usr/bin/env bash
# v8 pipeline driver. Same paper-faithful method as v4/v5/v6/v7:
#   usearch v11.0.667 -allpairs_local (Docker linux/amd64), permissive search,
#   >=30% aaid AND e<1e-5 post-filter, single-linkage components, both orientations.
# Input: v8/cleaned_pazy_final.tsv (snapshot of data/cleaned_pazy_final_v8-align.tsv,
#   the curated 484-node PAZy set; several v7 fragments re-sequenced, 6 rows dropped).
set -euo pipefail
cd "$(dirname "$0")/.."            # -> petadex-alignment/
ROOT="$PWD"
DATE=2026-07-20
IMG=debian:bookworm-slim

usearch() {  # $1 = fasta (repo-relative), $2 = out pairs (repo-relative)
  docker run --rm --platform linux/amd64 -v "$ROOT":/d -w /d "$IMG" \
    /d/bin/usearch11 -allpairs_local "/d/$1" -userout "/d/$2" \
    -userfields query+target+id+evalue+bits -acceptall 2>>"$ROOT/v8/usearch.log"
}

echo "== Step 0: engine sanity (213 v1 nodes -> expect 42 components) =="
python3 v8/step0_v8.py
usearch v8/step0_check/v1_213.fasta      v8/step0_check/v1_213_pairs.tsv
usearch v8/step0_check/v1_213_rev.fasta  v8/step0_check/v1_213_rev_pairs.tsv
python3 v8/step23_graph_v5.py --prefix v8/step0_check/v1_213 \
    --outdir v8/step0_check --tag step0 --date "$DATE"

echo "== Step 1: assemble md5-unique node set =="
python3 v8/step1_nodes_v8.py

echo "== Step 2: all-vs-all usearch (both orientations) =="
usearch v8/combined.fasta      v8/combined_pairs.tsv
usearch v8/combined_rev.fasta  v8/combined_rev_pairs.tsv

echo "== Step 3: partition + deliverables =="
python3 v8/step23_graph_v8.py --prefix v8/combined --outdir v8 --tag v8 --date "$DATE"
python3 v8/annotate_source_v8.py
echo "== done =="
