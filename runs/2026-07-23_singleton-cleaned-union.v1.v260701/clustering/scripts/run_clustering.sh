#!/usr/bin/env bash
# Reference-seeded clustering (closed-then-open reference OTU picking) @ 90% identity.
# Engine: USEARCH v11.0.667, run in a Docker linux/amd64 container (arm64 macOS host).
# Dataset: singleton-cleaned union v1.v2 v260701 (13 curated seeds + 598 no-seed = 611 records).
# NOTE: this seed set uses PL18 (WP_015787089.1) and NOT PL485 — it differs from the
#       seed set described in the top-level readme.md's Results (which had PL485, no PL18).
#
# Layout (all paths relative to RUN_DIR):
#   inputs/  source FASTAs   work/  intermediate USEARCH records+logs   results/  deliverables
# Reproduce: ./scripts/run_clustering.sh  then  python3 scripts/build_membership.py
set -euo pipefail

RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"   # -> <run>/clustering/
REPO="$(cd "$RUN_DIR/../../.." && pwd -P)"                      # -> repo root
BIN_DIR="$REPO/bin"
ID=0.90

SEEDS="inputs/curated_seed_set.fasta"                                 # 13 curated reference seeds
NOSEED="inputs/singleton-cleaned-union.v1.v2.no-seeds.fasta"          # 598 remaining

# usearch wrapper: mounts the binary dir at /b and the run dir at /d
u() {
  docker run --rm --platform linux/amd64 \
    -v "${BIN_DIR}:/b" -v "${RUN_DIR}:/d" debian:stable-slim \
    /b/usearch11 "$@"
}

cd "$RUN_DIR"
mkdir -p work results

echo "== Phase 1: seed self-cluster diagnostic (non-redundancy @ ${ID}) =="
u -cluster_fast /d/${SEEDS} -id ${ID} -sort length \
  -centroids /d/work/seed-centroids.fasta \
  -uc /d/work/01-seed.uc \
  -log /d/work/01-seed.log

echo "== Phase 2: closed-reference recruitment (no-seed vs seed centroids) =="
u -usearch_global /d/${NOSEED} -db /d/work/seed-centroids.fasta -id ${ID} \
  -maxaccepts 0 -maxrejects 0 -top_hit_only \
  -uc /d/work/02-closedref.uc \
  -matched /d/work/matched.fasta \
  -notmatched /d/work/notmatched.fasta \
  -log /d/work/02-closedref.log

echo "== Phase 3: open-reference de-novo cluster of unmatched remainder =="
u -cluster_fast /d/work/notmatched.fasta -id ${ID} -sort length \
  -centroids /d/work/denovo-centroids.fasta \
  -uc /d/work/03-denovo.uc \
  -log /d/work/03-denovo.log

echo "== Phase 4a: emit all-centroids (seed centroids + de-novo centroids) =="
cat work/seed-centroids.fasta work/denovo-centroids.fasta \
  > results/singleton-cleaned-union.v1.v2.all-centroids.fasta

echo "== raw USEARCH steps complete; now run: python3 scripts/build_membership.py =="
