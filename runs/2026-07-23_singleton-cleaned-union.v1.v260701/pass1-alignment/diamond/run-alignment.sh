#!/usr/bin/env bash
# All-vs-all DIAMOND alignment of the singleton-cleaned-union v1.v2 node set,
# tuned to mirror the usearch v4-v8 method (usearch11 -allpairs_local -acceptall,
# >=30% aaid AND e<1e-5 post-filter, single-linkage components).
#
# Runs on the EXACT SAME node set as the usearch run: seqs.faa is a copy of the
# usearch pipeline's combined.fasta (md5-unique nodes, m#### ids), so component /
# edge sets are directly comparable node-for-node. DIAMOND is a Linux ELF -> Docker.
set -euo pipefail
cd "$(dirname "$0")"
REPO="$PWD"; while [[ ! -d "$REPO/lib" || ! -d "$REPO/bin" ]]; do \
  REPO="$(dirname "$REPO")"; [[ "$REPO" == "/" ]] && { echo "repo root not found" >&2; exit 1; }; done
REL="${PWD#$REPO/}"                                    # this dir, repo-root-relative
DIAMOND="$REPO/bin/diamond.sh"
OUT=results
mkdir -p "$OUT"

# 1. Build DB from the shared node set (seqs.faa already staged = usearch combined.fasta)
"$DIAMOND" makedb --in "$OUT/seqs.faa" --db "$OUT/seqs" 2> "$OUT/diamond.log"

# 2. All-vs-all local alignment, permissive (filtering happens in cluster.py, not here)
"$DIAMOND" blastp \
  --query "$OUT/seqs.faa" --db "$OUT/seqs.dmnd" \
  --out "$OUT/allpairs.tsv" \
  --outfmt 6 qseqid sseqid pident evalue bitscore \
  --ultra-sensitive \
  --max-target-seqs 0 \
  --evalue 10 \
  --masking 0 \
  --comp-based-stats 0 \
  --threads 4 2>> "$OUT/diamond.log"

# 3. Post-filter (id>=30 AND e<1e-5) + single-linkage components
python3 cluster.py

echo "== done: results in $OUT/ =="
