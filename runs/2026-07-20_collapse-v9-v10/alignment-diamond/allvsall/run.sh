#!/usr/bin/env bash
# All-vs-all DIAMOND alignment of substrings-collapsed set, tuned to mirror the
# usearch v4-v8 method (usearch11 -allpairs_local -acceptall, >=30% aaid AND
# e<1e-5 post-filter, single-linkage components).
#
# Parameter choices deliberately track usearch; see the header of cluster.py and
# the notes below. DIAMOND is a Linux ELF -> run through Docker linux/amd64.
set -euo pipefail
cd "$(dirname "$0")/.."                 # -> petadex-alignment-diamond/
IN=substrings-collapsed/plasticome-substrings-collapsed.tsv
OUT=allvsall

D() { docker run --rm --platform linux/amd64 -v "$PWD":/data -w /data debian:stable-slim ./diamond "$@"; }

# 1. TSV -> FASTA (identifier as header, aa_sequence as body)
python3 - "$IN" "$OUT/seqs.faa" <<'PY'
import csv, sys
src, out = sys.argv[1], sys.argv[2]
with open(src) as f, open(out, "w") as w:
    for row in csv.DictReader(f, delimiter="\t"):
        seq = row["aa_sequence"].strip().replace(" ", "")
        if seq:
            w.write(f">{row['identifier'].strip()}\n{seq}\n")
PY

# 2. Build DB
D makedb --in "$OUT/seqs.faa" --db "$OUT/seqs"

# 3. All-vs-all local alignment, permissive (filter happens in cluster.py)
D blastp \
  --query "$OUT/seqs.faa" --db "$OUT/seqs.dmnd" \
  --out "$OUT/allpairs.tsv" \
  --outfmt 6 qseqid sseqid pident evalue bitscore \
  --ultra-sensitive \
  --max-target-seqs 0 \
  --evalue 10 \
  --masking 0 \
  --comp-based-stats 0 \
  --threads 4

# 4. Post-filter (id>=30 AND e<1e-5) + single-linkage components
python3 "$OUT/cluster.py"
