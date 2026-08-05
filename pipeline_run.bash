#!/usr/bin/env bash
# Plasticome end-to-end driver: source-data -> union -> clusters -> components -> FASTA.
#
# One run directory, four numbered deliverables, one sidecar per step:
#
#   runs/<date>.<run-name>/
#   ├── 01-union.tsv               610 rows   v1.1 u v260701, sequences backfilled
#   ├── 02-clusters.tsv            411 rows   reference-seeded 90% clusters
#   ├── 02-clusters.intermediates/              edge graph + seed/closed-ref/de-novo scratch
#   ├── 03-alignment.tsv           411 rows   + component_id
#   ├── 03-alignment.intermediates/             all-pairs alignment + graph scratch
#   ├── 04-<run-name>.fasta        411 records centroid FASTA
#   └── summary.md                            counts and provenance for the above
#
# Per-step logs (aligner stdout/stderr, phase logs) land in that step's sidecar.
# summary.md is rewritten on every invocation from what is on disk, so it also
# describes a partial run (steps not built are named as such).
#
# Usage:
#   ./run_pipeline.sh <run-name>              # new run under runs/<today>.<run-name>/
#   ./run_pipeline.sh --run-dir runs/<dir>    # resume / re-run an existing directory
#   ./run_pipeline.sh <run-name> --from 3     # only steps 3-4
#
# Options:
#   --run-dir DIR     operate on an existing run directory instead of creating one
#   --from N          first step to run  (default: 1)
#   --to N            last step to run   (default: 4)
#   --only N          shorthand for --from N --to N
#   --force           re-run steps whose deliverable already exists
#   --id PCT          clustering identity threshold (default: 0.90)
#   --engine NAME     diamond (default) or usearch, for BOTH steps 2 and 3
#   --v1 FILE         v1.1 table        (default: source-data/plasticome.v1.1/...)
#   --v2 FILE         v260701 table     (default: source-data/plasticome.v260701/...)
#   --v2-seq MODE     how to read v260701's aa_sequence column: from-accession
#                     (default) keeps the attached sequence only on rows with no
#                     accession and fetches the rest in step 1; as-given takes the
#                     column verbatim, for a v2 table already retrieved from accession
#   --seeds FILE      curated seeds     (default: source-data/plasticome-curated-seeds/...)
#   --dry-run         print the commands without running them
#   -h, --help        this message
#
# Env:
#   PY=...            python 3.9+ (default: the plasticome conda env, else python3)
#   NCBI_API_KEY=...  lifts step 1's NCBI rate limit from 3 to 10 req/s
#
# Steps 2 and 3 shell out to DIAMOND for the clustering and the all-vs-all
# alignment respectively. DIAMOND runs natively when one is on PATH (`brew install
# diamond`), and otherwise falls back to the bundled bin/diamond, which is a Linux
# x86-64 ELF and needs Docker linux/amd64 up on this arm64 host. --engine usearch
# selects the previous engine for both steps and always requires Docker.
# Everything is resumable: a step whose deliverable already exists is skipped
# unless --force, so an interrupted run picks up where it stopped.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$REPO"

# ---------------------------------------------------------------- defaults ---
# The conda env is where the pipeline's deps live; fall back to whatever python3
# is on PATH so the script is still usable on a fresh checkout.
PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python}"
[ -x "$PY" ] || PY="$(command -v python3 || true)"

V1="$REPO/source-data/plasticome.v1.1/plasticome.v1.1.no-seq.tsv"
V2="$REPO/source-data/plasticome.v260701/cleaned_pazy-260701-singletons.tsv"
SEEDS="$REPO/source-data/plasticome-curated-seeds/plasticome-curated-seeds.tsv"

# The singletons table carries the sequence attached to the PAZy record, which is
# the right one only where the record has no accession; everywhere else step 1
# blanks it and refetches from the accession. A v2 table whose sequences were
# already retrieved from accession wants --v2-seq as-given instead.
V2SEQ="from-accession"

RUN_NAME=""; RUN_DIR=""; FROM=1; TO=4; FORCE=""; DRY=""; ID="0.90"; ENGINE="diamond"

die() { printf 'error: %s\n' "$*" >&2; exit 1; }
usage() { sed -n '2,/^set -euo/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//; $d'; }

# ------------------------------------------------------------------- args ----
while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help)  usage; exit 0 ;;
    --run-dir)  RUN_DIR="${2:?--run-dir needs a path}"; shift 2 ;;
    --from)     FROM="${2:?--from needs a step number}"; shift 2 ;;
    --to)       TO="${2:?--to needs a step number}"; shift 2 ;;
    --only)     FROM="${2:?--only needs a step number}"; TO="$FROM"; shift 2 ;;
    --force)    FORCE=1; shift ;;
    --dry-run)  DRY=1; shift ;;
    --id)       ID="${2:?--id needs a value}"; shift 2 ;;
    --engine)   ENGINE="${2:?--engine needs diamond or usearch}"; shift 2 ;;
    --v1)       V1="${2:?--v1 needs a path}"; shift 2 ;;
    --v2)       V2="${2:?--v2 needs a path}"; shift 2 ;;
    --v2-seq)   V2SEQ="${2:?--v2-seq needs as-given or from-accession}"; shift 2 ;;
    --seeds)    SEEDS="${2:?--seeds needs a path}"; shift 2 ;;
    -*)         die "unknown option: $1 (try --help)" ;;
    *)          [ -z "$RUN_NAME" ] || die "unexpected argument: $1"
                RUN_NAME="$1"; shift ;;
  esac
done

case "$FROM$TO" in *[!1-4]*) die "--from/--to must be between 1 and 4" ;; esac
[ "$FROM" -le "$TO" ] || die "--from ($FROM) is after --to ($TO)"
case "$V2SEQ" in
  as-given|from-accession) ;;
  *) die "--v2-seq must be as-given or from-accession (got: $V2SEQ)" ;;
esac
case "$ENGINE" in
  diamond|usearch) ;;
  *) die "--engine must be diamond or usearch (got: $ENGINE)" ;;
esac

# ------------------------------------------------------------- run directory --
if [ -n "$RUN_DIR" ]; then
  [ -z "$RUN_NAME" ] || die "give either <run-name> or --run-dir, not both"
  RUN_DIR="$(cd "$(dirname "$RUN_DIR")" && pwd -P)/$(basename "$RUN_DIR")"
  [ -d "$RUN_DIR" ] || die "no such run directory: $RUN_DIR"
  # Strip the leading date so 04's filename matches the run's identity.
  RUN_NAME="$(basename "$RUN_DIR")"; RUN_NAME="${RUN_NAME#[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9].}"
else
  [ -n "$RUN_NAME" ] || die "need a <run-name> or --run-dir (try --help)"
  case "$RUN_NAME" in */*) die "run name must not contain '/'" ;; esac
  RUN_DIR="$REPO/runs/$(date +%F).$RUN_NAME"
fi

# Docker mounts the repo root at /d, so anything outside it is invisible to usearch.
case "$RUN_DIR" in "$REPO"/*) ;; *) die "run directory must sit under $REPO" ;; esac

# Deliverable paths. A run created before this driver existed may have named a
# step differently -- the 07-30 union-spec run calls step 1
# "01-plasticome_v1.v260701-union.tsv", not "01-union.tsv". Adopt such a file
# when the step prefix matches exactly one candidate, so resuming an older run
# reuses its deliverable instead of building a second one beside it.
# Array-free for bash 3.2, which is what /bin/bash is on macOS.
adopt() {  # $1 = canonical path, $2 = glob matching the step's deliverable
  local canon="$1" pat="$2" f found="" n=0
  if [ -e "$canon" ]; then printf '%s' "$canon"; return; fi
  for f in "$RUN_DIR"/$pat; do
    [ -e "$f" ] || continue          # unmatched globs stay literal; skip them
    found="$f"; n=$((n + 1))
  done
  if [ "$n" -eq 1 ]; then printf '%s' "$found"; else printf '%s' "$canon"; fi
}

S1="$(adopt "$RUN_DIR/01-union.tsv"        '01-*.tsv')"
S2="$(adopt "$RUN_DIR/02-clusters.tsv"     '02-*.tsv')"
S3="$(adopt "$RUN_DIR/03-alignment.tsv"    '03-*.tsv')"
S4="$(adopt "$RUN_DIR/04-$RUN_NAME.fasta"  '04-*.fasta')"

# --------------------------------------------------------------- preflight ---
want() { [ "$1" -ge "$FROM" ] && [ "$1" -le "$TO" ]; }

[ -n "$PY" ] && { [ -x "$PY" ] || command -v "$PY" >/dev/null; } || die "no usable python (set PY=)"

if want 1; then
  [ -f "$V1" ]    || die "missing v1.1 table: $V1"
  [ -f "$V2" ]    || die "missing v260701 table: $V2"
fi
if want 2; then
  [ -f "$SEEDS" ] || die "missing curated seeds: $SEEDS"
fi
if want 2 || want 3; then
  # Docker is only needed for the Linux-ELF binaries: always for usearch, and for
  # diamond only when there is no native one on PATH.
  NEED_DOCKER=1
  if [ "$ENGINE" = diamond ]; then
    if command -v diamond >/dev/null; then
      NEED_DOCKER=""
    else
      [ -f "$REPO/bin/diamond" ] || die "no diamond on PATH and missing $REPO/bin/diamond"
    fi
  else
    [ -f "$REPO/bin/usearch11" ] || die "missing $REPO/bin/usearch11"
  fi
  [ -z "$NEED_DOCKER" ] || [ -n "$DRY" ] || docker info >/dev/null 2>&1 \
    || die "Docker is not running — bin/$ENGINE is a Linux ELF binary and needs linux/amd64"
fi
# Steps 3 and 4 consume the previous step's deliverable; if it is not being built
# this invocation, it has to already be on disk.
want 3 && ! want 2 && [ ! -f "$S2" ] && die "step 3 needs $S2 (run --from 2)"
want 4 && ! want 3 && [ ! -f "$S3" ] && die "step 4 needs $S3 (run --from 3)"

[ -n "$DRY" ] || mkdir -p "$RUN_DIR"

# ----------------------------------------------------------------- helpers ---
# Echo the command, then run it unless --dry-run. Each argument is shell-quoted
# with %q rather than printed raw: the lib/ directories contain spaces, so a raw
# echo produces a line that looks copy-pasteable but splits into the wrong argv.
run() {
  local a q=""
  for a in "$@"; do q="$q $(printf '%q' "$a")"; done
  printf '  $%s\n' "$q"
  [ -n "$DRY" ] || "$@"
}

# Skip a step whose deliverable is already on disk. Steps 2 and 3 cost minutes of
# usearch time, so re-running a four-step pipeline to redo only the last one is
# the common case, not the exception.
skip() {  # $1 = step number, $2 = deliverable
  if [ -e "$2" ] && [ -z "$FORCE" ]; then
    printf '== step %s: SKIP — %s exists (--force to rebuild)\n' "$1" "${2#"$REPO"/}"
    return 0
  fi
  return 1
}

banner() { printf '\n== step %s: %s\n' "$1" "$2"; }

rows() {  # data rows in a TSV, or records in a FASTA; "-" if absent
  [ -f "$1" ] || { printf -- '-'; return; }
  case "$1" in
    *.fasta) grep -c '^>' "$1" ;;
    *)       awk 'END{print NR-1}' "$1" ;;
  esac
}

# ------------------------------------------------------------------- steps ---
printf '== plasticome pipeline\n'
printf '   run      %s\n' "${RUN_DIR#"$REPO"/}"
printf '   steps    %s-%s%s%s\n' "$FROM" "$TO" "${FORCE:+ (force)}" "${DRY:+ (dry run)}"
printf '   python   %s\n' "$PY"
(want 2 || want 3) && printf '   engine   %s\n' "$ENGINE"
want 1 && printf '   v2       %s (--v2-seq %s)\n' "${V2#"$REPO"/}" "$V2SEQ"

if want 1 && ! skip 1 "$S1"; then
  banner 1 "union — v1.1 ∪ v260701, then backfill sequences from accession"
  [ -n "${NCBI_API_KEY:-}" ] || printf '  note: NCBI_API_KEY unset — fetches are rate-limited to 3/s\n'
  run "$PY" "$REPO/lib/01 union/build_union.py" "$V1" "$V2" --v2-seq "$V2SEQ" -o "$S1"
  # fetch_sequences.py fills the table in place; it is a no-op when every row
  # already carries a sequence, so it is safe to re-run.
  run "$PY" "$REPO/lib/01 union/fetch_sequences.py" "$S1"
fi

if want 2 && ! skip 2 "$S2"; then
  banner 2 "reference-seeded clustering at ${ID} identity  [$ENGINE]"
  run "$PY" "$REPO/lib/02 clustering/cluster_reference_seeded.py" \
      "$S1" "$SEEDS" --id "$ID" --engine "$ENGINE" -o "$S2"
fi

if want 3 && ! skip 3 "$S3"; then
  banner 3 "all-vs-all alignment of the centroids → component_id  [$ENGINE]"
  # V1 is passed explicitly even though run_from_clusters.sh now defaults it
  # correctly: a missing overlay is NOT fatal there — step1_nodes.py only warns and
  # leaves v1_component/cath blank — so the failure mode is a silently poorer
  # result rather than an error. Naming it here keeps that impossible to miss.
  # OUT_TSV likewise overrides its `03 alignment.tsv` (space) default.
  # UNION lets the v1 overlay key on sequence md5 even though the v1.1 table is now
  # a .no-seq export: without it the overlay falls back to accession-only matching
  # and annotates far fewer clusters (161 -> 114 on the 08-05 node set).
  run env PY="$PY" ENGINE="$ENGINE" V1="$V1" UNION="$S1" OUT_TSV="$S3" \
      "$REPO/lib/03 alignment/run_from_clusters.sh" "$S2"
fi

if want 4 && ! skip 4 "$S4"; then
  banner 4 "centroid FASTA"
  run "$PY" "$REPO/lib/fasta/clusters_to_fasta.py" "$S3" -o "$S4"
fi

# ------------------------------------------------------------------ summary --
# summary.md is rebuilt from whatever is on disk on every invocation, including
# one that skipped every step -- it is a read-only description of the run, so
# there is nothing to preserve and a stale one would be worse than none.
printf '\n== summary\n'
run "$PY" "$REPO/lib/summary/summarize_run.py" --run-dir "$RUN_DIR" \
    --s1 "$S1" --s2 "$S2" --s3 "$S3" --s4 "$S4" \
    --v1 "$V1" --v2 "$V2" --v2-seq "$V2SEQ" --seeds "$SEEDS" \
    --id "$ID" --date "$(date +%F)"
printf '\n   %-38s %8s\n' "deliverable" "rows"
for f in "$S1" "$S2" "$S3" "$S4"; do
  printf '   %-38s %8s\n' "$(basename "$f")" "$(rows "$f")"
done
# Step 1 is 610 on the singletons table (473 v260701 rows + 137 unmatched v1.1).
# 607 / 411 / 411 / 411 was the previous default input,
# cleaned_pazy-260701_retrieving_from_accession.tsv; steps 2-4 have not been
# re-measured on the singletons table, so only step 1's count is stated here.
printf '\n   expected step 1 on the reference inputs: 610\n'
printf '   done: %s\n' "${RUN_DIR#"$REPO"/}"
