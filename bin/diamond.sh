#!/usr/bin/env bash
# Wrapper to run the Linux DIAMOND binary via Docker (linux/amd64) on arm64 macOS.
# Mounts the plasticome project root as /data so the binary and your data are both
# visible, and sets the container working directory to match your current shell dir
# so relative --query/--db/--out paths just work.
#
# Usage (from anywhere inside the plasticome tree):
#   diamond.sh version
#   diamond.sh blastp --query q.faa --db ref.dmnd --out hits.tsv
set -euo pipefail

# Directory this script lives in (…/petadex-alignment-diamond).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# plasticome project root = parent of the diamond folder.
ROOT_DIR="$(dirname "$SCRIPT_DIR")"

# Absolute path to the binary and to the caller's cwd, then make both relative to ROOT.
BIN="$SCRIPT_DIR/diamond"
CWD="$(pwd -P)"

case "$CWD/" in
  "$ROOT_DIR"/*) ;;
  *)
    echo "diamond.sh: run this from inside $ROOT_DIR (current: $CWD)" >&2
    exit 1
    ;;
esac

# Path of caller's cwd relative to ROOT_DIR, mapped under /data in the container.
REL_CWD="${CWD#"$ROOT_DIR"}"
CONTAINER_WD="/data${REL_CWD}"
REL_BIN="${BIN#"$ROOT_DIR"}"

exec docker run --rm --platform linux/amd64 \
  -v "$ROOT_DIR":/data \
  -w "$CONTAINER_WD" \
  debian:stable-slim \
  "/data${REL_BIN}" "$@"
