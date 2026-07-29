#!/usr/bin/env python3
"""Shared configuration + the usearch(Docker) helper for the PETadex-alignment
consolidated pipeline.

One place for every path and constant the four pipeline scripts and the central
notebook agree on. The scientific constants (ID_MIN, EVALUE_MAX, usearch flags,
both-orientation search) are the paper-faithful method carried forward from
v4-v8; do not change them without a version bump and a fresh REPORT.

usearch is a Linux ELF binary (bin/usearch11) and only runs here through Docker
linux/amd64 on this arm64 macOS host -- see the project memory note. `usearch()`
wraps that call so both run.sh and the notebook invoke it identically.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

# ---- layout ---------------------------------------------------------------
# This module lives in lib/alignment/ and is run-agnostic: it resolves the repo
# root, then reads inputs from / writes outputs into whichever run is selected.
# Point PLASTICOME_RUN at another runs/<date>_<name>/ to re-run this pipeline
# against a different run's data without editing anything here.
ROOT = Path(__file__).resolve().parents[2]       # repo root == the Docker /d mount
LIB = ROOT / "lib" / "alignment"
BIN = ROOT / "bin"
SOURCES = ROOT / "sources"

RUN = Path(os.environ.get("PLASTICOME_RUN", ROOT / "runs" / "2026-07-20_collapse-v9-v10"))
DATA = RUN / "inputs"
OUTPUTS = RUN / "alignment-usearch"
STEP0_DIR = OUTPUTS / "step0_check"

INPUT_TSV = DATA / "cleaned_pazy_final.tsv"       # frozen 484-node curated PAZy set
V1_CSV = SOURCES / "plasticome.v1.csv"            # v1 overlay (descriptive only)

# ---- provenance -----------------------------------------------------------
# md5 of the frozen input this pipeline is pinned to (v8 curated set).
INPUT_MD5 = "1d13e83691cc767c7ba9635c9bd2ed60"
TAG = "petadex"
DATE = "2026-07-21"

# ---- scientific constants (paper-faithful, v4-v8) -------------------------
ID_MIN = 30.0          # >= 30% amino-acid identity
EVALUE_MAX = 1e-5      # AND e-value < 1e-5
EVALUE_SCALE = 1.0     # threshold e-value as reported at this run's db (no down-scaling)

# ---- usearch / Docker -----------------------------------------------------
DOCKER_IMG = "debian:bookworm-slim"
DOCKER_PLATFORM = "linux/amd64"
USEARCH_BIN = "bin/usearch11"                     # relative to ROOT (mounted at /d)
USEARCH_VERSION = "v11.0.667"
USEARCH_LOG = OUTPUTS / "usearch.log"


def _rel(p: Path) -> str:
    """Path relative to ROOT, as usearch sees it inside the /d mount."""
    return str(Path(p).resolve().relative_to(ROOT))


def usearch(fasta: Path, out_pairs: Path) -> Path:
    """Run `usearch11 -allpairs_local` on one FASTA, writing best-HSP userout.

    Emits query+target+id+evalue+bits with -acceptall (post-filtering happens in
    step23_graph.py, never as a search parameter). Both file args must live under
    ROOT so the single -v mount reaches them. Appends usearch chatter to
    OUTPUTS/usearch.log. Returns out_pairs.
    """
    fasta, out_pairs = Path(fasta), Path(out_pairs)
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    cmd = [
        "docker", "run", "--rm", "--platform", DOCKER_PLATFORM,
        "-v", f"{ROOT}:/d", "-w", "/d", DOCKER_IMG,
        f"/d/{USEARCH_BIN}", "-allpairs_local", f"/d/{_rel(fasta)}",
        "-userout", f"/d/{_rel(out_pairs)}",
        "-userfields", "query+target+id+evalue+bits", "-acceptall",
    ]
    with USEARCH_LOG.open("a") as log:
        subprocess.run(cmd, check=True, stderr=log)
    return out_pairs


if __name__ == "__main__":
    # `python config.py` -> quick environment + provenance sanity print.
    import hashlib
    print(f"ROOT          {ROOT}")
    print(f"INPUT_TSV     {INPUT_TSV}  (exists={INPUT_TSV.exists()})")
    if INPUT_TSV.exists():
        got = hashlib.md5(INPUT_TSV.read_bytes()).hexdigest()
        print(f"input md5     {got}  ({'OK' if got == INPUT_MD5 else 'MISMATCH!'})")
    print(f"usearch bin   {(BIN / 'usearch11').exists()}")
    print(f"criterion     pct_id>={ID_MIN} AND evalue<{EVALUE_MAX:g}")
