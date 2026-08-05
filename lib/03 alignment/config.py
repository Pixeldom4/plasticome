#!/usr/bin/env python3
"""Shared configuration + the aligner helpers for the PETadex-alignment
consolidated pipeline.

One place for every path and constant the four pipeline scripts and the central
notebook agree on. The scientific constants (ID_MIN, EVALUE_MAX, the permissive
search / post-filter split) are the paper-faithful method carried forward from
v4-v8; do not change them without a version bump and a fresh REPORT.

Two aligners, one edge criterion:

  `diamond()`  DIAMOND blastp all-vs-all -- what the pipeline runs (step 3 of
               run_pipeline.bash). Search settings mirror the usearch method
               parameter for parameter: permissive at search time, everything
               decided by the >=30% aaid AND e<1e-5 post-filter. Runs natively
               when a `diamond` is on PATH, else bin/diamond via Docker.
  `usearch()`  usearch11 -allpairs_local -- kept for the Step-0 sanity fixture
               (213 v1 nodes -> the paper's 42 components / 3,178 edges), which
               is a usearch measurement and stays one. bin/usearch11 is a Linux
               ELF and only runs through Docker linux/amd64 on this arm64 macOS
               host -- see the project memory note.

`align()` dispatches on ALIGNER so a caller that does not care gets the
pipeline's engine.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

# ---- layout ---------------------------------------------------------------
# This module lives in lib/03 alignment/ and is run-agnostic: it resolves the repo
# root, then reads inputs from / writes outputs into whichever run is selected.
# Point PLASTICOME_RUN at another runs/<date>.<name>/ to re-run this pipeline
# against a different run's data without editing anything here.
ROOT = Path(__file__).resolve().parents[2]       # repo root == the Docker /d mount
LIB = ROOT / "lib" / "03 alignment"
BIN = ROOT / "bin"
SOURCES = ROOT / "source-data"                   # was `sources/` before 2026-07-30

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

# ---- aligners -------------------------------------------------------------
ALIGNER = os.environ.get("ALIGNER", "diamond")     # what align() dispatches to

DOCKER_IMG = "debian:bookworm-slim"
DOCKER_PLATFORM = "linux/amd64"
USEARCH_BIN = "bin/usearch11"                     # relative to ROOT (mounted at /d)
USEARCH_VERSION = "v11.0.667"
USEARCH_LOG = OUTPUTS / "usearch.log"

DIAMOND_BIN = "bin/diamond"                       # relative to ROOT (Linux ELF copy)
DIAMOND_THREADS = int(os.environ.get("DIAMOND_THREADS", "4"))
DIAMOND_LOG = OUTPUTS / "diamond.log"
# Permissive at search time; the post-filter in step23_graph.py is the only gate.
#   --ultra-sensitive     closest heuristic approximation of usearch's exhaustive allpairs
#   --max-target-seqs 0   DIAMOND keeps 25 hits/query by default, which silently truncates
#                         a dense all-vs-all -- the single most important flag here
#   --evalue 10           mirrors -acceptall: never drop a candidate edge at search time
#   --masking 0           usearch does not low-complexity mask; neither may we
#   --comp-based-stats 0  composition-based bitscore correction has no usearch counterpart
DIAMOND_FLAGS = ["--ultra-sensitive", "--max-target-seqs", "0", "--evalue", "10",
                 "--masking", "0", "--comp-based-stats", "0"]


def _rel(p: Path) -> str:
    """Path relative to ROOT, as a container sees it inside the /d mount."""
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


def _diamond_run(args: list, log_path: Path) -> None:
    """Invoke DIAMOND natively if it is on PATH, else bin/diamond under Docker.

    Path arguments are passed as Path objects and rewritten to the container's view
    when Docker is used; flags stay plain strings. Both forms need every file under
    ROOT, which is the single mount.
    """
    native = shutil.which("diamond")
    if native:
        cmd = [native] + [str(a) for a in args]
    else:
        cmd = ["docker", "run", "--rm", "--platform", DOCKER_PLATFORM,
               "-v", f"{ROOT}:/d", "-w", "/d", DOCKER_IMG, f"/d/{DIAMOND_BIN}"] \
            + [f"/d/{_rel(a)}" if isinstance(a, Path) else str(a) for a in args]
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)


def diamond(fasta: Path, out_pairs: Path, log_path: Path | None = None) -> Path:
    """All-vs-all `diamond blastp` of one FASTA against itself.

    Emits the same five fields in the same order as usearch(), so step23_graph.py
    reads either aligner's output unchanged: query, target, pct id, e-value, bits.

    Query == database is what makes this all-vs-all, and DIAMOND reports a pair in
    both directions within the one run -- so unlike usearch's -allpairs_local there is
    nothing asymmetric to correct for, and the reversed-FASTA second pass is not run.
    The max-bits-per-unordered-pair reduction in step23_graph.py still applies.

    Returns out_pairs. The .dmnd database is left beside it as an intermediate.
    """
    fasta, out_pairs = Path(fasta), Path(out_pairs)
    log_path = Path(log_path) if log_path else DIAMOND_LOG
    db = out_pairs.with_suffix(".dmnd")
    out_pairs.parent.mkdir(parents=True, exist_ok=True)
    _diamond_run(["makedb", "--in", fasta, "--db", db,
                  "--threads", DIAMOND_THREADS, "--quiet"], log_path)
    _diamond_run(["blastp", "--query", fasta, "--db", db, "--out", out_pairs,
                  "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore"]
                 + DIAMOND_FLAGS + ["--threads", DIAMOND_THREADS], log_path)
    return out_pairs


def align(fasta: Path, out_pairs: Path, engine: str | None = None) -> Path:
    """Run the pipeline's aligner (ALIGNER, default diamond) on one FASTA."""
    engine = engine or ALIGNER
    if engine == "diamond":
        return diamond(fasta, out_pairs)
    if engine == "usearch":
        return usearch(fasta, out_pairs)
    raise ValueError(f"unknown aligner: {engine}")


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
