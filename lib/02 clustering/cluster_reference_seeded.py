#!/usr/bin/env python3
"""Reference-seeded clustering of a plasticome union table — one row per cluster.

Implements `reference-seed-clustering.md` (closed-then-open reference OTU picking)
end to end, from the two TSVs to a cluster table:

  Phase 1  seed self-cluster diagnostic  — are the curated seeds non-redundant at $id?
  Phase 2  closed-reference recruitment  — union rows that hit a seed join that seed
  Phase 3  open-reference (de novo)      — the unmatched remainder is clustered fresh
  Phase 4  merge & emit                  — one row per cluster, incl. the centroid sequence

Every curated seed stays a representative; the criteria applied to everything else are
unchanged.

Engines
-------
`--engine diamond` (default) — DIAMOND v2.x, per `reference-seed-clustering.md`:
one all-vs-all `blastp` edge graph is built ONCE at `--id` and reused by every phase,
phases 1 and 3 are `greedy-vertex-cover` over that graph, and phase 2 recruits along it.
DIAMOND runs natively (a `diamond` on PATH, else the bundled `bin/diamond` through
Docker `linux/amd64`, since that copy is a Linux x86-64 ELF).

`--engine usearch` — USEARCH v11 `cluster_fast` / `usearch_global`, `-sort length`; the
engine the canonical 413-cluster run used. `bin/usearch11` is a Linux x86-64 ELF, so on
arm64 macOS it only runs inside Docker `linux/amd64` (Docker Desktop must be running);
pass `--no-docker` on a Linux x86-64 host.

Coverage differs between the two by necessity, so `--member-cov` defaults per engine:

  diamond  0.90 — DIAMOND identity is measured over a LOCAL HSP, so identity alone
                  would merge two sequences that share one well-conserved 40-residue
                  stretch. Coverage is what makes the 90% threshold mean "these are the
                  same protein". Applied as `--query-or-subject-cover` on the edge graph
                  and `--member-cover` in greedy-vertex-cover / recruitment. Pass
                  `--member-cov 0` to drop the rule (not recommended).
  usearch  unset — `usearch_global`/`cluster_fast` identity is already over a
                  near-global alignment, and recruiting on identity alone is what the
                  canonical 90% run did. `--member-cov` adds the doc's rule (`-query_cov`,
                  applied identically in all three phases) and recruits fewer sequences to
                  the seeds — on the v1/v260701 union at 90%, 17 instead of 28.

Greedy ordering is engine-inherent and worth restating in any report: `cluster_fast`
orders by length (longest = centroid), `greedy-vertex-cover` by graph degree.

Inputs
------
UNION_TSV   union table  : enzyme_name, accession, pazy_id, aa_sequence, source
SEEDS_TSV   curated seeds: enzyme_name, accession, pazy_id, justification,
                           aa_sequence, and optionally plasticome_id
Extra columns in either file are ignored.

A curated seed is an ANNOTATION on a union row, not a record of its own: the union row
carrying the seed's exact sequence is PROMOTED to be that seed's centroid, keeping its
`U####` label and gaining the seed's curated metadata. Any *further* union row that is
that same seed (same accession, version-insensitive, or an identical sequence) is folded
into its cluster as a member with `dup` in place of a percent identity. So every input
row is accounted for exactly once and the universe equals the union — building each seed
as a separate record instead double-counted it, since every seed is already a union row.

A seed with no exact-sequence union row has nothing to promote: it stays its own `S###`
record and genuinely adds one to the universe.

Output
------
One row per 90% cluster, ordered by `cluster_id` ascending:

  cluster_id            1..N. Reference (seed) clusters are numbered FIRST, in seed
                        order, so with the 13 curated seeds they are always 1-13 and
                        every id above that is a de-novo cluster. This ordering is the
                        only marker of a cluster's origin — there is no `origin`
                        column — so it is load-bearing, not cosmetic. (Under
                        `--on-seed-collapse=drop` a redundant seed is demoted and the
                        reference block shrinks accordingly; the default `keep` holds
                        the block at one cluster per seed.)
  size                  members incl. the centroid
  rep_*                 the representative: label, enzyme_name, accession, pazy_id,
                        plasticome_id (seeds only), source, seq_len
  member_*              the non-centroid members, "; "-joined and index-aligned:
                        labels, enzyme_names, accessions, pct ids to the centroid
  rep_aa_sequence       amino-acid sequence of the representative centroid

Labels are `U####|accession` (union rows, including the promoted seed centroids) and
`S###|accession` (only a seed that is absent from the union); the numeric prefix keeps
them unique when accessions repeat or are blank. A reference cluster is marked by its
`rep_plasticome_id` / `rep_source=seed` and by its position in the id order, not by the
letter its label starts with.

Examples
--------
  python3 cluster_reference_seeded.py union.tsv seeds.tsv -o clusters.tsv

  python3 cluster_reference_seeded.py union.tsv seeds.tsv -o clusters.tsv \\
      --id 0.95 --centroids-out all-centroids.fasta --workdir runs/<run>/clustering/work

  python3 cluster_reference_seeded.py union.tsv seeds.tsv -o clusters.tsv \\
      --engine usearch          # the pre-2026-08 engine, identity-only recruitment
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

UNION_COLS = ["enzyme_name", "accession", "pazy_id", "aa_sequence", "source"]
SEED_COLS = ["enzyme_name", "accession", "pazy_id", "aa_sequence"]
# plasticome_id is provenance, not a join key: seeds match union rows by exact
# sequence, and the merge below keeps the union's id wherever the seed's is blank.
# A seed table that omits the column therefore loses nothing, so it is read as
# optional and defaulted to "" rather than required.
SEED_OPTIONAL_COLS = ["plasticome_id"]

OUT_COLS = [
    "cluster_id", "size",
    "rep_label", "rep_enzyme_name", "rep_accession", "rep_pazy_id",
    "rep_plasticome_id", "rep_source", "rep_seq_len",
    "member_labels", "member_enzyme_names", "member_accessions", "member_pct_ids",
    "rep_aa_sequence",
]
JOIN = "; "


# --------------------------------------------------------------------------- io

def read_table(path: Path, required: list[str], what: str,
               optional: list[str] = ()) -> list[dict]:
    """Read a TSV into dicts, keeping the FIRST column of any repeated header name.

    csv.DictReader keeps the *last*, which silently blanks every row when a sheet
    carries a stray empty duplicate (the v260701 sheet has done exactly that with a
    second `aa_sequence`), so the header is indexed by hand here.

    `optional` columns are read when the header carries them and defaulted to ""
    when it does not, so callers can index them without a KeyError either way.
    """
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"error: {what} table {path} is empty")
    header = rows[0]
    idx: dict[str, int] = {}
    for i, col in enumerate(header):
        idx.setdefault(col.strip(), i)
    missing = [c for c in required if c not in idx]
    if missing:
        sys.exit(f"error: {what} table {path} is missing column(s): {', '.join(missing)}\n"
                 f"       found: {', '.join(header)}")
    dup = sorted({c for c, n in Counter(h.strip() for h in header).items() if n > 1})
    if dup:
        print(f"warning: {what} table has duplicate column(s) {dup}; using the first of each")
    absent = [c for c in optional if c not in idx]
    out = []
    for row in rows[1:]:
        if not any(cell.strip() for cell in row):
            continue
        rec = {c: (row[i].strip() if i < len(row) else "") for c, i in idx.items()}
        rec.update({c: "" for c in absent})
        out.append(rec)
    return out


def normalize(seq: str) -> str:
    """Strip every non-letter (gaps, `*` stops, whitespace) and upper-case."""
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def acc_key(acc: str) -> str:
    """Accession match key: version-stripped and upper-cased; blank/`-` -> ''."""
    a = (acc or "").strip()
    if a in ("", "-"):
        return ""
    return a.split(".")[0].upper()


def pct(value: str) -> float:
    """Sort key for a .uc identity field; `dup` (folded seed row) ranks as 100%."""
    if value == "dup":
        return 100.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def safe(text: str) -> str:
    """Label-safe token: USEARCH truncates a label at the first whitespace."""
    return re.sub(r"\s+", "_", (text or "").strip()) or "noacc"


def write_fasta(path: Path, records: list[dict], wrap: int = 60) -> None:
    with open(path, "w") as fh:
        for r in records:
            fh.write(f">{r['label']}\n")
            seq = r["seq"]
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")


# ---------------------------------------------------------------------- engines
#
# Both engines expose the same two operations, so phases 1-4 read identically
# whichever one is selected:
#
#   cluster(records, tag)  -> (centroid labels in cluster order,
#                              {centroid label: [(member label, pct id), ...]})
#   recruit(queries, refs) -> {query label: (ref label, pct id)}
#
# Everything engine-specific — the search itself, the greedy ordering, where the
# coverage rule is enforced — lives below this line.


class Runner:
    """Runs a bundled Linux binary: natively when it can, Docker `linux/amd64` when not.

    Both `bin/usearch11` and `bin/diamond` are Linux x86-64 ELF, so on this arm64 host
    they need a container. DIAMOND additionally ships as a normal macOS binary, so a
    `diamond` on PATH is preferred and skips Docker entirely (see `native`).
    """

    def __init__(self, binary: Path, workdir: Path, use_docker: bool, image: str,
                 native: str | None = None):
        self.binary, self.workdir, self.image = binary, workdir, image
        self.native = native                      # PATH binary; when set, no Docker
        self.use_docker = use_docker and not native

    def path(self, p: Path) -> str:
        """Rewrite a workdir path into the container's view of it."""
        rel = Path(p).resolve().relative_to(self.workdir.resolve())
        return f"/w/{rel}" if self.use_docker else str(Path(p).resolve())

    def describe(self) -> str:
        if self.native:
            return f"native {self.native}"
        return f"docker {self.image} {self.binary.name}" if self.use_docker else str(self.binary)

    def run(self, args: list[str], log_name: str, check: bool = True):
        if self.use_docker:
            cmd = ["docker", "run", "--rm", "--platform", "linux/amd64",
                   "-v", f"{self.binary.parent.resolve()}:/b",
                   "-v", f"{self.workdir.resolve()}:/w",
                   self.image, f"/b/{self.binary.name}"] + args
        else:
            cmd = [self.native or str(self.binary.resolve())] + args
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if log_name:
            (self.workdir / log_name).write_text(proc.stdout + proc.stderr)
        if check and proc.returncode != 0:
            hint = ""
            if self.use_docker and "docker" in (proc.stderr or "").lower():
                hint = ("\nhint: Docker Desktop must be running (the bundled binaries are "
                        "Linux x86-64 ELF and cannot run natively on arm64 macOS).")
            where = f"; see {self.workdir / log_name}" if log_name else ""
            sys.exit(f"error: {self.binary.name} step failed (exit {proc.returncode})"
                     f"{where}\n{proc.stderr.strip()}{hint}")
        return proc


# ------------------------------------------------------------------ usearch ---

def parse_uc(path: Path):
    """Parse a USEARCH .uc record file.

    Returns (centroids, hits) where centroids maps cluster number -> centroid label
    and hits is a list of (query, target, pct_id). `N` (no hit) lines are ignored;
    callers derive the unmatched set from the input labels instead.
    """
    centroids: dict[int, str] = {}
    hits: list[tuple[str, str, str]] = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        if f[0] == "S":
            centroids[int(f[1])] = f[8]
        elif f[0] == "H":
            hits.append((f[8], f[9], f[3]))
    return centroids, hits


class UsearchEngine:
    """`cluster_fast` per phase, `usearch_global` for recruitment. No shared substrate:
    every phase re-searches, which is what USEARCH offers."""

    name = "usearch"

    def __init__(self, runner: Runner, work: Path, idp: str, member_cov: float | None):
        self.r, self.work, self.idp = runner, work, idp
        # One coverage rule, applied identically in every phase — mixing rules between the
        # recruitment and de-novo stages shifts boundary sequences between clusters on
        # coverage definition alone (see reference-seed-clustering.md, "Coverage semantics").
        self.cov = ["-query_cov", f"{member_cov:g}"] if member_cov else []

    def prepare(self, records: list[dict]) -> None:
        pass                                     # nothing to build up front

    def version(self) -> str:
        """USEARCH has no `version` subcommand; it banners itself on every run."""
        for log in sorted(self.work.glob("*.log")):
            head = log.read_text()[:200].strip().splitlines()
            if head and head[0].startswith("usearch"):
                return head[0].split(",")[0]     # drop the host's RAM/core chatter
        return "usearch (version unknown)"

    def cluster(self, records: list[dict], tag: str):
        fa = self.work / f"{tag}.fasta"
        write_fasta(fa, records)
        self.r.run(["-cluster_fast", self.r.path(fa), "-id", self.idp, "-sort", "length"]
                   + self.cov +
                   ["-centroids", self.r.path(self.work / f"{tag}-centroids.fasta"),
                    "-uc", self.r.path(self.work / f"{tag}.uc")], f"{tag}.log")
        centroids, hits = parse_uc(self.work / f"{tag}.uc")
        members: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for q, t, pid in hits:
            members[t].append((q, pid))
        return [centroids[n] for n in sorted(centroids)], members

    def recruit(self, queries: list[dict], refs: list[dict]):
        qfa, dbfa = self.work / "02-queries.fasta", self.work / "02-refdb.fasta"
        write_fasta(qfa, queries)
        write_fasta(dbfa, refs)
        self.r.run(["-usearch_global", self.r.path(qfa), "-db", self.r.path(dbfa),
                    "-id", self.idp, "-maxaccepts", "0", "-maxrejects", "0", "-top_hit_only"]
                   + self.cov +
                   ["-uc", self.r.path(self.work / "02-closedref.uc"),
                    "-notmatched", self.r.path(self.work / "02-notmatched.fasta")],
                   "02-closedref.log")
        _, hits = parse_uc(self.work / "02-closedref.uc")
        return {q: (t, pid) for q, t, pid in hits}


# ------------------------------------------------------------------ diamond ---

class DiamondEngine:
    """One all-vs-all edge graph, built once at `--id`, reused by all three phases.

    Node ids handed to DIAMOND are ORDINALS (`--max-oid`), not labels: the label format
    is `U0001|WP_012345.1`, and `|` is exactly the character DIAMOND's seqid parser
    splits accessions on. Ordinals sidestep that parser completely and keep the mapping
    back to labels in this file, where it is checkable.

    Phase 1/3 are `greedy-vertex-cover` over the induced subgraph (edges reindexed to the
    subset, so a phase never sees a node it is not clustering). Phase 2 walks the same
    edges in Python: DIAMOND has no closed-reference mode, and the rule is one line —
    recruit `v` to the seed that covers it best.
    """

    name = "diamond"

    def __init__(self, runner: Runner, work: Path, ident: float, member_cov: float | None,
                 threads: int):
        self.r, self.work, self.threads = runner, work, threads
        self.pid_min = ident * 100.0
        self.cov_min = (member_cov or 0.0) * 100.0
        self.labels: list[str] = []              # oid -> label
        self.oid: dict[str, int] = {}
        # (q, s) -> (pident, qcovhsp, scovhsp, bits); q covered qcov% by s and vice versa
        self.edges: dict[tuple[int, int], tuple[float, float, float, float]] = {}

    def version(self) -> str:
        out = self.r.run(["version"], "", check=False).stdout.strip().splitlines()
        return out[0] if out else "diamond (version unknown)"

    # -- Phase 0: the substrate ------------------------------------------------
    def prepare(self, records: list[dict]) -> None:
        self.labels = [r["label"] for r in records]
        self.oid = {lab: i for i, lab in enumerate(self.labels)}
        faa, db = self.work / "00-substrate.faa", self.work / "00-substrate"
        with open(faa, "w") as fh:
            for i, rec in enumerate(records):
                fh.write(f">n{i}\n")
                for j in range(0, len(rec["seq"]), 60):
                    fh.write(rec["seq"][j:j + 60] + "\n")
        self.r.run(["makedb", "--in", self.r.path(faa), "--db", self.r.path(db),
                    "--threads", str(self.threads)], "00-makedb.log")
        out = self.work / "00-edges.tsv"
        # Permissive on e-value, strict on identity and coverage: the edge graph IS the
        # 90% criterion, so every downstream phase inherits it and cannot drift from it.
        cover = ["--query-or-subject-cover", f"{self.cov_min:g}"] if self.cov_min else []
        self.r.run(["blastp", "--query", self.r.path(faa), "--db", self.r.path(db.with_suffix(".dmnd")),
                    "--out", self.r.path(out),
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "bitscore",
                    "--id", f"{self.pid_min:g}"] + cover +
                   ["--ultra-sensitive", "--max-target-seqs", "0", "--evalue", "10",
                    "--masking", "0", "--comp-based-stats", "0",
                    "--threads", str(self.threads)], "00-blastp.log")
        for line in open(out):
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                q, s = int(f[0][1:]), int(f[1][1:])
                val = (float(f[2]), float(f[3]), float(f[4]), float(f[5]))
            except ValueError:
                continue
            if q == s:
                continue
            cur = self.edges.get((q, s))
            if cur is None or val[3] > cur[3]:   # best HSP per ordered pair, by bits
                self.edges[(q, s)] = val
        print(f"substrate: {len(records)} sequences, {len(self.edges)} directed edges "
              f"at >={self.pid_min:g}% id"
              + (f" and >={self.cov_min:g}% query-or-subject coverage" if self.cov_min else ""))

    def _cover_of(self, member: int, rep: int):
        """How well `rep` covers `member`, over both stored directions of the pair.

        Returns (cover%, bits, pident) or None. An edge is stored per direction, and
        the member's coverage is `qcovhsp` when it is the query and `scovhsp` when it
        is the subject — the same quantity read off either row.
        """
        best = None
        for key, cov_i in (((member, rep), 1), ((rep, member), 2)):
            e = self.edges.get(key)
            if e is None:
                continue
            cand = (e[cov_i], e[3], e[0])
            if best is None or cand > best:
                best = cand
        return best

    # -- Phases 1 and 3: greedy vertex cover -----------------------------------
    def cluster(self, records: list[dict], tag: str):
        sub = [self.oid[r["label"]] for r in records]
        local = {o: i for i, o in enumerate(sub)}
        efile, ofile = self.work / f"{tag}-edges.tsv", self.work / f"{tag}-gvc.tsv"
        with open(efile, "w") as fh:
            for (q, s), (pid, qcov, scov, bits) in self.edges.items():
                if q in local and s in local:
                    fh.write(f"{local[q]}\t{local[s]}\t{qcov:g}\t{scov:g}\t{bits:g}\n")
        if not sub:
            return [], defaultdict(list)
        self.r.run(["greedy-vertex-cover",
                    "--db", self.r.path(self.work / "00-substrate.dmnd"),
                    "--edges", self.r.path(efile), "--out", self.r.path(ofile),
                    "--member-cover", f"{self.cov_min:g}",
                    "--max-oid", str(len(sub) - 1),
                    "--threads", str(self.threads)], f"{tag}.log")
        members: dict[str, list[tuple[str, str]]] = defaultdict(list)
        centroid_of: dict[int, int] = {}
        for line in open(ofile):
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            centroid_of[sub[int(f[1])]] = sub[int(f[0])]
        # Centroids in input order: `greedy-vertex-cover` emits by internal degree order,
        # which is not stable to report against. Input order is the union's own order.
        order = [o for o in sub if centroid_of.get(o) == o]
        for o in sub:
            c = centroid_of.get(o, o)
            if c == o:
                continue
            hit = self._cover_of(o, c)
            members[self.labels[c]].append((self.labels[o], f"{hit[2]:.1f}" if hit else "-"))
        missing = [o for o in sub if o not in centroid_of]
        if missing:                              # GVC emits every node; belt and braces
            sys.exit(f"error: greedy-vertex-cover left {len(missing)} node(s) unassigned in {tag}")
        return [self.labels[o] for o in order], members

    # -- Phase 2: closed-reference recruitment ---------------------------------
    def recruit(self, queries: list[dict], refs: list[dict]):
        ref_oids = [self.oid[r["label"]] for r in refs]
        out = {}
        for q in queries:
            qo = self.oid[q["label"]]
            best, best_ref = None, None
            for ro in ref_oids:
                hit = self._cover_of(qo, ro)     # (cover, bits, pident)
                if hit is None or hit[0] < self.cov_min:
                    continue
                if best is None or hit > best:   # contention: cover, then bits, then id
                    best, best_ref = hit, ro
            if best_ref is not None:
                out[q["label"]] = (self.labels[best_ref], f"{best[2]:.1f}")
        return out


# ------------------------------------------------------------------------- main

def parse_args():
    p = argparse.ArgumentParser(
        description="Reference-seeded clustering of a union table; one row per cluster.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("union", type=Path, help="union TSV (enzyme_name, accession, pazy_id, aa_sequence, source)")
    p.add_argument("seeds", type=Path, help="curated seed TSV (enzyme_name, accession, pazy_id, justification, aa_sequence; optional plasticome_id)")
    p.add_argument("-o", "--out", type=Path, default=Path("clusters.tsv"), help="output TSV (default: ./clusters.tsv)")
    p.add_argument("--id", type=float, default=0.90, help="identity threshold (default: 0.90)")
    p.add_argument("--engine", choices=["diamond", "usearch"], default="diamond",
                   help="clustering engine (default: diamond)")
    p.add_argument("--member-cov", type=float, default=None,
                   help="member coverage rule: a sequence only joins a cluster if it is at "
                        "least this fraction covered by its representative. Defaults per "
                        "engine — 0.90 for diamond (whose identity is over a local HSP, so "
                        "coverage is what makes the threshold mean anything), unset for "
                        "usearch (identity alone, as in the canonical 90%% run). Applied "
                        "identically in all three phases; 0 disables it")
    p.add_argument("--workdir", type=Path, default=None,
                   help="intermediate FASTA/edges/logs "
                        "(default: <out stem>.intermediates/ beside the output)")
    p.add_argument("--centroids-out", type=Path, default=None, help="also write all representatives as FASTA")
    p.add_argument("--on-seed-collapse", choices=["keep", "drop"], default="keep",
                   help="seeds redundant with another seed at --id: keep both as designated "
                        "centroids (default) or drop the redundant one into the survivor's cluster")
    p.add_argument("--engine-bin", type=Path, default=None,
                   help="engine binary (default: <repo>/bin/diamond or <repo>/bin/usearch11)")
    p.add_argument("--threads", type=int, default=4, help="diamond threads (default: 4)")
    p.add_argument("--no-docker", action="store_true",
                   help="run the bundled Linux binary natively (Linux x86-64 host). A "
                        "`diamond` on PATH is used natively regardless")
    p.add_argument("--docker-image", default="debian:stable-slim", help="image used to run the engine (default: debian:stable-slim)")
    return p.parse_args()


def build_engine(args, work: Path):
    """Resolve binary, container policy and coverage default for the selected engine."""
    repo = Path(__file__).resolve().parents[2]
    if args.member_cov is None:
        args.member_cov = 0.90 if args.engine == "diamond" else None

    if args.engine == "diamond":
        binary = args.engine_bin or repo / "bin" / "diamond"
        # The bundled copy is a Linux ELF; a `diamond` on PATH is a real macOS binary and
        # is both faster and Docker-free, so it wins unless a binary was named explicitly.
        native = None if args.engine_bin else shutil.which("diamond")
        if not native and not binary.exists():
            sys.exit(f"error: no diamond on PATH and no binary at {binary} (pass --engine-bin)")
        runner = Runner(binary, work, not args.no_docker, args.docker_image, native)
        engine = DiamondEngine(runner, work, args.id, args.member_cov, args.threads)
    else:
        binary = args.engine_bin or repo / "bin" / "usearch11"
        if not binary.exists():
            sys.exit(f"error: usearch binary not found: {binary} (pass --engine-bin)")
        runner = Runner(binary, work, not args.no_docker, args.docker_image)
        engine = UsearchEngine(runner, work, f"{args.id:g}", args.member_cov)

    if args.no_docker and not runner.native and (
            platform.system() != "Linux" or platform.machine() not in ("x86_64", "amd64")):
        print(f"warning: --no-docker on {platform.system()}/{platform.machine()}; "
              f"the bundled binaries are Linux x86-64 ELF")
    return engine, runner


def load_records(args):
    """Read both tables into labelled records; return (refs, remaining, dropped, n_union).

    `refs` is one record per curated seed, in seed order — the reference centroids.
    Each is the UNION ROW carrying that seed's exact sequence, promoted in place, so
    the seed contributes an annotation rather than a second copy of a sequence the
    union already holds. `n_union` is the number of union rows with a sequence, which
    is what the universe must come to once nothing is double-counted.
    """
    seed_rows = read_table(args.seeds, SEED_COLS, "seed", SEED_OPTIONAL_COLS)
    union_rows = read_table(args.union, UNION_COLS, "union")

    seeds, empty_seed = [], 0
    for i, r in enumerate(seed_rows, 1):
        seq = normalize(r["aa_sequence"])
        if not seq:
            empty_seed += 1
            continue
        seeds.append({
            "label": f"S{i:03d}|{safe(r['accession'] or r['plasticome_id'])}",
            "enzyme_name": r["enzyme_name"], "accession": r["accession"],
            "pazy_id": r["pazy_id"], "plasticome_id": r["plasticome_id"],
            "source": "seed", "seq": seq,
        })
    if empty_seed:
        print(f"warning: {empty_seed} seed row(s) had no sequence and were skipped")

    union, empty_union = [], 0
    for i, r in enumerate(union_rows, 1):
        seq = normalize(r["aa_sequence"])
        if not seq:
            empty_union += 1
            continue
        union.append({
            "label": f"U{i:04d}|{safe(r['accession'])}",
            "enzyme_name": r["enzyme_name"], "accession": r["accession"],
            "pazy_id": r["pazy_id"], "plasticome_id": "",
            "source": r["source"], "seq": seq,
        })
    if empty_union:
        print(f"warning: {empty_union} union row(s) had no sequence and were skipped")

    # Promote, in seed order, the union row carrying each seed's exact sequence: that
    # row *is* the seed, so it becomes the reference centroid instead of the seed being
    # added beside it. Claims are exclusive, so two seeds sharing a sequence take two
    # different union rows rather than both landing on the first one.
    by_seq_union: dict[str, list[dict]] = defaultdict(list)
    for rec in union:
        by_seq_union[rec["seq"]].append(rec)
    claimed: set[int] = set()
    refs, seed_only = [], []
    for s in seeds:
        promoted = next((r for r in by_seq_union.get(s["seq"], ()) if id(r) not in claimed), None)
        if promoted is None:
            refs.append(s)              # not in the union: the seed is its own record
            seed_only.append(s)
            continue
        claimed.add(id(promoted))
        # Curated values win; the union row fills whatever the seed leaves blank, so
        # neither table's annotation is lost in the merge.
        for col in ("enzyme_name", "accession", "pazy_id", "plasticome_id"):
            promoted[col] = s[col] or promoted[col]
        promoted["source"] = "seed"
        refs.append(promoted)
    if seed_only:
        print(f"warning: {len(seed_only)} seed(s) have no exact-sequence union row "
              f"({', '.join(s['plasticome_id'] or s['label'] for s in seed_only)}); each "
              f"stays its own record, so the universe exceeds the union by that many")

    by_acc = {acc_key(r["accession"]): r["label"] for r in reversed(refs) if acc_key(r["accession"])}
    by_seq = {r["seq"]: r["label"] for r in reversed(refs)}

    remaining, dropped = [], []
    for rec in union:
        if id(rec) in claimed:
            continue                    # already a reference centroid
        # A further union row that is also a seed leaves REMAINING so the pool stays
        # disjoint from the seed centroids, but it is still a record of the universe:
        # it is folded into that seed's cluster as a member rather than dropped.
        anchor = by_acc.get(acc_key(rec["accession"])) or by_seq.get(rec["seq"])
        if anchor:
            dropped.append((rec, anchor))
        else:
            remaining.append(rec)
    return refs, remaining, dropped, len(union)


def main():
    args = parse_args()

    # Sidecar for the intermediates, named for the deliverable it produced:
    # `02-clusters.tsv` -> `02-clusters.intermediates/`. The suffix is *replaced*,
    # not appended (the older scheme was `<out>.work/`), so the folder reads as a
    # peer of the output rather than as a second extension on it.
    work = args.workdir or args.out.with_suffix(".intermediates")
    work.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    engine, runner = build_engine(args, work)

    refs, remaining, dropped, n_union = load_records(args)
    if not refs:
        sys.exit("error: no usable seed sequences")
    by_label = {r["label"]: r for r in refs + remaining + [d for d, _ in dropped]}
    print(f"universe: {len(refs)} seed centroids + {len(remaining)} remaining "
          f"+ {len(dropped)} union row(s) folded into their seed = {len(by_label)} records")
    extra = len(by_label) - n_union
    print(f"          {n_union} union row(s) with a sequence"
          + (f" + {extra} seed(s) absent from the union"
             if extra else "; no record is counted twice"))

    write_fasta(work / "seeds.fasta", refs)
    write_fasta(work / "remaining.fasta", remaining)
    idp = f"{args.id:g}"
    cov_txt = (f"{args.member_cov:g} (member covered by its representative)"
               if args.member_cov else "not enforced (identity only)")
    print(f"engine: {engine.name} ({runner.describe()})   id={idp}   member coverage: {cov_txt}")

    # --- Phase 0: the shared substrate (diamond only) ----------------------------
    engine.prepare(refs + remaining)

    # --- Phase 1: are the seeds non-redundant at $id? ----------------------------
    print(f"\n== Phase 1: seed self-cluster diagnostic @ {idp} ==")
    _, seed_members = engine.cluster(refs, "01-seed")
    collapsed = {q: (t, pid) for t, mem in seed_members.items() for q, pid in mem}
    if collapsed:
        print(f"{len(collapsed)} of {len(refs)} seed(s) are redundant at {idp}:")
        for q, (t, pid) in collapsed.items():
            print(f"   {q} -> {t} ({pid}% id)")
        print(f"policy --on-seed-collapse={args.on_seed_collapse}: "
              + ("all seeds kept as designated centroids (near-duplicate reps accepted)"
                 if args.on_seed_collapse == "keep" else
                 "redundant seeds demoted to members of the surviving seed's cluster"))
    else:
        print(f"all {len(refs)} seeds are non-redundant at {idp}; each is a reference centroid")

    if args.on_seed_collapse == "keep":
        demoted = {}
        seed_db = refs
    else:
        demoted = collapsed
        seed_db = [r for r in refs if r["label"] not in demoted]
    ref_centroids = [r["label"] for r in seed_db]

    # --- Phase 2: closed-reference recruitment -----------------------------------
    print(f"\n== Phase 2: closed-reference recruitment of {len(remaining)} remaining ==")
    recruited = engine.recruit(remaining, seed_db)
    unmatched = [r for r in remaining if r["label"] not in recruited]
    print(f"recruited to a seed: {len(recruited)}   unmatched: {len(unmatched)}")

    # --- Phase 3: open-reference (de novo) ---------------------------------------
    print(f"\n== Phase 3: de-novo clustering of the {len(unmatched)} unmatched ==")
    # seeds/remaining/notmatched are the engine-independent record of what each phase was
    # handed, in a stable order; the engine writes whatever else it needs beside them.
    write_fasta(work / "notmatched.fasta", unmatched)
    dn_order, dn_members = engine.cluster(unmatched, "03-denovo")
    print(f"de-novo clusters: {len(dn_order)}")

    # --- Phase 4: merge & emit ----------------------------------------------------
    members: dict[str, list[tuple[str, str]]] = defaultdict(list)   # centroid -> [(label, pid)]
    for q, (t, pid) in recruited.items():
        members[t].append((q, pid))
    for q, (t, pid) in demoted.items():
        members[t].append((q, pid))
    for rec, anchor in dropped:                      # union rows that are a seed
        if anchor in demoted:                        # its seed was itself demoted
            anchor = demoted[anchor][0]
        members[anchor].append((rec["label"], "dup"))

    # Reference clusters are numbered first, in seed order, so cluster_id 1..N_seeds
    # is always the curated seeds and everything above them is de novo. That is the
    # only remaining marker of origin — the column was dropped — so the two groups
    # must stay in this order.
    clusters = []
    for lab in ref_centroids:
        clusters.append(("reference", lab, members.get(lab, [])))
    for lab in dn_order:
        clusters.append(("de_novo", lab, dn_members.get(lab, [])))

    rows = []
    for i, (origin, lab, mem) in enumerate(clusters, start=1):
        rep = by_label[lab]
        mem = sorted(mem, key=lambda m: (-pct(m[1]), m[0]))   # closest to the centroid first
        rows.append({
            "cluster_id": i, "_origin": origin, "size": len(mem) + 1,
            "rep_label": lab, "rep_enzyme_name": rep["enzyme_name"],
            "rep_accession": rep["accession"], "rep_pazy_id": rep["pazy_id"],
            "rep_plasticome_id": rep["plasticome_id"], "rep_source": rep["source"],
            "rep_seq_len": len(rep["seq"]),
            "member_labels": JOIN.join(m for m, _ in mem),
            "member_enzyme_names": JOIN.join(by_label[m]["enzyme_name"] for m, _ in mem),
            "member_accessions": JOIN.join(by_label[m]["accession"] for m, _ in mem),
            "member_pct_ids": JOIN.join(pid for _, pid in mem),
            "rep_aa_sequence": rep["seq"],
        })
    rows.sort(key=lambda r: r["cluster_id"])   # earliest to latest; seeds first

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    if args.centroids_out:
        args.centroids_out.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(args.centroids_out, [by_label[r["rep_label"]] for r in rows])

    # --- report (reference-seed-clustering.md, "Reporting checklist") -------------
    placed = sum(r["size"] for r in rows)
    sizes = Counter(r["size"] for r in rows)
    singl = sum(1 for r in rows if r["size"] == 1)
    ref_rows = [r for r in rows if r["_origin"] == "reference"]
    big = max(rows, key=lambda r: r["size"])
    print(f"\n== Phase 4: {args.out} ==")
    print(f"clusters              : {len(rows)}  (reference {len(ref_rows)} + de_novo {len(rows) - len(ref_rows)})")
    print(f"records placed        : {placed} / {len(by_label)}"
          + ("  ** MISMATCH **" if placed != len(by_label) else "  (all accounted for)"))
    print(f"seeds as centroids    : {len(ref_centroids)} / {len(refs)}")
    print(f"recruited to seeds    : {len(recruited)} (+{len(dropped)} folded seed duplicates)"
          f"   de-novo members: {len(unmatched)}")
    print(f"singletons            : {singl} ({100 * singl / len(rows):.1f}% of clusters)")
    print(f"largest cluster       : {big['size']} -> {big['cluster_id']} {big['rep_label']}")
    print(f"size distribution     : {dict(sorted(sizes.items()))}")
    print(f"engine                : {engine.version()} ({runner.describe()})")
    print(f"work dir              : {work}")

    # Engine + parameter provenance, item 3 of the doc's reporting checklist. Written
    # beside the intermediates so summarize_run.py can state which engine produced the
    # table without re-deriving it from logs.
    (work / "provenance.json").write_text(json.dumps({
        "engine": engine.name, "engine_version": engine.version(),
        "invocation": runner.describe(),
        "id": args.id, "member_cov": args.member_cov,
        "greedy_ordering": "graph degree (greedy-vertex-cover)" if engine.name == "diamond"
                           else "length (cluster_fast -sort length)",
        "on_seed_collapse": args.on_seed_collapse,
        "n_records": len(by_label), "n_clusters": len(rows),
        "n_reference_clusters": len(ref_rows), "n_recruited": len(recruited),
        "n_singletons": singl,
    }, indent=2))

    if placed != len(by_label):
        sys.exit("error: cluster membership does not cover the universe")


if __name__ == "__main__":
    main()
