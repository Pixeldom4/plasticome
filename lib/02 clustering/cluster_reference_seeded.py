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

Engine: USEARCH v11 `cluster_fast` / `usearch_global`, `-sort length`. The bundled
`bin/usearch11` is a Linux x86-64 ELF, so on arm64 macOS it runs inside a Docker
`linux/amd64` container (Docker Desktop must be running); pass `--no-docker` on a
Linux x86-64 host.

`usearch_global` recruits on **identity alone**, which is what the canonical 90% run did.
`--member-cov` adds the doc's member-coverage rule (`-query_cov`, applied identically in
all three phases) and recruits fewer sequences to the seeds — on the v1/v260701 union at
90%, 17 instead of 28.

Inputs
------
UNION_TSV   union table  : enzyme_name, accession, pazy_id, aa_sequence, source
SEEDS_TSV   curated seeds: plasticome_id, enzyme_name, accession, pazy_id,
                           justification, aa_sequence
Extra columns in either file are ignored.

`SEEDS` and `REMAINING` must be disjoint: a union row that *is* a seed (same accession,
version-insensitive, or an identical sequence) leaves the pool and is folded into that
seed's cluster as a member with `dup` in place of a percent identity — so every input
row is still accounted for in the output.

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

Labels are `S###|accession` (seeds) and `U####|accession` (union rows); the numeric
prefix keeps them unique when accessions repeat or are blank.

Examples
--------
  python3 cluster_reference_seeded.py union.tsv seeds.tsv -o clusters.tsv

  python3 cluster_reference_seeded.py union.tsv seeds.tsv -o clusters.tsv \\
      --id 0.95 --centroids-out all-centroids.fasta --workdir runs/<run>/clustering/work
"""

from __future__ import annotations

import argparse
import csv
import platform
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

UNION_COLS = ["enzyme_name", "accession", "pazy_id", "aa_sequence", "source"]
SEED_COLS = ["plasticome_id", "enzyme_name", "accession", "pazy_id", "aa_sequence"]

OUT_COLS = [
    "cluster_id", "size",
    "rep_label", "rep_enzyme_name", "rep_accession", "rep_pazy_id",
    "rep_plasticome_id", "rep_source", "rep_seq_len",
    "member_labels", "member_enzyme_names", "member_accessions", "member_pct_ids",
    "rep_aa_sequence",
]
JOIN = "; "


# --------------------------------------------------------------------------- io

def read_table(path: Path, required: list[str], what: str) -> list[dict]:
    """Read a TSV into dicts, keeping the FIRST column of any repeated header name.

    csv.DictReader keeps the *last*, which silently blanks every row when a sheet
    carries a stray empty duplicate (the v260701 sheet has done exactly that with a
    second `aa_sequence`), so the header is indexed by hand here.
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
    out = []
    for row in rows[1:]:
        if not any(cell.strip() for cell in row):
            continue
        out.append({c: (row[i].strip() if i < len(row) else "") for c, i in idx.items()})
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


# ---------------------------------------------------------------------- usearch

class Usearch:
    """USEARCH runner — native on Linux x86-64, Docker `linux/amd64` elsewhere."""

    def __init__(self, binary: Path, workdir: Path, use_docker: bool, image: str):
        self.binary, self.workdir, self.use_docker, self.image = binary, workdir, use_docker, image

    def path(self, p: Path) -> str:
        """Rewrite a workdir path into the container's view of it."""
        rel = Path(p).resolve().relative_to(self.workdir.resolve())
        return f"/w/{rel}" if self.use_docker else str(Path(p).resolve())

    def run(self, args: list[str], log_name: str) -> None:
        if self.use_docker:
            cmd = ["docker", "run", "--rm", "--platform", "linux/amd64",
                   "-v", f"{self.binary.parent.resolve()}:/b",
                   "-v", f"{self.workdir.resolve()}:/w",
                   self.image, f"/b/{self.binary.name}"] + args
        else:
            cmd = [str(self.binary.resolve())] + args
        proc = subprocess.run(cmd, capture_output=True, text=True)
        (self.workdir / log_name).write_text(proc.stdout + proc.stderr)
        if proc.returncode != 0:
            hint = ""
            if self.use_docker and "docker" in (proc.stderr or "").lower():
                hint = ("\nhint: Docker Desktop must be running (the bundled usearch binaries are "
                        "Linux x86-64 ELF and cannot run natively on arm64 macOS).")
            sys.exit(f"error: USEARCH step failed (exit {proc.returncode}); see "
                     f"{self.workdir / log_name}\n{proc.stderr.strip()}{hint}")


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


# ------------------------------------------------------------------------- main

def parse_args():
    p = argparse.ArgumentParser(
        description="Reference-seeded clustering of a union table; one row per cluster.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("union", type=Path, help="union TSV (enzyme_name, accession, pazy_id, aa_sequence, source)")
    p.add_argument("seeds", type=Path, help="curated seed TSV (plasticome_id, enzyme_name, accession, pazy_id, justification, aa_sequence)")
    p.add_argument("-o", "--out", type=Path, default=Path("clusters.tsv"), help="output TSV (default: ./clusters.tsv)")
    p.add_argument("--id", type=float, default=0.90, help="identity threshold (default: 0.90)")
    p.add_argument("--member-cov", type=float, default=None,
                   help="member coverage rule: a sequence only joins a cluster if it is at "
                        "least this fraction covered by its representative (usearch -query_cov, "
                        "applied identically in all three phases). Default: unset — usearch "
                        "recruits on identity alone, as in the canonical 90%% run")
    p.add_argument("--workdir", type=Path, default=None,
                   help="intermediate FASTA/.uc/logs "
                        "(default: <out stem>.intermediates/ beside the output)")
    p.add_argument("--centroids-out", type=Path, default=None, help="also write all representatives as FASTA")
    p.add_argument("--on-seed-collapse", choices=["keep", "drop"], default="keep",
                   help="seeds redundant with another seed at --id: keep both as designated "
                        "centroids (default) or drop the redundant one into the survivor's cluster")
    p.add_argument("--usearch-bin", type=Path, default=None, help="usearch binary (default: <repo>/bin/usearch11)")
    p.add_argument("--no-docker", action="store_true", help="run usearch natively (Linux x86-64 host)")
    p.add_argument("--docker-image", default="debian:stable-slim", help="image used to run usearch (default: debian:stable-slim)")
    return p.parse_args()


def load_records(args):
    """Read both tables into labelled records; return (seeds, remaining, dropped)."""
    seed_rows = read_table(args.seeds, SEED_COLS, "seed")
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

    by_acc = {acc_key(s["accession"]): s["label"] for s in reversed(seeds) if acc_key(s["accession"])}
    by_seq = {s["seq"]: s["label"] for s in reversed(seeds)}

    remaining, dropped, empty_union = [], [], 0
    for i, r in enumerate(union_rows, 1):
        seq = normalize(r["aa_sequence"])
        if not seq:
            empty_union += 1
            continue
        rec = {
            "label": f"U{i:04d}|{safe(r['accession'])}",
            "enzyme_name": r["enzyme_name"], "accession": r["accession"],
            "pazy_id": r["pazy_id"], "plasticome_id": "",
            "source": r["source"], "seq": seq,
        }
        # A union row that *is* a seed keeps SEEDS/REMAINING disjoint by leaving the
        # pool, but it is still a record of the universe: it is folded into that
        # seed's cluster as a member rather than dropped from the output.
        anchor = by_acc.get(acc_key(r["accession"])) or by_seq.get(seq)
        if anchor:
            dropped.append((rec, anchor))
        else:
            remaining.append(rec)
    if empty_union:
        print(f"warning: {empty_union} union row(s) had no sequence and were skipped")
    return seeds, remaining, dropped


def main():
    args = parse_args()
    repo = Path(__file__).resolve().parents[2]
    usearch_bin = args.usearch_bin or repo / "bin" / "usearch11"
    if not usearch_bin.exists():
        sys.exit(f"error: usearch binary not found: {usearch_bin} (pass --usearch-bin)")

    # Sidecar for the intermediates, named for the deliverable it produced:
    # `02-clusters.tsv` -> `02-clusters.intermediates/`. The suffix is *replaced*,
    # not appended (the older scheme was `<out>.work/`), so the folder reads as a
    # peer of the output rather than as a second extension on it.
    work = args.workdir or args.out.with_suffix(".intermediates")
    work.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    use_docker = not args.no_docker
    if args.no_docker and (platform.system() != "Linux" or platform.machine() not in ("x86_64", "amd64")):
        print(f"warning: --no-docker on {platform.system()}/{platform.machine()}; "
              f"the bundled usearch binaries are Linux x86-64 ELF")
    u = Usearch(usearch_bin, work, use_docker, args.docker_image)

    seeds, remaining, dropped = load_records(args)
    if not seeds:
        sys.exit("error: no usable seed sequences")
    by_label = {r["label"]: r for r in seeds + remaining + [d for d, _ in dropped]}
    print(f"universe: {len(seeds)} seeds + {len(remaining)} remaining "
          f"+ {len(dropped)} union row(s) folded into their seed = {len(by_label)} records")

    seeds_fa, remaining_fa = work / "seeds.fasta", work / "remaining.fasta"
    write_fasta(seeds_fa, seeds)
    write_fasta(remaining_fa, remaining)
    idp = f"{args.id:g}"
    # One coverage rule, applied identically in every phase — mixing rules between the
    # recruitment and de-novo stages shifts boundary sequences between clusters on
    # coverage definition alone (see reference-seed-clustering.md, "Coverage semantics").
    cov = ["-query_cov", f"{args.member_cov:g}"] if args.member_cov is not None else []
    print(f"engine: usearch ({'docker ' + args.docker_image if use_docker else 'native'})   "
          f"id={idp}   member coverage: {f'{args.member_cov:g}' if cov else 'not enforced (identity only)'}")

    # --- Phase 1: are the seeds non-redundant at $id? ----------------------------
    print(f"\n== Phase 1: seed self-cluster diagnostic @ {idp} ==")
    u.run(["-cluster_fast", u.path(seeds_fa), "-id", idp, "-sort", "length"] + cov +
          ["-centroids", u.path(work / "seed-centroids.fasta"),
           "-uc", u.path(work / "01-seed.uc")], "01-seed.log")
    _, seed_hits = parse_uc(work / "01-seed.uc")
    collapsed = {q: (t, pid) for q, t, pid in seed_hits}
    if collapsed:
        print(f"{len(collapsed)} of {len(seeds)} seed(s) are redundant at {idp}:")
        for q, (t, pid) in collapsed.items():
            print(f"   {q} -> {t} ({pid}% id)")
        print(f"policy --on-seed-collapse={args.on_seed_collapse}: "
              + ("all seeds kept as designated centroids (near-duplicate reps accepted)"
                 if args.on_seed_collapse == "keep" else
                 "redundant seeds demoted to members of the surviving seed's cluster"))
    else:
        print(f"all {len(seeds)} seeds are non-redundant at {idp}; each is a reference centroid")

    if args.on_seed_collapse == "keep":
        ref_centroids = [s["label"] for s in seeds]
        seed_db, demoted = seeds_fa, {}
    else:
        demoted = collapsed
        ref_centroids = [s["label"] for s in seeds if s["label"] not in demoted]
        seed_db = work / "seed-centroids.fasta"

    # --- Phase 2: closed-reference recruitment -----------------------------------
    print(f"\n== Phase 2: closed-reference recruitment of {len(remaining)} remaining ==")
    u.run(["-usearch_global", u.path(remaining_fa), "-db", u.path(seed_db), "-id", idp,
           "-maxaccepts", "0", "-maxrejects", "0", "-top_hit_only"] + cov +
          ["-uc", u.path(work / "02-closedref.uc"),
           "-notmatched", u.path(work / "notmatched.fasta")], "02-closedref.log")
    _, ref_hits = parse_uc(work / "02-closedref.uc")
    recruited = {q: (t, pid) for q, t, pid in ref_hits}
    unmatched = [r for r in remaining if r["label"] not in recruited]
    print(f"recruited to a seed: {len(recruited)}   unmatched: {len(unmatched)}")

    # --- Phase 3: open-reference (de novo) ---------------------------------------
    print(f"\n== Phase 3: de-novo clustering of the {len(unmatched)} unmatched ==")
    write_fasta(work / "notmatched.fasta", unmatched)   # authoritative, and stable order
    u.run(["-cluster_fast", u.path(work / "notmatched.fasta"), "-id", idp, "-sort", "length"] + cov +
          ["-centroids", u.path(work / "denovo-centroids.fasta"),
           "-uc", u.path(work / "03-denovo.uc")], "03-denovo.log")
    dn_centroids, dn_hits = parse_uc(work / "03-denovo.uc")
    print(f"de-novo clusters: {len(dn_centroids)}")

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
    dn_by_num = defaultdict(list)
    for q, t, pid in dn_hits:
        dn_by_num[t].append((q, pid))

    # Reference clusters are numbered first, in seed order, so cluster_id 1..N_seeds
    # is always the curated seeds and everything above them is de novo. That is the
    # only remaining marker of origin — the column was dropped — so the two groups
    # must stay in this order.
    clusters = []
    for lab in ref_centroids:
        clusters.append(("reference", lab, members.get(lab, [])))
    for num in sorted(dn_centroids):
        lab = dn_centroids[num]
        clusters.append(("de_novo", lab, dn_by_num.get(lab, [])))

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
    print(f"seeds as centroids    : {len(ref_centroids)} / {len(seeds)}")
    print(f"recruited to seeds    : {len(recruited)} (+{len(dropped)} folded seed duplicates)"
          f"   de-novo members: {len(unmatched)}")
    print(f"singletons            : {singl} ({100 * singl / len(rows):.1f}% of clusters)")
    print(f"largest cluster       : {big['size']} -> {big['cluster_id']} {big['rep_label']}")
    print(f"size distribution     : {dict(sorted(sizes.items()))}")
    print(f"work dir              : {work}")
    if placed != len(by_label):
        sys.exit("error: cluster membership does not cover the universe")


if __name__ == "__main__":
    main()
