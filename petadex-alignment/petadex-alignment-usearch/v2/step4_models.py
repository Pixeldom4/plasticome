#!/usr/bin/env python3
"""Step 4 - build and validate one domain model (HMM) per component.

Per component: MSA (mafft --auto) -> hmmbuild -> hmmsearch back against the
members, mirroring the paper's HMM-coverage filtering. Emits a manifest with
member count, mean coverage, and model_source.

A component whose HMM covers little of its members is not a domain family and
cannot seed an expansion - those are flagged rather than silently shipped.
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

MAFFT = os.environ.get("MAFFT", "/watson/rnalab/js/programs/mafft")
COVERAGE_FLOOR = 0.50  # mean member coverage below this -> flagged

SOURCE = {"inherited": "inherited", "merged": "merge-rebuilt",
          "novel": "novel-built", "split": "split-rebuilt"}


def run(cmd, **kw):
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        sys.exit(f"FAILED: {' '.join(map(str, cmd))}\n{p.stderr[-2000:]}")
    return p


def read_fasta(path):
    d, k, buf = {}, None, []
    for line in Path(path).open():
        if line.startswith(">"):
            if k:
                d[k] = "".join(buf)
            k, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if k:
        d[k] = "".join(buf)
    return d


def domtbl_coverage(path):
    """best-domain envelope coverage per target sequence, plus hmm coverage"""
    best = {}
    for line in Path(path).open():
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 23:
            continue
        tgt, tlen, qlen = f[0], int(f[2]), int(f[5])
        score = float(f[13])
        hmm_from, hmm_to = int(f[15]), int(f[16])
        env_from, env_to = int(f[19]), int(f[20])
        if tgt not in best or score > best[tgt][0]:
            best[tgt] = (score, (env_to - env_from + 1) / tlen,
                         (hmm_to - hmm_from + 1) / qlen)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignment", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    if not (shutil.which(MAFFT) or Path(MAFFT).exists()):
        sys.exit(f"mafft not found at {MAFFT}")

    seqs = read_fasta(args.fasta)
    rows = list(csv.DictReader(Path(args.assignment).open()))
    comps = defaultdict(list)
    kind, merged_from = {}, {}
    for r in rows:
        comps[r["component_id"]].append(r)
        kind[r["component_id"]] = r["component_kind"]
        merged_from[r["component_id"]] = r["merged_from"]

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "_work"
    work.mkdir(exist_ok=True)

    manifest = []
    for cid in sorted(comps, key=lambda c: (-len(comps[c]), c)):
        members = comps[cid]
        mf = work / f"{cid}.fa"
        with mf.open("w") as fh:
            for m in members:
                fh.write(f">{m['node_id']}\n{seqs[m['node_id']]}\n")

        afa = work / f"{cid}.afa"
        if len(members) == 1:
            afa.write_text(mf.read_text())
        else:
            p = run([MAFFT, "--auto", "--anysymbol", "--quiet", str(mf)])
            afa.write_text(p.stdout)

        hmm = out / f"{cid}.hmm"
        run(["hmmbuild", "--amino", "-n", f"PETadex_{cid}", str(hmm), str(afa)])

        dom = work / f"{cid}.domtbl"
        run(["hmmsearch", "--max", "--domtblout", str(dom), str(hmm), str(mf)])
        cov = domtbl_coverage(dom)

        n = len(members)
        detected = len(cov)
        seq_cov = [v[1] for v in cov.values()]
        hmm_cov = [v[2] for v in cov.values()]
        mean_seq = sum(seq_cov) / len(seq_cov) if seq_cov else 0.0
        mean_hmm = sum(hmm_cov) / len(hmm_cov) if hmm_cov else 0.0
        hmm_len = int([l.split()[1] for l in hmm.read_text().splitlines()
                       if l.startswith("LENG")][0])

        flags = []
        if detected < n:
            flags.append(f"undetected={n - detected}")
        if mean_seq < COVERAGE_FLOOR:
            flags.append("low_member_coverage")
        if n == 1:
            flags.append("singleton_model")

        manifest.append({
            "component_id": cid, "model_source": SOURCE[kind[cid]],
            "n_members": n, "n_v1": sum(1 for m in members if m["is_new"] == "0"),
            "n_new": sum(1 for m in members if m["is_new"] == "1"),
            "merged_from": merged_from[cid], "hmm_length": hmm_len,
            "n_detected": detected,
            "mean_member_coverage": round(mean_seq, 4),
            "mean_hmm_coverage": round(mean_hmm, 4),
            "min_member_coverage": round(min(seq_cov), 4) if seq_cov else 0.0,
            "flags": ";".join(flags),
        })
        print(f"{cid:6s} {SOURCE[kind[cid]]:14s} n={n:3d} hmm_len={hmm_len:4d} "
              f"cov={mean_seq:.3f} {';'.join(flags)}", file=sys.stderr)

    mpath = out / "manifest.csv"
    with mpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0]))
        w.writeheader()
        w.writerows(manifest)
    print(f"\nwrote {mpath} ({len(manifest)} models)", file=sys.stderr)


if __name__ == "__main__":
    main()
