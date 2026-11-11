#!/usr/bin/env python3
"""
Step 4: check which components can share a tree at all.

Grafting two component trees under one backbone claims they have a common
ancestor. Step 1 proposed fold groups from the CATH superfamily of the
catalytic domain. This step tests that proposal three independent ways, every
component against every other:

  TM-score   global structural superposition of the representative models
             (Foldseek in TM-align mode), normalised by the shorter chain.
             0.5 and above is the usual same-fold line.
  Foldseek   3Di plus amino acid local alignment E-value. Sensitive to
             homology, not just to overall shape.
  HHalign    HMM against HMM probability, profiles built from each
             component's own alignment. No structures involved.

TM-score alone is not enough. Alpha/beta hydrolases, SGNH hydrolases and
subtilisins are all three layer alpha/beta/alpha sandwiches (CATH 3.40.50), and
they reach TM 0.5 to 0.63 against each other without being homologous. So the
rule has two halves:

  inside a group   every component must reach TM >= 0.5 with another member,
                   or HHalign probability >= 90 when it has no structure
  between groups   a pair is flagged when Foldseek E <= 1e-3 or HHalign
                   probability >= 90, which would mean the groups should merge

Writes work/homology/pairs.tsv (every pair, all three scores),
work/homology/groups.tsv (per component verdict) and
out/homology_check.tsv (the group by group summary for the methods).
"""

import collections
import itertools
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

from common import OUT, THREADS, WORK, cdir, read_tsv, run, tool, write_tsv

HDIR = os.path.join(WORK, "homology")
TM_MIN = 0.5
FS_E = 1e-3
HH_PROB = 90.0
FMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits,qtmscore,ttmscore,prob"


def foldseek(structs, out, tmalign):
    tmp = os.path.join(HDIR, "tmp")
    cmd = [tool("foldseek"), "easy-search", structs, structs, out, tmp,
           "--exhaustive-search", "1", "-e", "10", "--threads", str(THREADS),
           "--format-output", FMT, "-v", "1"]
    if tmalign:
        cmd += ["--alignment-type", "1", "--tmscore-threshold", "0"]
    run(cmd, os.path.join(HDIR, "log.txt"))
    shutil.rmtree(tmp, ignore_errors=True)


def load_m8(path):
    """Best TM (normalised by the shorter chain) and best E-value per pair."""
    tm, ev = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, b = int(f[0][1:3]), int(f[1][1:3])
            if a == b:
                continue
            k = (min(a, b), max(a, b))
            tm[k] = max(tm.get(k, 0.0), float(f[8]), float(f[9]))
            ev[k] = min(ev.get(k, 10.0), float(f[6]))
    return tm, ev


def hhalign_all(comps):
    """HHalign probability for every pair of component profiles."""
    hdir = os.path.join(HDIR, "hhm")
    os.makedirs(hdir, exist_ok=True)
    devnull = subprocess.DEVNULL
    for c in comps:
        d = cdir(c)
        src = os.path.join(d, "reference.afa")
        if not os.path.exists(src):
            src = os.path.join(d, "seqs.faa")
        subprocess.run([tool("hhmake"), "-i", src, "-o",
                        os.path.join(hdir, f"c{c:02d}.hhm"), "-M", "50",
                        "-name", f"c{c:02d}", "-v", "0"],
                       check=True, stdout=devnull, stderr=devnull)

    def one(pair):
        a, b = pair
        r = subprocess.run([tool("hhalign"), "-i", os.path.join(hdir, f"c{a:02d}.hhm"),
                            "-t", os.path.join(hdir, f"c{b:02d}.hhm"),
                            "-o", "/dev/stdout", "-v", "0"],
                           capture_output=True, text=True, check=True)
        m = re.search(r"^\s*1\s.{30}\s+([\d.]+)\s", r.stdout, re.M)
        return pair, float(m.group(1)) if m else 0.0

    with ThreadPoolExecutor(THREADS) as ex:
        return dict(ex.map(one, itertools.combinations(sorted(comps), 2)))


def main():
    os.makedirs(HDIR, exist_ok=True)
    comps = {int(r["component"]): r for r in read_tsv(os.path.join(WORK, "components.tsv"))}
    has_struct = {int(r["component"]) for r in
                  read_tsv(os.path.join(WORK, "structures.tsv")) if r["source"] != "none"}
    group = {c: int(r["fold_group"]) for c, r in comps.items()}

    structs = os.path.join(HDIR, "structures")
    shutil.rmtree(structs, ignore_errors=True)
    os.makedirs(structs)
    for c in has_struct:
        shutil.copy(os.path.join(WORK, "structures", f"c{c:02d}.pdb"), structs)
    foldseek(structs, os.path.join(HDIR, "tmalign.m8"), tmalign=True)
    foldseek(structs, os.path.join(HDIR, "foldseek.m8"), tmalign=False)
    tm, _ = load_m8(os.path.join(HDIR, "tmalign.m8"))
    _, ev = load_m8(os.path.join(HDIR, "foldseek.m8"))
    hh = hhalign_all(comps)

    pairs = []
    for a, b in itertools.combinations(sorted(comps), 2):
        both = a in has_struct and b in has_struct
        pairs.append([a, b, group[a], group[b], int(group[a] == group[b]),
                      f"{tm.get((a, b), 0.0):.3f}" if both else "NA",
                      f"{ev.get((a, b), 10.0):.2e}" if both else "NA",
                      f"{hh[(a, b)]:.1f}"])
    write_tsv(os.path.join(HDIR, "pairs.tsv"),
              ["a", "b", "group_a", "group_b", "same_group", "tm", "foldseek_e",
               "hhalign_prob"], pairs)

    def val(a, b, d, default):
        return d.get((min(a, b), max(a, b)), default)

    members = collections.defaultdict(list)
    for c, g in group.items():
        members[g].append(c)

    # --- per component -------------------------------------------------------
    rows, bad = [], []
    for c in sorted(comps):
        g = group[c]
        mates = [o for o in members[g] if o != c]
        others = [o for o in comps if group[o] != g]
        s_mates = [o for o in mates if o in has_struct] if c in has_struct else []
        s_others = [o for o in others if o in has_struct] if c in has_struct else []
        in_tm = max((val(c, o, tm, 0.0) for o in s_mates), default=None)
        in_hh = max((val(c, o, hh, 0.0) for o in mates), default=None)
        out_tm = max((val(c, o, tm, 0.0) for o in s_others), default=None)
        out_e = min(((val(c, o, ev, 10.0), o) for o in s_others), default=(None, None))
        out_hh = max(((val(c, o, hh, 0.0), o) for o in others), default=(None, None))
        if not mates:
            verdict = "alone"
        elif (in_tm is not None and in_tm >= TM_MIN) or (in_hh or 0) >= HH_PROB:
            verdict = "confirmed"
        else:
            verdict = "unconfirmed"
            bad.append(c)
        leak = []
        if out_e[0] is not None and out_e[0] <= FS_E:
            leak.append(f"foldseek c{out_e[1]:02d}")
        if out_hh[0] is not None and out_hh[0] >= HH_PROB:
            leak.append(f"hhalign c{out_hh[1]:02d}")
        fmt = lambda v, p: "NA" if v is None else format(v, p)
        rows.append([c, g, comps[c]["fold_key"], comps[c]["n_tips"],
                     "yes" if c in has_struct else "no",
                     fmt(in_tm, ".3f"), fmt(in_hh, ".1f"), fmt(out_tm, ".3f"),
                     fmt(out_e[0], ".2e"), fmt(out_hh[0], ".1f"), verdict,
                     ";".join(leak) or "-"])
    write_tsv(os.path.join(HDIR, "groups.tsv"),
              ["component", "fold_group", "fold_key", "n_tips", "structure",
               "best_tm_in_group", "best_hh_in_group", "best_tm_other_groups",
               "best_foldseek_e_other_groups", "best_hh_other_groups", "verdict",
               "cross_group_signal"], rows)

    # --- per group -----------------------------------------------------------
    summary = []
    for g in sorted(members):
        cs = sorted(members[g])
        st = [c for c in cs if c in has_struct]
        tms = [val(a, b, tm, 0.0) for a, b in itertools.combinations(st, 2)]
        best = [max(val(c, o, tm, 0.0) for o in st if o != c) for c in st] if len(st) > 1 else []
        cross = [r for r in rows if r[1] == g and r[11] != "-"]
        summary.append([f"g{g:02d}", comps[cs[0]]["fold_key"], comps[cs[0]]["family"],
                        len(cs), sum(int(comps[c]["n_tips"]) for c in cs), len(st),
                        f"{min(best):.2f}" if best else "NA",
                        f"{sorted(tms)[len(tms) // 2]:.2f}" if tms else "NA",
                        sum(1 for r in rows if r[1] == g and r[10] == "unconfirmed"),
                        ";".join(f"c{r[0]:02d}:{r[11]}" for r in cross) or "-"])
    os.makedirs(OUT, exist_ok=True)
    write_tsv(os.path.join(OUT, "homology_check.tsv"),
              ["fold_group", "fold_key", "family", "components", "tips",
               "with_structure", "weakest_best_tm_in_group", "median_tm_in_group",
               "unconfirmed_components", "cross_group_signal"], summary)

    same = [float(p[5]) for p in pairs if p[4] == 1 and p[5] != "NA"]
    diff = [float(p[5]) for p in pairs if p[4] == 0 and p[5] != "NA"]
    print(f"TM-score, same group: median {sorted(same)[len(same) // 2]:.2f} "
          f"(n={len(same)}); different groups: median "
          f"{sorted(diff)[len(diff) // 2]:.2f}, max {max(diff):.2f} (n={len(diff)})")
    for r in summary:
        if r[3] > 1:
            print(f"  {r[0]} {r[1]:<14} {r[3]:>2} components: weakest link TM {r[6]}, "
                  f"{r[8]} unconfirmed")
    print("unconfirmed members: " + (", ".join(f"c{c:02d}" for c in bad) or "none"))
    leaks = [r for r in rows if r[11] != "-"]
    print("cross-group homology signal: " +
          (", ".join(f"c{r[0]:02d} ({r[11]})" for r in leaks) or "none"))


if __name__ == "__main__":
    main()
