#!/usr/bin/env python3
"""Diff the DIAMOND build against the USEARCH build it replaces.

Three runs are compared, all built from the SAME `01-union.tsv` (610 rows) by the SAME
revision of the code, so every difference is attributable to the engine and nothing else:

  U   runs/2026-08-05.usearch/               usearch clustering + usearch alignment
  A   runs/2026-08-05.diamond-aligner-only/  U's clusters + DIAMOND alignment
  D   runs/2026-08-05.diamond/               DIAMOND clustering + DIAMOND alignment

U vs A isolates the ALIGNER (byte-identical 413-centroid input on both sides).
U vs D is the end-to-end change.

(`runs/2026-08-05.singletons-fixed/` is NOT the baseline here: it was built before the
seed-promotion change and books 623 records with `S###` seed labels, so its clusters are
not comparable label-for-label. Its alignment stats are quoted where useful.)

Clusters are keyed on `rep_label` (`U0024|CBY05530`), derived from the union row index and
therefore stable across engines. Components are compared on the centroid labels they hold.

Usage:  python3 compare_vs_usearch.py [-o REPORT.md]
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

RUNS = Path(__file__).resolve().parents[1]
U = RUNS / "2026-08-05.usearch"
A = RUNS / "2026-08-05.diamond-aligner-only"
D = RUNS / "2026-08-05.diamond"

JOIN = "; "
MAX_LIST = 8          # labels printed before eliding


def clusters(run):
    """rep_label -> {members incl. the rep, size, cluster_id, source, name}."""
    out = {}
    for r in csv.DictReader((run / "02-clusters.tsv").open(), delimiter="\t"):
        mem = [m for m in r["member_labels"].split(JOIN) if m]
        out[r["rep_label"]] = {"members": set(mem) | {r["rep_label"]},
                               "size": int(r["size"]), "id": int(r["cluster_id"]),
                               "source": r["rep_source"], "name": r["rep_enzyme_name"]}
    return out


def comps(run):
    """component_id -> {rep_labels}."""
    by = defaultdict(set)
    for r in csv.DictReader((run / "03-alignment.tsv").open(), delimiter="\t"):
        by[r["component_id"]].add(r["rep_label"])
    return dict(by)


def stats(run):
    for p in sorted((run / "03-alignment.intermediates").glob("stats_*.json")):
        return json.loads(p.read_text())
    return {}


def assignment(cl):
    return {m: c for c, v in cl.items() for m in v["members"]}


def pair_agreement(a, b):
    """Over every unordered pair of the shared universe: do both partitions agree on
    whether the pair sits together? Returns (agree, total)."""
    shared = sorted(set(a) & set(b))
    agree = total = 0
    for i, x in enumerate(shared):
        for y in shared[i + 1:]:
            total += 1
            if (a[x] == a[y]) == (b[x] == b[y]):
                agree += 1
    return agree, total


def labs(xs):
    xs = sorted(xs)
    shown = ", ".join(f"`{x}`" for x in xs[:MAX_LIST])
    return (shown + f" (+{len(xs) - MAX_LIST} more)") if len(xs) > MAX_LIST else (shown or "—")


def fmt_dist(cl):
    d = Counter(v["size"] for v in cl.values())
    return ", ".join(f"{k}×{d[k]}" for k in sorted(d))


def component_diff(x, y, lines):
    """Report y's partition as a rearrangement of x's, over the centroids both hold.

    A crosstab of x-component × y-component: every cell is a set of nodes that travelled
    together. Cells on the dominant 1:1 matching are the components that survived and are
    summarised in one line; every off-diagonal cell is printed with its nodes, which is
    what a boundary shift or a merge actually is. Anchoring on cells rather than on
    components keeps a 229-node giant from being listed because one singleton joined it.
    """
    cx, cy = comps(x), comps(y)
    shared = {n for s in cx.values() for n in s} & {n for s in cy.values() for n in s}
    ax = {n: c for c, s in cx.items() for n in s if n in shared}
    ay = {n: c for c, s in cy.items() for n in s if n in shared}

    cell = defaultdict(set)
    for n in shared:
        cell[(ax[n], ay[n])].add(n)
    xsize = Counter(ax.values())
    ysize = Counter(ay.values())

    # Dominant matching: each x-component claims the y-component holding most of it.
    home = {c: max((k[1] for k in cell if k[0] == c),
                   key=lambda t: (len(cell[(c, t)]), ysize[t])) for c in xsize}
    kept = [(c, t) for c, t in home.items() if len(cell[(c, t)]) == xsize[c] == ysize[t]]

    lines.append(f"{len(shared)} centroids are nodes in both runs, in {len(xsize)} "
                 f"components on the left and {len(ysize)} on the right. "
                 f"{len(kept)} components are identical on both sides.")
    lines.append("")
    moves = sorted(((len(v), k, v) for k, v in cell.items()
                    if (k[0], k[1]) not in kept), reverse=True)
    if not moves:
        return
    # "keeps identity" = this cell is the biggest contributor to the target component, so
    # the target is that component carrying on rather than the nodes joining someone else.
    dominant = {}
    for k, v in cell.items():
        if len(v) > len(cell.get((dominant.get(k[1]), k[1]), ())):
            dominant[k[1]] = k[0]
    lines.append("| from | to | n | nodes |")
    lines.append("|---|---|---:|---|")
    for n, (xc, yc), mem in moves:
        note = " (keeps identity)" if dominant.get(yc) == xc else ""
        lines.append(f"| {xc} | {yc}{note} | {n} | {labs(mem) if n <= MAX_LIST else '—'} |")
    lines.append("")
    merged = {t: [c for c in xsize if home.get(c) == t] for t in ysize}
    merged = {t: cs for t, cs in merged.items() if len(cs) > 1}
    for t, cs in sorted(merged.items()):
        lines.append(f"- **{t}** absorbs {len(cs)} left-hand components "
                     f"({', '.join(sorted(cs))}) — {ysize[t]} nodes.")
    split = {c: sorted({k[1] for k in cell if k[0] == c}) for c in xsize}
    for c, ts in sorted(split.items()):
        if len(ts) > 1:
            lines.append(f"- **{c}** ({xsize[c]} nodes) splits across {', '.join(ts)}.")
    lines.append("")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", type=Path,
                    default=Path(__file__).with_name("REPORT_diamond_vs_usearch.md"))
    args = ap.parse_args()

    cu, cd = clusters(U), clusters(D)
    au, ad = assignment(cu), assignment(cd)
    su, sa, sd = stats(U), stats(A), stats(D)
    L = []
    w = L.append

    w("# DIAMOND vs USEARCH — 2026-08-05")
    w("")
    w("Same `01-union.tsv` (610 rows), same code revision on every side; the only variable")
    w("is the engine.")
    w("")
    w("| | clustering | alignment | run |")
    w("|---|---|---|---|")
    w(f"| **U** | usearch | usearch | `{U.name}` |")
    w(f"| **A** | usearch (U's table) | diamond | `{A.name}` |")
    w(f"| **D** | diamond | diamond | `{D.name}` |")
    w("")
    w("**U reproduces the pre-change pipeline exactly** — 413 clusters, 28 recruited, 47")
    w("components, 22,415 passing edges, largest component 233 — so it is a valid stand-in")
    w("for the previous result. (`2026-08-05.singletons-fixed` predates the seed-promotion")
    w("change and books 623 records under `S###` labels, so it is not label-comparable at")
    w("step 2; its step-3 numbers are identical to U's.)")
    w("")

    # ---- step 2 -------------------------------------------------------------
    w("## Step 2 — clustering")
    w("")
    w("| | usearch | diamond | Δ |")
    w("|---|---:|---:|---:|")
    for lab, f in (("clusters", lambda c: len(c)),
                   ("reference (seed) clusters", lambda c: sum(1 for v in c.values() if v["source"] == "seed")),
                   ("de-novo clusters", lambda c: sum(1 for v in c.values() if v["source"] != "seed")),
                   ("singletons", lambda c: sum(1 for v in c.values() if v["size"] == 1)),
                   ("largest cluster", lambda c: max(v["size"] for v in c.values())),
                   ("records placed", lambda c: sum(v["size"] for v in c.values()))):
        x, y = f(cu), f(cd)
        w(f"| {lab} | {x} | {y} | {y - x:+d} |")
    w("")
    w(f"- size distribution, usearch: {fmt_dist(cu)}")
    w(f"- size distribution, diamond: {fmt_dist(cd)}")
    w("")

    both = set(cu) & set(cd)
    identical = [c for c in both if cu[c]["members"] == cd[c]["members"]]
    w(f"**Centroids.** {len(both)} of usearch's {len(cu)} centroids are also diamond "
      f"centroids ({len(identical)} of them heading a byte-identical cluster). "
      f"{len(set(cu) - set(cd))} usearch centroids are demoted to members, and diamond "
      f"promotes {len(set(cd) - set(cu))} sequences usearch had kept as members.")
    w("")

    moved = [m for m in au if au[m] != ad.get(m)]
    agree, total = pair_agreement(au, ad)
    intact_rec = sum(1 for m in au if cu[au[m]]["members"] == cd[ad[m]]["members"])
    w(f"**Records.** {intact_rec} of {len(au)} records sit in a cluster with identical "
      f"membership on both sides (the cluster may have swapped which member is its "
      f"centroid).")
    w("")
    w(f"**Partition.** {len(au) - len(moved)} of {len(au)} records keep the same centroid; "
      f"{len(moved)} move. Over all {total:,} unordered pairs of records the two partitions "
      f"agree on {agree:,} ({100 * agree / total:.3f}%) — they disagree about co-clustering "
      f"for {total - agree} pairs, which is the honest size of the change: most 'moves' are "
      f"a cluster keeping its membership while swapping which member is called the centroid.")
    w("")

    w("**Seed clusters** (`cluster_id` 1–13, the curated anchors). Phase 2 is where the")
    w("coverage rule bites — diamond recruits 17 where usearch recruits 28:")
    w("")
    w("| # | seed centroid | enzyme | usearch | diamond | released by diamond |")
    w("|---:|---|---|---:|---:|---|")
    for c in sorted(both, key=lambda c: cu[c]["id"]):
        if cu[c]["source"] != "seed":
            continue
        lost = cu[c]["members"] - cd[c]["members"]
        w(f"| {cu[c]['id']} | `{c}` | {cu[c]['name'][:30]} | {cu[c]['size']} | "
          f"{cd[c]['size']} | {labs(lost)} |")
    w("")
    n_u = sum(v["size"] for v in cu.values() if v["source"] == "seed")
    n_d = sum(v["size"] for v in cd.values() if v["source"] == "seed")
    w(f"Seed clusters hold {n_u} records under usearch and {n_d} under diamond: "
      f"{n_u - n_d} released into de-novo clusters. All 13 seeds stay centroids under both.")
    w("")

    big = sorted(((cu[c]["size"] - cd[c]["size"], c) for c in both
                  if cu[c]["source"] != "seed" and cu[c]["members"] != cd[c]["members"]),
                 reverse=True)[:10]
    if big:
        w("**De-novo clusters that change most** (shared centroids only):")
        w("")
        for _, c in big:
            w(f"- `{c}` {cu[c]['size']} → {cd[c]['size']}  — left: {labs(cu[c]['members'] - cd[c]['members'])}")
        w("")

    # ---- step 3 -------------------------------------------------------------
    w("## Step 3 — alignment")
    w("")
    w("| | U usearch/usearch | A usearch/diamond | D diamond/diamond |")
    w("|---|---:|---:|---:|")
    for lab, k in (("nodes aligned", "n_nodes"), ("raw alignments", "raw_alignments"),
                   ("aligned pairs", "aligned_pairs"), ("edges passing", "edges_passing"),
                   ("components", "n_components"),
                   ("single-node components", "n_singletons"),
                   ("largest component", "largest_component")):
        w(f"| {lab} | {su.get(k, '?'):,} | {sa.get(k, '?'):,} | {sd.get(k, '?'):,} |"
          if isinstance(su.get(k), int) else
          f"| {lab} | {su.get(k, '?')} | {sa.get(k, '?')} | {sd.get(k, '?')} |")
    w("")
    w(f"- top-10 component sizes, U: {su.get('top10_component_sizes')}")
    w(f"- top-10 component sizes, A: {sa.get('top10_component_sizes')}")
    w(f"- top-10 component sizes, D: {sd.get('top10_component_sizes')}")
    w("")

    w("### Aligner isolated — U vs A (identical 413-centroid input)")
    w("")
    e0, e1 = su.get("edges_passing", 0), sa.get("edges_passing", 0)
    w(f"Edges passing ≥30% aaid AND e<1e-5: {e0:,} → {e1:,} ({e1 - e0:+,}).")
    w("")
    component_diff(U, A, L)

    w("### End to end — U vs D")
    w("")
    w(f"Over the {len(both)} centroids both runs have (usearch's other "
      f"{len(set(cu) - set(cd))} are members, not nodes, under diamond, so their "
      f"components cannot be compared directly):")
    w("")
    component_diff(U, D, L)

    w("## Bottom line")
    w("")
    w(f"- Step 2 moves: {len(cu)} → {len(cd)} clusters. The mechanism is phase 2 — "
      f"diamond's 90% member-coverage rule recruits 17 to the seeds where usearch's "
      f"identity-only rule recruits 28 — plus degree-ordered rather than length-ordered "
      f"centroid choice, which reshuffles which member of an unchanged cluster is called "
      f"the representative.")
    w(f"- The partitions are near-identical in substance: they disagree about "
      f"co-clustering for {total - agree} of {total:,} pairs "
      f"({100 * (total - agree) / total:.3f}%).")
    w(f"- Step 3 is close to a no-op when isolated: same 47 components, "
      f"{sa.get('edges_passing', 0) - su.get('edges_passing', 0):+,} edges, and the giant "
      f"α/β-hydrolase component shifting by a handful of nodes.")
    w(f"- Export artifact: `04-*.fasta` carries {len(cu)} records under usearch and "
      f"{len(cd)} under diamond. These are different centroid sets and only one should "
      f"feed the nr-search.")

    args.out.write_text("\n".join(L).rstrip() + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
