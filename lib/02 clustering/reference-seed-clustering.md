# Reference-seeded clustering (generalized pipeline)

Reference-guided greedy clustering — **closed-then-open reference OTU picking**. A curated
**seed set** is fixed as authoritative cluster representatives; the remaining sequences are
first recruited to those seeds (closed reference), and whatever is left is clustered de novo
(open reference). This guarantees every curated seed is a representative **without changing
the clustering criteria applied to the rest of the set**.

Use it when you have hand-picked reference sequences that must survive as centroids (e.g.
named enzymes that anchor downstream stages), but you still want unbiased clustering of
everything else.

Two engines are documented and interchangeable: **USEARCH** (`cluster_fast` / `usearch_global`)
and **DIAMOND** (`greedy-vertex-cover`). Concrete runs:
`plasticome-union-v1.1/reference-seed-clustering/` (usearch) and
`plasticome-clustering-diamond/reference-seed-clustering-diamond/` (diamond).

---

## Inputs (parameters)

| Parameter | Meaning | Example |
|---|---|---|
| `SEEDS` | FASTA of curated reference sequences (the anchors) | 13 pazy seeds |
| `REMAINING` | FASTA of everything else (universe minus seeds, no overlap) | 597 no-seed |
| `$id` | identity threshold | 0.90 (usearch) / 90 (diamond) |
| `$cov` | member coverage rule | member ≥90% covered by its representative |
| engine | USEARCH or DIAMOND | — |

`SEEDS ∪ REMAINING` is the whole universe; the two files must be disjoint. Labels must be
unique and stop characters (`*`) stripped.

---

## Pipeline

### Phase 0 — build the substrate once
Concatenate `SEEDS + REMAINING`. Build a **single** all-vs-all similarity substrate at `$id`
and reuse it for every phase — never rebuild per phase (that invites parameter drift).
- USEARCH: no persistent graph; each step re-searches (see command table).
- DIAMOND: build one edge graph — `diamond blastp` all-vs-all, `--id $id`,
  `--query-or-subject-cover 90`, output `qseqid sseqid pident qcovhsp scovhsp bitscore`.

### Phase 1 — seed diagnostic (are the seeds non-redundant?)
Cluster the seed set **against itself** at `$id`. Report how many seeds remain centroids:
- **All seeds distinct** → all are reference centroids; proceed.
- **Some collapse** (one seed ≥`$id` to another) → the seed set is redundant at this
  threshold. **Decide per collapse**: either drop the redundant seed (curation) or keep both
  as designated centroids (accepting near-duplicate reps). Document the choice.

  *Example:* at 90% the provided 14-seed set collapsed BhrPETase (GBD22443.1) into LCC
  (AEV21261.1) at 92.2% id → BhrPETase was dropped, leaving 13 non-redundant seeds.

### Phase 2 — closed-reference recruitment
Map `REMAINING` against the seed centroids at `$id`. Each remaining sequence that matches a
seed (best hit, subject to the coverage rule) **joins that seed's cluster**. Matched
sequences leave the pool.
- USEARCH: `usearch_global REMAINING -db SEEDS -id $id -top_hit_only`.
- DIAMOND: recruit `v` to seed `s` where an edge `s→v` exists with `scovhsp ≥ $cov`
  (v is ≥90% covered by the seed); on contention assign to the seed with highest
  `(scov, bitscore, pident)`.

### Phase 3 — open-reference (de novo)
Cluster the **unmatched remainder** at `$id` into fresh clusters.
- USEARCH: `cluster_fast unmatched -id $id -sort length`.
- DIAMOND: `greedy-vertex-cover` on the residual node list, same `--member-cover $cov`, same
  edge graph, passing only the residual ids as `-d`.

### Phase 4 — merge & emit
Concatenate reference clusters (seeds + recruits) and de-novo clusters. Emit:
- `all-centroids.fasta` — the representatives (all seeds + all de-novo centroids).
- `cluster_membership.tsv` — `cluster_id, size, role (centroid|member), label,
  representative, pident, membercov`.
- verify `Σ cluster members == |SEEDS ∪ REMAINING|`.

Cluster ids: `1..N`, ordered earliest to latest, with the reference (seed) clusters
numbered **first** in seed order — so with the 13 curated seeds, ids 1–13 are always
the reference clusters and 14+ are de novo. There is no `origin` column: that block
ordering is the sole marker of a cluster's origin, so the two groups must be emitted
in this order.

---

## Engine command map

| Step | USEARCH v11 | DIAMOND v2.x |
|---|---|---|
| substrate | (implicit per search) | `diamond blastp … --query-or-subject-cover 90` (once) |
| 1 seed self-cluster | `cluster_fast SEEDS -id $id -sort length` | `greedy-vertex-cover` on seeds |
| 2 closed-ref recruit | `usearch_global -db SEEDS -id $id -top_hit_only` | edge `s→v`, `scovhsp ≥ $cov` |
| 3 de-novo | `cluster_fast unmatched -id $id -sort length` | `greedy-vertex-cover` on remainder, `--member-cover $cov` |

Greedy ordering differs and is worth stating in any report: USEARCH `cluster_fast` orders by
**length** (`-sort length`); DIAMOND `greedy-vertex-cover` orders by **graph degree**. Neither
is "length-disadvantaging" in the other's sense — the length-sort artifact is a USEARCH
concern; degree-based centroid choice is a diamond/GVC concern.

---

## Coverage semantics — and a caveat that propagates

Two different coverage rules commonly appear in the same pipeline and must be reconciled:
- **Edge/graph stage** often uses `--query-or-subject-cover 90` (**OR**): an edge exists if
  *either* sequence is ≥90% covered.
- **Member assignment** uses member-only coverage (`scovhsp ≥ 90`): the *recruited* sequence
  must be ≥90% covered by its representative.

A short sequence fully covered by a long one produces an **edge** but can fail **member
assignment** — so it is recruited by a long centroid but not by a short one. This shifts
boundary sequences between clusters purely on the coverage definition, independent of
identity. Pick one rule (e.g. `--mutual-cover`, or align both stages on subject coverage)
and apply it consistently, because reference clusters here define seed sets for downstream
stages and the inconsistency propagates.

Engines also differ in what they recruit: usearch `usearch_global` recruits on **identity
alone**; diamond enforces `--member-cover`. Expect diamond to recruit fewer sequences to
seed clusters and push the rest to de-novo (in the pazy run: 25 vs 36 recruited).

---

## Reporting checklist

1. **Are all curated seeds representative centroids?** (Phase 1 result; note any collapses
   and the curation decision.)
2. **General cluster stats:** total records; total clusters (reference + de-novo); count
   recruited to seeds vs de-novo; singletons (% of clusters); largest cluster; size
   distribution; where a key sequence (e.g. IsPETase) landed.
3. **Engine/parameter provenance:** engine + version, `$id`, coverage rule, greedy ordering,
   and (usearch on arm64) the Docker `linux/amd64` wrapper.

---

## Variant: forcing a *specific* representative inside the de-novo set

If the goal is narrower — guarantee one particular sequence becomes a representative that
greedy clustering would otherwise absorb — seed **only that anchor** (and optionally its
≥95% variant family) rather than the whole curated set, then run open-reference on the rest.
Report the **membership delta vs the unseeded baseline** (how many nodes actually change
cluster): a small delta means the intervention is cosmetic (it only relabels the
representative of an existing group); a large delta means greedy centroid choice was doing
real structural work that the seeding overrides. See
`plasticome-clustering-diamond/ispetase-anchor-experiment/` for the arms design (A unseeded /
B single anchor / C variant family / D full seed set) across identity thresholds.

---

## Reproducibility notes

- DIAMOND runs natively. USEARCH `usearch11/12` are Linux x86-64 ELF binaries — on arm64
  macOS run them inside a Docker `linux/amd64` container (Docker Desktop must be up).
- Fix the identity threshold, coverage rule, and greedy ordering explicitly; defaults differ
  between tools and versions.
- Deduplicate or note exact-duplicate sequences: they co-cluster and can inflate a
  centroid's degree (in the pazy set, BAB86909's three duplicate copies inflated its GVC
  degree).
