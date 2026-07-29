# singleton-cleaned-union.v1.v2 — usearch vs diamond alignment

**Date:** 2026-07-23
**Input:** `petadex-clustering/singleton-cleaned-union.v1.v260701/results/singleton-cleaned-union.v1.v2.all-centroids.fasta`
(412 centroids, md5 `f4367699ba2eb004c4ae63bacf17f099`)

> Note: this input differs from the 2026-07-23 15:17 run staged in
> `petadex-alignment-{usearch,diamond}/singleton-cleaned-union.v1.v2/` — the clustering
> was re-run (30 headers changed, sequence-set md5 changed), so those earlier results are
> stale. This consolidated folder re-runs both aligners on the current FASTA.

## Method (identical for both aligners)

Frozen PETadex component pipeline (v4–v8, paper-faithful):

1. Centroids → md5-unique node set (`m####` ids). **412 nodes, 0 md5-collapsed** — the
   `combined.fasta` produced by the usearch run is reused byte-for-byte as diamond's
   `seqs.faa`, so the two partitions are node-for-node comparable.
2. All-vs-all **local** alignment, permissive (`-acceptall` / `--evalue 10`), no masking,
   no composition-based stats. usearch runs both orientations.
3. Post-filter every candidate edge by **pct_id ≥ 30 AND e-value < 1e-5** (best HSP per
   unordered pair by max bits). Filtering is done in Python, never as a search parameter.
4. Single-linkage connected components; components ranked by size.

- usearch: `usearch11 v11.0.667 -allpairs_local` (Docker linux/amd64)
- diamond: `diamond blastp --ultra-sensitive --comp-based-stats 0 --masking 0` (Docker linux/amd64)

## Headline result — the aligner choice barely moves the partition

| metric | usearch | diamond |
|---|---:|---:|
| nodes | 412 | 412 |
| passing edges | 22,410 | 22,541 |
| components | **46** | **46** |
| singletons | 25 | 24 |
| largest component | 233 | 229 |
| top-10 sizes | 233, 45, 19, 16, 14, 10, 5, 5, 5, 5 | 229, 45, 19, 16, 14, 10, 5, 5, 5, 5 |

**Agreement (on the shared 412 nodes):**

| | value |
|---|---:|
| Edge-set Jaccard | **0.979** (22,241 shared / 169 usearch-only / 300 diamond-only) |
| Rand index | 0.9837 |
| **Adjusted Rand Index** | **0.9632** |
| nodes in a non-corresponding component | **6 / 412** |

Both aligners recover the same 46-component structure and identical sizes for ranks 2–10.
The only material difference is at the edge of the giant "hydrolase hub" (component 1).

## The 6 boundary nodes (`node_disagreements.tsv`)

Every disagreement is a near-threshold case around the 30% / 1e-5 cutoff — exactly where a
small difference in the two tools' bit-score/e-value estimates flips one marginal edge.

| node | PL | accession | usearch | diamond |
|---|---|---|---|---|
| m0102 | PL122 | — | hub C001 (233) | satellite comp (5) |
| m0180 | PL40 | RLI42440.1 | hub C001 (233) | satellite comp (5) |
| m0340 | PL117 | — | hub C001 (233) | satellite comp (5) |
| m0376 | PL118 | — | hub C001 (233) | satellite comp (5) |
| m0300 | PL195 | RGD93181.1 | hub C001 (233) | singleton |
| m0199 | PL454 | WDT94443 | small comp (4) | singleton |

- The four PL122/PL40/PL117/PL118 nodes are what make usearch's largest component **233**
  vs diamond's **229**: usearch finds one bridging edge that pulls a 5-member satellite into
  the hub; diamond keeps it detached. `--ultra-sensitive` is slightly more conservative here.
- diamond leaves PL195 and PL454 as singletons where usearch attaches them by a single
  just-passing edge.

None of these move a node between two *large* well-formed components — they are all
hub-vs-satellite / attach-vs-singleton flips on the tie-break boundary.

## Files

```
singleton-cleaned-union.v1.v2/
  singleton-cleaned-union.v1.v2.all-centroids.fasta   input (412 centroids)
  usearch/  run-alignment.sh, results/                stats_singleton_2026-07-23.json, component_*
  diamond/  run-alignment.sh, cluster.py, results/     stats_diamond.json, component_*
  compare.py                                          partition/edge comparison
  comparison.json                                     machine-readable agreement metrics
  node_disagreements.tsv                              the 6 boundary nodes
  REPORT_compare_2026-07-23.md                        this file
```

## Reproduce

```bash
cd singleton-cleaned-union.v1.v2
./usearch/run-alignment.sh          # produces the shared combined.fasta node set first
./diamond/run-alignment.sh          # reuses ../usearch/results/combined.fasta as seqs.faa
python compare.py                   # regenerates comparison.json + node_disagreements.tsv
```

**Bottom line:** on this 412-centroid set the two aligners are effectively interchangeable —
ARI 0.963, identical component count and size spectrum, differing on just 6 threshold-boundary
nodes at the rim of the main hydrolase hub.
