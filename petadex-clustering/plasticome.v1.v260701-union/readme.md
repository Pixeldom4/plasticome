# Reference-seeded clustering — plasticome union v1 v260701 @ 90% identity

Closed-then-open reference clustering (reference-guided greedy OTU picking): the **curated
seed set is fixed as authoritative reference centroids**, the remaining sequences are mapped
to those centroids (closed reference), and the unmatched remainder is clustered de novo (open
reference). This guarantees every curated seed survives as a representative without altering
the clustering criteria applied to the rest of the set. See
`../reference-seed-clustering.md` for the generalized pipeline.

**Tool:** USEARCH v11.0.667 (Docker `linux/amd64`; the binary is a Linux x86-64 ELF and will
not run natively on arm64 macOS). All steps at **90% amino-acid identity**, `-sort length`.

## Inputs

| File | Records | Role |
|---|---|---|
| `inputs/putative-curated-seed-set.fasta` | 13 | curated seeds (reference centroids) |
| `inputs/plasticome.v1.v260701-union.no-seeds.fasta` | 598 | union minus the seeds (`REMAINING`) |
| `inputs/plasticome.v1.v260701-union.fasta` | 611 | full universe (`SEEDS ∪ REMAINING`) |

Pre-run validation: the two source files are **disjoint** (0 overlapping labels), their union
equals the 611-record universe exactly, all labels are unique, and no `*` stop characters are
present. The 13 seeds are already non-redundant (BhrPETase/PL17 removed upstream), so no seed
curation was needed here.

## Pipeline

1. **Diagnostic** — `cluster_fast` the seed set alone (`-id 0.90 -sort length`) to test seed
   non-redundancy → **13/13 centroids, 0 collapses**: all seeds mutually non-redundant at 90%.
   (`work/01-seed.uc`, `work/01-seed.log`, `work/seed-centroids.fasta`)
2. **Closed-reference** — `usearch_global` the 598 no-seed vs the 13 seed centroids
   (`-id 0.90 -maxaccepts 0 -maxrejects 0 -top_hit_only`). Hits join the seed's cluster.
   → **37 matched**, **561 unmatched**. (`work/02-closedref.uc`, `work/matched.fasta`, `work/notmatched.fasta`)
3. **Open-reference (de novo)** — `cluster_fast` the 561 unmatched (`-id 0.90 -sort length`)
   → **399 new clusters**. (`work/03-denovo.uc`, `work/denovo-centroids.fasta`)
4. Merge into `results/cluster_membership.tsv`; concatenate reference + de novo centroids into
   `results/plasticome.v1.v260701-union.all-centroids.fasta` (412 total).

Reproduce end-to-end: `./scripts/run_clustering.sh` then `python3 scripts/build_membership.py`.

## Results

- **611 records → 412 clusters** (13 reference + 399 de novo).
- All **13 curated seeds are representative centroids** (guaranteed by construction; each also
  verified mutually non-redundant at 90% in the Phase 1 diagnostic).
- **37 no-seed sequences recruited** into **10 of the 13** seed clusters; **3 seed clusters
  stayed singletons** (seed only): R005 PL239, R006 PL245, R011 PL59.
- Singletons: **273** (66.3% of clusters). Largest cluster: **21** (R002, PL24 / CBY05530).
- Conservation verified: Σ cluster members = 611 = |universe|; every union label appears in
  `cluster_membership.tsv` exactly once (0 missing, 0 extra, 0 duplicate).
- IsPETase (A0A0K8P6T7 / GAP38373.1) is **not present** in this universe, so it anchors no cluster here.

Cluster size distribution (size : count):

| size | 1 | 2 | 3 | 4 | 5 | 6 | 21 |
|---|---|---|---|---|---|---|---|
| clusters | 273 | 114 | 14 | 5 | 3 | 2 | 1 |

Reference cluster sizes (cluster centroid : size):

| cluster | seed centroid | size |
|---|---|---|
| R002 | PL24 / CBY05530 | 21 |
| R009 | PL53 / WP_012476895.1 | 6 |
| R007 | PL51 / WP_012476897.1 | 4 |
| R001 | PL16 / AEV21261.1 (LCC) | 3 |
| R012 | PL48 / WP_011061486.1 | 3 |
| R000 | PL1 / WP_054022242.1 | 2 |
| R003 | PL18 / WP_015787089.1 | 2 |
| R004 | PL43 / A0A075B5G4 | 2 |
| R008 | PL52 / CAA24927.1 | 2 |
| R010 | PL85 / WP_011614907.1 | 2 |
| R005 | PL239 / SDZ16714.1 | 1 |
| R006 | PL245 / WP_116180173.1 | 1 |
| R011 | PL59 / WP_037320891.1 | 1 |

## Folder layout

```
plasticome.v1.v260701-union/
├── readme.md                        this file
├── inputs/                          source FASTAs (unmodified)
│   ├── putative-curated-seed-set.fasta               13 curated seeds
│   ├── plasticome.v1.v260701-union.no-seeds.fasta    598 remaining
│   └── plasticome.v1.v260701-union.fasta             611 full universe
├── scripts/
│   ├── run_clustering.sh            3-phase USEARCH driver (Docker wrapper) + all-centroids emit
│   └── build_membership.py          merge/parse → results/cluster_membership.tsv
├── work/                            intermediate USEARCH records, logs, interim FASTAs
│   ├── seed-centroids.fasta         13 seed centroids used as the closed-reference DB
│   ├── matched.fasta / notmatched.fasta   no-seed split by closed-reference (37 / 561)
│   ├── denovo-centroids.fasta       399 de novo centroids
│   ├── 01-seed.uc / 02-closedref.uc / 03-denovo.uc   raw USEARCH cluster records
│   └── *.log                        USEARCH run logs
└── results/                         deliverables
    ├── plasticome.v1.v260701-union.all-centroids.fasta   412 representatives (13 ref + 399 de novo)
    └── cluster_membership.tsv       cluster_id, origin, size, role, label, centroid, pct_id_to_centroid
```

Cluster ids: `R###` = reference (seed) cluster, `D###` = de novo cluster.

## Notes

- USEARCH `cluster_fast` orders by decreasing length (`-sort length`); longest sequences
  become centroids first. `usearch_global` recruits on **identity alone** (no member-coverage
  rule) — a DIAMOND `--member-cover 90` run would recruit fewer to the seeds and push more to
  de novo (see `../reference-seed-clustering.md`, "Coverage semantics").
- USEARCH `bin/usearch11` is a Linux x86-64 ELF at
  `petadex-alignment/petadex-alignment-usearch/bin/` — Docker-only on this arm64 workstation.
