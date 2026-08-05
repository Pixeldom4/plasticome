# 2026-08-05.test-usearch

Summary written 2026-08-05 by `lib/summary/summarize_run.py`.

## Inputs

- v1.1: `source-data/plasticome.v1.1/plasticome.v1.1.no-seq.tsv`
- v260701: `source-data/plasticome.v260701/cleaned_pazy-260701-singletons.tsv`
- seeds: `source-data/plasticome-curated-seeds/plasticome-curated-seeds.tsv`
- v260701 sequence rule: `from-accession`
- clustering identity: 0.90

## Deliverables

| step | file | rows |
| --- | --- | --- |
| 1 | `01-union.tsv` | 610 |
| 2 | `02-clusters.tsv` | 413 |
| 3 | `03-alignment.tsv` | 413 |
| 4 | `04-test-usearch.fasta` | 413 |

## 1. Union

- 610 rows — v260701 398 / both 75 / v1.1 137
- sequences: 610 present, 0 blank
- accession lookups: 566 resolved, 0 unresolved
  - resolved by: NCBI 474, UniProt 53, RCSB 4, nuccore CDS 3, MGnify 6, UniParc 26

## 2. Clusters (0.90 identity)

- 413 clusters over 610 sequences
- 275 singleton clusters, largest 21
- engine: usearch v11.0.667_i86linux64 — centroids ordered by length (cluster_fast -sort length), member coverage not enforced (identity only)
- representatives — v260701 325 / both 55 / v1.1 20 / seed 13
- representative length: min 186 / median 315 / max 734 aa

## 3. Components

- 47 components over 413 clusters
- 26 single-cluster components; largest: C001 (233), C005 (45), C003 (19)
- aligner: usearch v11.0.667 -allpairs_local -acceptall, both orientations
- edge criterion: pct_id>=30.0 AND evalue<1e-05 @ 145227-residue db
- 30859 aligned pairs, 22415 edges passing
- 167 clusters carry a v1.1 overlay; 0 v1 labels split across components

## 4. FASTA

- 413 centroid records in `04-test-usearch.fasta`
