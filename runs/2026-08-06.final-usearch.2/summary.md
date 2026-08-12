# 2026-08-06.final-usearch.2

Summary written 2026-08-11 by `lib/summary/summarize_run.py`.

## Inputs

- v1.1: `source-data/plasticome.v1.1/plasticome.v1.1.no-seq.tsv`
- v260701: `source-data/plasticome.v260701/cleaned_pazy-260701-singletons.tsv`
- seeds: `source-data/plasticome-curated-seeds/plasticome-curated-seeds.tsv`
- v260701 sequence rule: `from-accession`
- clustering identity: 0.90

## Deliverables

| step | file | rows |
| --- | --- | --- |
| 1 | `01-union.tsv` | 609 |
| 2 | `02-clusters.tsv` | 411 |
| 3 | `03-alignment.tsv` | 411 |
| 4 | `04-final-usearch.2.fasta` | 411 |

## 1. Union

- 609 rows — v260701 397 / both 76 / v1.1 136
- sequences: 609 present, 0 blank
- accession lookups: 565 resolved, 0 unresolved
  - resolved by: NCBI 473, UniProt 53, RCSB 4, nuccore CDS 3, MGnify 6, UniParc 26

## 2. Clusters (0.90 identity)

- 411 clusters over 609 sequences
- 272 singleton clusters, largest 21
- engine: usearch v11.0.667_i86linux64 — centroids ordered by length (cluster_fast -sort length), member coverage not enforced (identity only)
- representatives — v260701 324 / both 56 / v1.1 18 / seed 13
- representative length: min 186 / median 315 / max 734 aa

## 3. Components

- 46 components over 411 clusters
- 25 single-cluster components; largest: C001 (232), C005 (45), C003 (19)
- aligner: usearch v11.0.667 -allpairs_local -acceptall, both orientations
- edge criterion: pct_id>=30.0 AND evalue<1e-05 @ 144560-residue db
- 30588 aligned pairs, 22203 edges passing
- 166 clusters carry a v1.1 overlay; 0 v1 labels split across components

## 4. FASTA

- 411 centroid records in `04-final-usearch.2.fasta`
