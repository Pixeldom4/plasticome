# 2026-08-05.test-diamond

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
| 2 | `02-clusters.tsv` | 404 |
| 3 | `03-alignment.tsv` | 404 |
| 4 | `04-test-diamond.fasta` | 404 |

## 1. Union

- 610 rows — v260701 398 / both 75 / v1.1 137
- sequences: 610 present, 0 blank
- accession lookups: 566 resolved, 0 unresolved
  - resolved by: NCBI 474, UniProt 53, RCSB 4, nuccore CDS 3, MGnify 6, UniParc 26

## 2. Clusters (0.90 identity)

- 404 clusters over 610 sequences
- 268 singleton clusters, largest 14
- engine: diamond version 2.2.4 — centroids ordered by graph degree (greedy-vertex-cover), member coverage 0.9
- representatives — v260701 229 / both 49 / v1.1 113 / seed 13
- representative length: min 186 / median 314 / max 734 aa

## 3. Components

- 48 components over 404 clusters
- 25 single-cluster components; largest: C001 (218), C005 (45), C003 (19)
- aligner: diamond blastp all-vs-all, --ultra-sensitive --max-target-seqs 0 --evalue 10 --masking 0 --comp-based-stats 0
- edge criterion: pct_id>=30.0 AND evalue<1e-05 @ 140831-residue db
- 27259 aligned pairs, 20736 edges passing
- 173 clusters carry a v1.1 overlay; 1 v1 labels split across components

## 4. FASTA

- 404 centroid records in `04-test-diamond.fasta`
