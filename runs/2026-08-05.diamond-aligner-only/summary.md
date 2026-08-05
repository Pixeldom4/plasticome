# 2026-08-05.diamond-aligner-only

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
| 4 | `04-diamond-aligner-only.fasta` | 413 |

## 1. Union

- 610 rows — v260701 398 / both 75 / v1.1 137
- sequences: 610 present, 0 blank

## 2. Clusters (0.90 identity)

- 413 clusters over 610 sequences
- 275 singleton clusters, largest 21
- representatives — v260701 325 / both 55 / v1.1 20 / seed 13
- representative length: min 186 / median 315 / max 734 aa

## 3. Components

- 47 components over 413 clusters
- 25 single-cluster components; largest: C001 (229), C005 (45), C003 (19)
- aligner: diamond blastp all-vs-all, --ultra-sensitive --max-target-seqs 0 --evalue 10 --masking 0 --comp-based-stats 0
- edge criterion: pct_id>=30.0 AND evalue<1e-05 @ 145227-residue db
- 29343 aligned pairs, 22536 edges passing
- 167 clusters carry a v1.1 overlay; 1 v1 labels split across components

## 4. FASTA

- 413 centroid records in `04-diamond-aligner-only.fasta`
