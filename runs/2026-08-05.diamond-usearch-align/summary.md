# 2026-08-05.diamond-usearch-align

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
| 4 | `04-diamond-usearch-align.fasta` | 404 |

## 1. Union

- 610 rows — v260701 398 / both 75 / v1.1 137
- sequences: 610 present, 0 blank

## 2. Clusters (0.90 identity)

- 404 clusters over 610 sequences
- 268 singleton clusters, largest 14
- representatives — v260701 229 / both 49 / v1.1 113 / seed 13
- representative length: min 186 / median 314 / max 734 aa

## 3. Components

- 47 components over 404 clusters
- 26 single-cluster components; largest: C001 (224), C005 (45), C003 (19)
- aligner: usearch v11.0.667 -allpairs_local -acceptall, both orientations
- edge criterion: pct_id>=30.0 AND evalue<1e-05 @ 140831-residue db
- 28755 aligned pairs, 20616 edges passing
- 173 clusters carry a v1.1 overlay; 0 v1 labels split across components

## 4. FASTA

- 404 centroid records in `04-diamond-usearch-align.fasta`
