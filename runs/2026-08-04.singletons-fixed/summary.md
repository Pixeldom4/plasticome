# 2026-08-04.singletons-fixed

Summary written 2026-08-04 by `lib/summary/summarize_run.py`.

## Inputs

- v1.1: `source-data/plasticome.v1.1/plasticome.v1.1.tsv`
- v260701: `source-data/plasticome.v260701/cleaned_pazy-260701-singletons.tsv`
- seeds: `source-data/plasticome-curated-seeds/plasticome-curated-seeds.tsv`
- v260701 sequence rule: `from-accession`
- clustering identity: 0.90

## Deliverables

| step | file | rows |
| --- | --- | --- |
| 1 | `01-union.tsv` | 610 |
| 2 | `02-clusters.tsv` | 403 |
| 3 | `03-alignment.tsv` | 403 |
| 4 | `04-singletons-fixed.fasta` | 403 |

## 1. Union

- 610 rows — v260701 398 / both 75 / v1.1 137
- sequences: 596 present, 14 blank
- accession lookups: 544 resolved, 14 unresolved
  - unresolved: `1ETH_A`, `5ZOA_A`, `7C2A_A`, `MBN1171866.1`, `MGYP000221121644`, `MGYP000271253724`, `MGYP000321434903`, `MGYP000532440779`, `MGYP001121581011`, `MGYP001477358452`, `MGYP001489421514`, `OP972509` (+2 more)

## 2. Clusters (0.90 identity)

- 403 clusters over 609 sequences
- 265 singleton clusters, largest 21
- representatives — v260701 314 / both 55 / v1.1 21 / seed 13
- representative length: min 186 / median 315 / max 734 aa

## 3. Components

- 46 components over 403 clusters
- 26 single-cluster components; largest: C001 (225), C005 (45), C003 (19)
- edge criterion: pct_id>=30.0 AND evalue<1e-05 @ 141600-residue db
- 29003 aligned pairs, 20795 edges passing
- 161 clusters carry a v1.1 overlay; 0 v1 labels split across components

## 4. FASTA

- 403 centroid records in `04-singletons-fixed.fasta`
