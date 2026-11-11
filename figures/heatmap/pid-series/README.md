# Pairwise identity heatmap at a series of clustering identities

The heatmap of `../03_sequence_heatmap.py` (every centroid against every
other, one cell per pair) drawn on centroids at 90, 80, 70, 60, 50, 40 and 30
percent identity.

```
PY=/opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python
cd .. && $PY 06_pid_series.py            # about 50 s, rebuilds this folder
$PY 06_pid_series.py --ids 95,85,75      # any other series
```

## What is here

- `overview.png`, all seven heatmaps on one sheet
- `summary.tsv`, one row per identity
- `idNN/figure_heatmap_sequences_idNN.{svg,pdf,png}`, the heatmap
- `idNN/clusters.tsv`, the cluster table, one row per cluster with its members
- `idNN/sequence_order_idNN.tsv` and `sequence_matrix_idNN.tsv`, the rows of
  the heatmap and the full matrix (`NA` for a pair that does not align)
- `idNN/clustering.log` and `idNN/work/`, the clustering run
- `union.tsv`, the 609 records as the clustering script reads them

## How the centroids are made

Each identity is one run of the project's own
`lib/02 clustering/cluster_reference_seeded.py` with its defaults: DIAMOND
engine, a member must be covered 90 percent by its representative, and the 13
curated seeds are kept as representatives whatever the identity.

| identity | centroids | alpha/beta hydrolase | share | pairs aligned |
|---:|---:|---:|---:|---:|
| 90 | 403 | 301 | 75% | 31.5% |
| 80 | 338 | 241 | 71% | 24.8% |
| 70 | 287 | 194 | 68% | 19.3% |
| 60 | 233 | 148 | 64% | 14.2% |
| 50 | 184 | 109 | 59% | 10.6% |
| 40 | 128 | 81 | 63% | 8.4% |
| 30 | 97 | 68 | 71% | 8.4% |

## Reading the series

- The coloured triangle of SSN component 0 shrinks at every step, because the
  centroids are what is left after its close relatives are merged. At 50
  percent and below the matrix is mostly grey.
- The hydrolase share of the centroids stays between 59 and 75 percent.
- All 16 folds and all 56 SSN components are present at every identity.

## Things that could be misread

- **The 90 percent figure here has 403 centroids, not 411.** The 411 are the
  rows marked `centroid = yes` in the union table, from an earlier clustering
  run. The script today gives 403, and 375 sequences are in both sets. The 411
  figure is still `../out_centroids/figure_heatmap_sequences.*`. This series
  uses one method for every identity so the panels can be compared.
- **Pairs above the clustering identity remain.** Two centroids at 50 percent
  can still be 60 percent identical, for two reasons. The curated seeds are
  never merged, so LCC, TfCut_2 and Cut190 stay separate rows at every
  identity. And a sequence joins a cluster only when its representative covers
  90 percent of it, so two proteins of different length stay apart.
- **97 centroids at 30 percent against 56 SSN components.** The SSN joins by
  single linkage with 50 percent coverage. The clustering here needs every
  member within reach of one representative at 90 percent coverage, which is
  stricter.
- **Rows are in the order of the tree of all sequences** (`../../tree/runs/all-sequences/results`),
  so a centroid keeps its place from one panel to the next. The 411 figure
  uses the centroid tree instead and is ordered a little differently.
- **One representative is missing at 30 percent** (96 drawn of 97). It is the
  28 residue arylacylamidase fragment, which the trees leave out.
- The pale yellow cells in the grey areas are chance alignments under 30
  percent, as in the 411 figure.
