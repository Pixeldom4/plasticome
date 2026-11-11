# Identity heatmap, plasticome v2

A triangular matrix of amino acid identity across the plasticome, after panel
c of the Serratus Coronaviridae figure (`coronaviridae-heatmap.webp`, the
reference image). In the Serratus panel the colour tiers are order, family and
genus. Here they are fold, family, subfamily and centroid.

There are two figures. The main one is step 3: all 411 centroids against each
other, one cell per pair, `out_centroids/figure_heatmap_sequences.*`. Step 2
is the earlier version, 21 groups with one median per cell.

Input is the same table as the SSN and the tree,
`../plasticome.v2 - plasticome-v1-v260701-union-annotated (1).csv`. Families
are the SSN components of `../ssn` and folds are the fold groups of `../tree`.

## Run

```
PY=/opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python

$PY 01_identity.py     # all against all alignment, about 35 s on 10 cores
$PY 03_sequence_heatmap.py   # out_centroids/figure_heatmap_sequences.*

$PY 02_heatmap.py                                  # grouped version, out/
PLASTICOME_TREE_SET=centroids $PY 02_heatmap.py    # grouped, out_centroids/
```

Step 1 runs once and serves both sets. Step 3 draws the centroids unless
`PLASTICOME_TREE_SET=all` is given, which draws the 492 distinct sequences
into `out/`. Step 2 follows `../tree` and defaults to all. Step 2 takes
`--min-seqs`, `--split`, `--sub-id` and `--no-numbers`. Step 3 takes
`--no-seed-labels`.

## Step 1, identity for every pair

The SSN was built from a DIAMOND search, which reports only the pairs it
finds. A median over all pairs of two groups needs every pair, so step 1
aligns all of them: 492 distinct sequences, 120,786 pairs, Smith-Waterman with
BLOSUM62 and gap open 11, extend 1 (the BLAST defaults), through Biopython.

The sequences are the tips the trees kept in `../tree/runs/all-sequences/intermediate/tips.tsv`, so
identical records count once and the one fragment is out.

Identity is identical columns over alignment columns, the pident of BLAST and
DIAMOND. Coverage is alignment columns over the length of the shorter
sequence, the definition `../ssn/build_ssn.py` uses.

Check against the SSN search: for the 40,486 pairs DIAMOND also reported, the
median difference in identity is 0.0 and the 95th percentile is 0.7 points.
DIAMOND had found 94 percent of the pairs that pass the coverage rule here.

## Step 3, all 411 centroids against each other

One cell per pair, 84,255 pairs, nothing pooled and nothing left out.

**Order.** Rows run fold by fold, largest first, and inside a fold in the
order of its tree in `../tree/runs/centroids/results` (`order_in_tree`). Each SSN
component is one clade of that tree, so a component is a run of neighbouring
rows, and close relatives sit together inside it. The strip on the left marks
the 56 components in two alternating greys. A hairline in the matrix marks
where a component ends and a wider one where a fold ends.

**Cells.** The identity of the pair from step 1. A pair that does not align
over half of the shorter sequence has no identity and is grey, darker when the
two share a fold. The 13 curated seeds are called out from the diagonal.

**The bars above the legend** are the distribution itself: the share of the
aligned pairs in each bin.

| | pairs | share |
|---|---:|---:|
| aligned | 27,272 | 32.4% |
| not aligned, same fold | 22,179 | 26.3% |
| not aligned, different fold | 34,804 | 41.3% |

| identity of the aligned pairs | share |
|---|---:|
| under 30 | 18.7% |
| 30 to 40 | 18.7% |
| 40 to 50 | 37.3% |
| 50 to 60 | 17.2% |
| 60 to 70 | 5.8% |
| 70 to 80 | 1.2% |
| 80 to 90 | 1.0% |
| 90 and over | 7 pairs |

The median over the aligned pairs is 43.6 percent. Component 0, with 225
centroids, holds 88 percent of them.

Three things to read off the figure:

- The coloured triangle in the middle is component 0. It has visible
  subfamily blocks at 50 to 90 percent on a background of 30 to 50, and its
  first rows (the top left of its block) align with the rest at under 30.
- The other 35 hydrolase components are the grey band above it. They share
  the fold and do not align with component 0 or with each other.
- Every other fold is a small triangle on the diagonal with light grey to its
  left.

Two things that could be misread:

- **The pale yellow specks in the grey areas are not homology.** 1,784
  aligned pairs join two different components, 624 of them two different
  folds, all between 17 and 30 percent. A local alignment of two proteins of
  similar length sometimes reaches half the shorter one by chance. Step 4 of
  `../tree` found no link between folds by structure or by profile.
- **7 centroid pairs are at 90 percent or over**, the highest 91.8. The
  centroids were cut at 90 percent by a different aligner, so a few pairs land
  just over the line here.

`sequence_order.tsv` has the row of every sequence and `sequence_matrix.tsv`
the full matrix, with `NA` for a pair that does not align.

## Step 6, the pairwise heatmap at a series of clustering identities

`06_pid_series.py` clusters the 609 records at 90, 80, 70, 60, 50, 40 and 30
percent with `lib/02 clustering/cluster_reference_seeded.py` and draws the
step 3 heatmap on the representatives of each run. Everything is in
`pid-series/`, one folder per identity plus `overview.png` and `summary.tsv`.
`pid-series/README.md` has the numbers and what to watch for. The main one:
the script gives 403 centroids at 90 percent, not the 411 marked in the union
table, and 375 sequences are in both.

Step 3 takes `--centroids <cluster table>`, `--out-dir`, `--tag` and `--title`
for this. Without them it draws the 411 as before, byte for byte.

## Step 4, diversity per fold

`04_fold_diversity.py` draws four figures for one statement: most distinct
sequences are alpha/beta hydrolases and the other 15 folds are a long tail.
A pairwise matrix cannot say this, because 41 percent of its pairs cross a
fold and have no identity to plot. These count clusters per fold instead.

- `figure_fold_cutoffs`: heatmap, fold by identity cut-off, clusters per cell,
  shaded by the fold's share of the clusters at that cut-off
- `figure_fold_bars`: ranked bars of centroids per fold, cut into SSN families
- `figure_fold_treemap`: area is centroids, fold then SSN family, in the fold
  tints of the pie
- `figure_fold_identity`: identity of the pairs inside each fold, as shares

Clusters at 30 percent are the SSN components. From 40 to 80 they are single
linkage clusters inside each component, on the pairs of step 1 that reach the
cut-off and cover half of the shorter sequence. At 90 they are the centroids.

| cut-off | clusters | alpha/beta hydrolase | share |
|---:|---:|---:|---:|
| 30 | 56 | 36 | 64% |
| 40 | 101 | 60 | 59% |
| 50 | 157 | 84 | 54% |
| 60 | 200 | 116 | 58% |
| 70 | 259 | 166 | 64% |
| 80 | 332 | 235 | 71% |
| 90 | 411 | 309 | 75% |

The hydrolases are the majority at every cut-off, by the least at 50 percent.
12 of the 16 folds have 5 centroids or fewer. `fold_cutoffs.tsv` has the
table for every fold.

Two limits on the wording. Diversity here is the number of distinct sequence
clusters; counted as unrelated folds, 15 of 16 are not hydrolases. And rare
means rarely reported, since PET dominates the source table.

## Step 5, the same statement as heatmaps only

`05_fold_heatmaps.py` redraws each step 4 figure as a heatmap. The 16 folds
are the rows in the same order every time, largest first. Run step 4 first.

- `figure_heat_cutoffs`: fold by identity cut-off, clusters per cell. The
  step 4 heatmap under this name.
- `figure_heat_families`: fold by SSN family rank, centroids per cell. The
  hydrolase row has 36 filled cells, no other row has more than 4.
- `figure_heat_tiles`: one tile per centroid, blocked by fold, coloured by its
  identity to the nearest other centroid of its fold. An alignment across
  folds is not counted, since those are chance. 8 centroids align with nothing
  in their fold, none of them a hydrolase. The median nearest neighbour is 74
  percent for hydrolases and 48 for the rest.
- `figure_heat_identity`: fold by identity bin, share of the fold's pairs.
- `figure_heat_substrate`: fold by substrate, centroids per cell. The
  hydrolases cover 9 of the 11 substrate groups and no other fold covers more
  than 3. They hold 238 of 242 PET centroids and none of the 54 PA ones.

Counts in the family and substrate heatmaps are shaded on a log scale so a
cell with 1 is still visible. The cut-off and identity heatmaps are shaded by
share, on a linear scale.

## Step 2, the grouped version: rows and cells

**A pair counts as aligned when the alignment covers half of the shorter
sequence**, the SSN rule. This matters because a local alignment always
returns something. Between unrelated proteins it returns a short stretch at 20
to 30 percent, which is not a measurement of anything.

**Cells are medians over all pairs**, one sequence from each group, with the
pairs that do not align ranked below every aligned pair. When half or more of
the pairs do not align, the cell has no median and is drawn grey. The diagonal
is the median inside a group.

**Rows are SSN components with 3 or more distinct sequences.** Component 0 is
cut into subfamilies first, by average linkage at 50 percent identity, and
each subfamily with 3 or more sequences gets a row.

Why component 0 is split and the rest are not: a matrix of the 56 components
as they are has a coloured diagonal and nothing else. Of the 30,427 pairs
between two components of the same fold, 4 percent pass the coverage rule, so
no pair of components has a median. That is what a component at 30 percent
means. The identity between groups that can be measured is inside component 0,
which holds 284 of the 492 sequences. `--split` takes a list if another
component should be cut too.

Rows are named for the curated seed enzymes in the group, or for the most
central named enzyme when there is no seed. The number in brackets is the
count of distinct sequences.

| | all sequences | centroids |
|---|---:|---:|
| sequences | 492 | 411 |
| rows | 21 | 21 |
| sequences in a row | 409 | 330 |
| off diagonal cells with a median | 21 of 210 | 21 of 210 |
| range of those medians | 25 to 48 | 39 to 50 |
| range of the diagonal | 29 to 85 | 27 to 66 |

The sequences not in a row are the 34 of component 0 that fall in subfamilies
of 1 or 2, and the 41 components with 1 or 2 sequences. They are listed when
step 2 runs.

## What the grouped figure shows

- Every cell with a median sits inside component 0. Its subfamilies are 25 to
  48 percent identical to each other and 54 to 85 percent inside themselves.
- No two components have a median between them, in the same fold or not. The
  15 alpha/beta hydrolase rows are one fold by structure (`../tree`, step 4)
  but only component 0 can be compared by sequence. The darker grey block is
  that statement.
- The other folds are separate islands, light grey against everything.
- Inside the other families the median runs from 29 to 54 percent. They are diverse:
  the nylon hydrolase and amidase families sit near 30.

Two cells to know about:

- **PhaZ5 (component 11) has no median on its own diagonal.** It has 5
  sequences and exactly half of its 10 pairs pass the coverage rule. It is
  one component only because the SSN joins by single linkage.
- **Row names move between the two sets.** The centroid run names a group for
  its most central centroid, so the largest subfamily is ESM121 in one figure
  and TEP041 in the other. `heatmap_groups.tsv` lists the members of each row.

## Colour

Bins are 10 points wide. Lightness falls from bin to bin across the whole
scale (OKLab L 0.94 down to 0.30, at least 0.06 per step), so the order reads
in greyscale and under colour blindness. The hue changes at the three
thresholds the project already uses, which is what the tier labels under the
legend mark:

| tier | identity | where the threshold comes from |
|---|---|---|
| ~Fold | under 30 | below the SSN component threshold |
| ~Family | 30 to 50 | SSN component, 30 percent |
| ~Subfamily | 50 to 90 | the `--sub-id` cut, 50 percent |
| ~Centroid | 90 and over | the 90 percent centroids |

The Serratus panel uses a separate hue per bin. That palette is not ordered in
lightness and several of its neighbours merge under deuteranopia, so it was
not copied. Each cell with a median also carries its number, so colour is not
the only channel. `--no-numbers` removes them.

## Outputs

In `out/` and `out_centroids/`:

- `figure_heatmap_sequences.{svg,pdf,png}`, step 3, in `out_centroids/`. The
  matrix is one embedded image in the SVG and PDF, at 6 pixels per cell; the
  labels, strip, brackets and legend are vector, with real text.
- `sequence_order.tsv` and `sequence_matrix.tsv`, step 3.
- `figure_heatmap.{svg,pdf,png}`, step 2. The SVG has real text and one element per
  cell with the id `cell-<row>-<col>`.
- `heatmap_matrix.tsv`, one row per cell: the two groups, the median, the
  number of pairs and the share of them that aligned.
- `heatmap_groups.tsv`, one row per row of the figure: name, fold, component,
  count, median inside the group and the member tips.

In `work/`: `pairs.tsv`, identity and coverage for all 120,786 pairs.

Two runs of step 2 or of step 3 give the same SVG, checked by md5.

## To decide

- Whether the chance alignments under 30 percent between components should be
  drawn at all. Greying every pair that crosses a component boundary would
  remove the specks, at the price of a rule that is no longer about the pair
  alone.
- The order inside component 0. Tree order is used so the figure matches the
  pie. Ordering by average linkage on identity would make the blocks squarer.
- Whether to keep step 2. Its subfamily cut at 50 percent is the one free
  choice in it: at 40 component 0 barely splits (248, 9, 4 and 4), at 60 it
  gives 17 rows of 3 or more.
