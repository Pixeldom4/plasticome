# Plasticome forest: component trees grafted on a structural backbone

The plasticome v2 union set cannot be aligned as one family, so it has no
single tree. This folder builds what the data supports instead: one tree per
SSN component, joined under a structure based backbone where the components
are homologous, and left as separate trees where they are not.

It implements the recommendation in the vault note "Plasticome Global Tree -
Grafting Recommendation" and answers the open questions in it. Input is the
same table as panel a-i,
`../plasticome.v2 - plasticome-v1-v260701-union-annotated (1).csv`, and the
components come straight from `../ssn/out/ssn_nodes_c30_d40_k3.tsv`.

## Answers to the note's open questions

**What is a component.** A connected component of the SSN at 30 percent
identity over 50 percent coverage, the family cluster definition of panel a-i.
There are 56.

**Has component level homology been checked.** Yes, in step 4, three ways. The
CATH annotation holds up. 36 components are alpha/beta hydrolases and are
homologous to each other. The other 20 fall into 14 unrelated folds plus one
sequence that could not be confirmed.

**Tree or forest.** A forest of 16 trees. 4 are grafts of several components
(alpha/beta hydrolase with 36, beta-lactamase like with 4, trypsin like with
2, cupredoxin with 2), 6 are a single component, and 6 are a single sequence.

**Can the hydrolase components be ordered relative to each other.** Mostly
not. Of 33 internal backbone branches, 4 reach UFBoot 95 and 3 of those are
recovered by an independent tree. The backbone in the figure is therefore
close to a star. This is the Holmes and Duchene criticism of the RdRp tree
showing up in this data, and the figure reports it instead of hiding it.

## Run

```
PY=/opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python

$PY 01_components.py        # components, duplicates, fold groups
$PY 02_component_trees.py   # Muscle5 ensemble + IQ-TREE per component (slow)
$PY 03_structures.py        # one AlphaFold DB or ESMFold model per component
$PY 04_homology.py          # TM-score, Foldseek, HHalign between components
$PY 05_backbone.py          # FoldMason + IQ-TREE backbone per fold group
$PY 06_graft.py             # graft, collapse unsupported backbone branches
$PY 07_draw.py              # out/figure_tree.{svg,pdf,png}, fan plus side trees
$PY 08_draw_pie.py          # out/figure_pie.{svg,pdf,png}, one pie sliced by fold
$PY 09_cath_table.py        # out/cath_domains.xlsx, CATH domains coloured by fold
```

Python runs in the `plasticome` env. The binaries live in a separate env,
`plasticome-tree`, so `plasticome` is not changed:

```
conda create -n plasticome-tree -c conda-forge -c bioconda \
  muscle=5.3 iqtree fasttree foldseek foldmason hhsuite mafft
```

`PLASTICOME_TREE_BIN` points the scripts at a different bin directory and
`PLASTICOME_TREE_THREADS` sets the thread count.

Step 2 is the only slow one. Component 0 has 284 sequences, and its 16
alignment ensemble takes about 2 minutes per alignment on 10 cores, then
about 80 minutes of IQ-TREE on 6 (928 search iterations). It reuses a finished `ensemble.efa` and resumes an
interrupted IQ-TREE search, so it can be stopped and rerun. On a laptop on
battery, run it under `caffeinate -i`, otherwise idle sleep stretches it to
several hours.

Step 3 makes network calls: UniProt (accession lookup by sequence checksum),
AlphaFold DB, and for 5 sequences the ESM Atlas folding API. Everything is
cached in `work/structures/cache`, so a rerun is offline. `--no-esmfold`
skips the folding API.

## Two sets: all sequences, or the 90 percent centroids

`PLASTICOME_TREE_SET` picks which sequences the trees are built from. Every
script reads it, and each set has its own directories, so the two runs sit
side by side.

| set | records | tips | directories |
|---|---:|---:|---|
| `all` (default) | 609 | 492 | `work/`, `out/` |
| `centroids` | 411 | 411 | `work_centroids/`, `out_centroids/` |

```
export PLASTICOME_TREE_SET=centroids
$PY 01_components.py && $PY 02_component_trees.py && ... && $PY 09_cath_table.py
```

`centroids` keeps the records with `centroid = yes` in the union table. All
411 are distinct sequences and none is a fragment. All 13 curated seeds are
centroids, so none is lost.

Components are the SSN components of the full set in both runs, so a component
number means the same family here as in panel a-i. Every component still
has at least one centroid. This matters for two of them: on a 30 percent graph
drawn from the centroids alone, component 0 would fall into 3 pieces and
component 9 into 2, because some of their links run through non centroid
members. They are kept whole.

44 of the 56 components contain only centroids. Their sequences are the same
in both sets, so step 2 copies the finished alignment and tree from `work/`
instead of redoing them. 12 components lose members and are rebuilt:
component 0 goes from 284 to 225 tips, the rest lose 1 to 5. The
representative structure of a component can differ between the sets, since it
is the most central tip that has a model and the centre moves when members
drop out. The backbone is therefore rebuilt too.

Everything below describes the `all` run unless it says otherwise. The
numbers of the `centroids` run are in "Centroid run" near the end.

## Step 1, components

609 records become 493 distinct sequences. The union set lists a protein twice
when both PAZy pulls carry it, 116 times in all. Identical sequences are
collapsed to one tip named for the smallest plasticome_id, and `records` in
the tip table lists every record a tip stands for. Where the records of one
tip disagree on substrate, the tip is `Combi.`.

The trees cover the whole union set, not only the 90 percent centroids. Of the
492 tips, 411 are centroids (`centroid = yes`) and 81 are not.

One sequence is left out of inference as a fragment: plasticome_id 56,
Arylacylamidase P80008.1, 28 residues in a component with a median of 479.
The rule is shorter than half the component median. That leaves 492 tips.

Fold groups are proposed from the CATH superfamily of the catalytic domain.
The table shows them as they stand after the step 4 check:

| fold group | CATH | components | tips |
|---|---|---:|---:|
| Alpha/beta hydrolase | 3.40.50.1820 | 36 | 377 |
| L-amino peptidase D-ALA esterase/amidase | 3.60.70.12 | 1 | 48 |
| Amidase signature domain | 3.90.1300.10 | 1 | 24 |
| Oxygenase MpaB family | none | 1 | 11 |
| Peptidase S8/S53 | 3.40.50.200 | 1 | 8 |
| Trypsin like serine proteases | 2.40.10.10 | 2 | 6 |
| DD-peptidase/beta-lactamase | 3.40.710.10 | 4 | 5 |
| Cytochrome c like | 1.10.760.10 | 1 | 3 |
| Cupredoxin | 2.60.40.420 | 2 | 2 |
| SGNH hydrolase | 3.40.50.1110 | 1 | 2 |
| six single sequences | various | 6 | 6 |

## Step 2, one tree per component

For the 14 components with 4 or more tips:

1. Muscle5 stratified ensemble, 16 alignments.
2. The alignment with the highest mean column confidence is the reference.
   Columns with more than 50 percent gaps or column confidence below 0.5 are
   dropped.
3. IQ-TREE, ModelFinder over LG, WAG, JTT and Q.pfam, 1000 ultrafast
   bootstraps.
4. Ensemble bootstrap. 100 alignments are made by drawing columns with
   replacement from all 16 alignments pooled, a FastTree LG+gamma tree is
   built on each, and each branch of the ML tree gets the share of the 100
   trees that contain it.

The two supports measure different things and both are kept in the tree as
`UFBoot/ensemble`. UFBoot resamples one fixed alignment and cannot see
alignment error. The ensemble bootstrap draws from alignments that disagree,
so it is the stricter of the two and is the one the figure uses. Component 0
shows the gap: median UFBoot 96, median ensemble support 52. Across the 14
components, 168 of 398 internal branches reach ensemble support 70.

Component 0 is also where alignment uncertainty bites hardest. Its reference
alignment has 1598 columns for 284 sequences, mean column confidence 0.57, and
148 columns pass the two filters. The tree inside the largest component rests
on those 148 columns.

| component | tips | columns kept | median UFBoot | median ensemble |
|---:|---:|---:|---:|---:|
| 0 | 284 | 148 | 96 | 52 |
| 1 | 48 | 263 | 66 | 24 |
| 2 | 24 | 361 | 99 | 99 |
| 3 | 18 | 195 | 77 | 77 |
| 5 | 11 | 403 | 100 | 100 |
| 4 | 10 | 335 | 55 | 49 |
| 6 | 8 | 366 | 100 | 100 |
| 7 | 7 | 413 | 80 | 70 |
| 9 | 7 | 161 | 75 | 84 |
| 8 | 5 | 404 | 99 | 99 |
| 10 | 5 | 607 | 97 | 98 |
| 11 | 5 | 536 | 74 | 55 |
| 12 | 4 | 331 | 98 | 95 |
| 14 | 4 | 243 | 82 | 59 |

Muscle5 has a `-resample` command for step 4. It is redone in numpy here
because it takes no seed and the output should be identical run to run.

The 9 smaller components get one alignment each. 3 tips: IQ-TREE LG+G4, no
support because there is one topology. 2 tips: a FastTree LG+gamma distance.

`work/component_trees.tsv` has the model, column counts and support summary
for every component.

## Step 3, representative structures

Components share under 30 percent identity with each other, so they are
compared on structure. Each component is represented by one model: the most
central tip (highest mean identity to the rest) that has an AlphaFold DB entry
with a matching sequence.

| source | components |
|---|---:|
| AlphaFold DB, sequence identical or contained | 48 |
| AlphaFold DB, 99.6 percent identical (porcine pancreatic lipase, 1ETH chain against P00591) | 1 |
| ESMFold through the ESM Atlas API | 5 |
| none | 2 |

The 5 ESMFold models are Est11_32, Ces39_5, APET-23, TA26 and MnP, mean pLDDT
83 to 95. The 2 without a model are IS11 (455 aa) and the PA4-degrading enzyme
(716 aa), both over the API's 400 residue limit. In 5 components the
representative is not the medoid because the medoid has no model; the largest
shift is component 0, rank 5 of 284.

## Step 4, which components can share a tree

Every component against every other, three measures:

- TM-score from a global superposition of the representatives (Foldseek in
  TM-align mode), normalised by the shorter chain
- Foldseek 3Di plus amino acid E-value
- HHalign probability between profiles built from each component's alignment

Result: median TM-score 0.65 inside a fold group and 0.27 between fold groups.
Every alpha/beta hydrolase component reaches TM 0.60 or better with another
one. No pair from different fold groups has a Foldseek E-value under 1e-3 or
an HHalign probability over 90, so no fold groups need merging.

Two things worth knowing:

**TM-score alone would have merged three superfamilies.** The SGNH hydrolase
and subtilisin components reach TM 0.5 to 0.63 against alpha/beta hydrolases.
All three are three layer alpha/beta/alpha sandwiches (CATH 3.40.50), so the
overall shape matches, but the best Foldseek E-value between them is 0.026 and
CATH puts them in different homologous superfamilies. They are kept apart.
This is why the rule uses TM-score to confirm membership inside a group and
Foldseek and HHalign to test for links between groups.

**HHalign has little power here.** The profiles are built only from the
sequences in this set, many from one or two sequences, and the median
probability inside a fold group is 2. It is useful for a component with no
structure and as a check for cross group signal, not as a general test.

**One component is not confirmed.** The PA4-degrading enzyme (component 46) is
annotated CATH 3.40.710.10 but has no structure and an HHalign probability of
0.3 against the other beta-lactamase like components. It is taken out of that
group and listed as a single sequence marked "homology unconfirmed". IS11, the
other one without a structure, is confirmed by HHalign at 98.

`out/homology_check.tsv` is the per group summary and
`work/homology/pairs.tsv` has all three scores for all 1540 pairs.

## Step 5, the backbone

One tip per component. FoldMason aligns the representative structures and
writes the alignment twice, as amino acids and as 3Di structural states.
Columns present in fewer than half the representatives are dropped. IQ-TREE
runs both alphabets as two partitions with a shared topology, the amino acid
model picked by ModelFinder and 3Di under the 3DI matrix of Puente-Lelievre et
al. (`models/3DI`, provenance in `models/README.md`), with 1000 ultrafast
bootstraps and the NNI correction. Branch lengths are then re-estimated on the
amino acid partition alone so the backbone is in the same units as the
component trees.

For the alpha/beta hydrolases, 324 of 866 columns survive and the median
pairwise identity on them is 10.6 percent. That is far below where an amino
acid alignment carries reliable signal, which is why the 3Di partition is
there and why every branch is cross checked against two trees that use the
data differently:

- `aa_only`: ML on the amino acid partition alone
- `tm_nj`: neighbour joining on 1 minus TM-score, no alignment columns at all

| | branches |
|---|---:|
| internal backbone branches, alpha/beta hydrolase | 33 |
| UFBoot 95 or more | 4 |
| also recovered by aa_only or tm_nj | 3 |
| shared between the combined tree and aa_only | 6 |

The three that survive are small: components 32, 33 and 45; components 9, 21,
26 and 28; components 34 and 48. The fourth, a clade of 7 components (18, 19,
23, 27, 36, 50, 55) with UFBoot 99, is recovered by neither cross check and is
collapsed by default. `06_graft.py --no-cross-check` keeps it.

Standard bootstrap was the first choice for the backbone and was dropped on
cost, over a minute per replicate. UFBoot is calibrated differently, so the
threshold is 95, not 70.

The beta-lactamase like backbone has 4 tips and one internal branch, MGS0010
with IS11, UFBoot 100 and recovered by the amino acid tree. IS11 has no
structure and joins the amino acid alignment by sequence (mafft --add
--keeplength), missing in the 3Di partition. The trypsin like and cupredoxin
groups have 2 components each, so their backbone is one branch.

`out/backbone_branches.tsv` lists every backbone branch with its UFBoot and
both cross checks.

## Step 6, grafting

1. Root the backbone at its midpoint.
2. Collapse every backbone branch below UFBoot 95 or not recovered by a cross
   check. 30 of 33 collapse for the hydrolases.
3. Root each component tree at its midpoint.
4. Replace each backbone tip with its component's rooted tree.

The stem that carries a component is the backbone's terminal branch for the
representative, less the distance from the component's own root to that
representative, so that stretch is not counted twice. `out/stems.tsv` has the
numbers.

Internal node names in the Newick files say what kind of branch sits above
the node:

| name | meaning |
|---|---|
| `B_<ufboot>_<aa>_<tm>` | backbone branch, with the two cross checks as 0 or 1 |
| `S_c<NN>` | stem of component NN, the root of its subtree |
| `C_<ufboot>_<ensemble>` | branch inside a component |
| `X` | no support defined: a root, or a tree of 2 or 3 tips |

## Step 7, the figure

The alpha/beta hydrolase tree is a fan, the other trees are rectangular beside
it, and single sequence fold groups are listed. Tips use the panel a-i
encoding, imported from `../ssn/build_ssn.py` so the two panels cannot drift.

- tip colour: substrate
- tip shape: diamond bacterial, hexagon eukaryotic, circle metagenomic
- ringed tip with a name: curated seed enzyme
- dashed black: backbone between components, with its UFBoot at the node
- dark grey: branch inside a component with ensemble support 70 or more, or a
  terminal branch
- light grey: branch inside a component with ensemble support under 70
- outer arc and shaded sector (fan only): one SSN component, numbered when it
  has 4 or more tips

**Backbone branches are not drawn to scale.** Their estimated lengths run from
3 to 11 substitutions per site. At 10 percent identity that is saturation, not
a measurement, and drawn to scale the stems would be several times longer than
any component tree. They are drawn at a fixed length and dashed. The estimates
stay in the Newick files and in `out/stems.tsv`. Every solid branch in the
figure is to scale and shares the one scale bar.

Three layout choices in the fan are for legibility and carry no meaning:

- The root sits on an inner circle instead of at a point. With 36 stems
  leaving one point, their dashes line up into rings that are not in the tree.
- Component 0 has 284 of the 377 tips and would take three quarters of the
  circle. Its tips are spaced at 0.6 of a row, which leaves the other 35
  components about 40 percent of the circle.
- A single sequence component gets two rows, since its marker sits near the
  centre where rows are narrow.

The canvas is 22.1 by 12.7 inches. The SVG has real `<text>` elements and one
named group per tree, same as panel a-i, so a tree can be moved as a unit in
Inkscape and its labels retyped. Steps 6 and 7 are byte reproducible: two runs
give the same `forest.nwk` and `figure_tree.svg`.

## Step 8, the same forest as one pie

`08_draw_pie.py` draws the same trees in a second layout,
`out/figure_pie.{svg,pdf,png}`. Every fold group is one slice of a single
circle, the slices are cut apart by a gap, and the background of a slice is
coloured by its fold. Nothing joins two slices. Each tree starts on the inner
circle on its own, because trees in different fold groups share no ancestor in
this analysis.

- slice: one fold group, 16 in all, largest first, clockwise from the top
- slice colour: the fold, one pale tint each. The fold is also named outside
  its slice with a matching swatch, so colour is never the only cue. The 6
  fold groups that are a single sequence share one grey, since 16 tints cannot
  be told apart
- fold labels carry the name only. The CATH accession and the component and
  sequence counts are in the fold table, a figure of its own
  (`out/figure_pie_table.{svg,pdf,png}`) so it can be placed independently,
  and also written to `out/fold_table.tsv` (with the record count and the homology
  verdict as extra columns). Two slices are DD-peptidase/beta-lactamase; the
  one that failed the homology check is labelled "(unconfirmed)" so the two
  can be told apart
- two tints inside a slice: neighbouring components alternate, so the 36
  components of the hydrolase slice stay visible
- tips, branches, seeds, the dashed backbone: as in step 7
- tip names are off by default and the 13 curated seeds are called out by
  name in every slice. `--tip-names` puts every tip name back along the rim in
  the small slices (not the hydrolase one, where 377 names do not fit)

All slices share one scale. A slice gets angle in proportion to its tips, with
three adjustments so small folds stay readable: tips in the small slices take
1.35 rows (with or without names, so the slices do not change width),
component 0 takes 0.6 of a row per tip, and no slice is narrower than 4 rows.
The hydrolase slice ends up at 55 percent of the circle.

The tints are background, behind substrate coloured tips, so they are kept
pale and are not a validated categorical palette. The hydrolase slice is
mostly PET blue tips and has a warm tint for contrast. To change a fold's
colour edit `FOLD_TINT` at the top of the script.

SVG groups, for editing in Inkscape:

```
slice-g00 ... slice-g15     one per slice: background, tree, tips, tip names
fold labels                 parent group
   fold-label-g00 ...       one per fold: leader line, swatch, name
legend                      parent group
   legend-branches
   legend-substrate
   legend-origin
   legend-scale
```

and in `figure_pie_table.svg`, one parent group `fold table` holding
`table-header`, one `table-row-gNN` per fold, and `table-total`.

The fold labels and the legend can each be moved whole through the parent, or
one piece at a time.

Step 7 and `out/figure_tree.*` are untouched by this script, and
`backup_v1_fan/` holds a copy of the step 7 script and its three figure files
as they were on 2026-10-01.

## Step 9, the CATH domains as a table

`09_cath_table.py` writes `out/cath_domains.xlsx`, with every row filled in
its fold's tint from the pie, and the same data as `out/cath_folds.tsv` and
`out/cath_domains.tsv`.

- **Folds**: one row per slice. Catalytic CATH accession and class, the fill
  colour as hex, components, sequences, records, sequences per substrate and
  per origin, sequences with a PDB structure, curated seeds, the accessory
  domains found beside the catalytic one, the homology verdict from step 4,
  the weakest TM-score link, and backbone branches kept.
- **Domains**: one row per CATH domain per fold it occurs in, 34 rows for 31
  distinct domains (16 catalytic, one of them the MpaB family with no
  accession, and 15 accessory). An accessory domain shared by
  two folds has two rows; the galactose-binding domain is the one case, in an
  alpha/beta hydrolase component and in the unconfirmed beta-lactamase
  sequence.
- **Notes**: what every column means.

The `domain` column of the union table joins names with ";", and some names
contain "; " themselves. A ";" with no space after it is taken as the
separator, which lines the names up one to one with the CATH ids for all 609
records. The script stops if that ever fails or if one accession turns up
with two names.

Two records have no CATH accession. The oxygenase MpaB family has a name but
no accession and is shown as `none`. The Arylacylamidase fragment has no
domain at all and is the "No domain annotated" row under the amidase fold.

The Total row and the share column are spreadsheet formulas. LibreOffice is
not installed on this machine, so they were not recalculated here and carry no
cached values until the file is opened in Excel, Numbers or Google Sheets.

## What the figure does and does not claim

- A fork drawn in a solid dark line inside a component is supported under
  alignment uncertainty.
- A fork drawn in the dashed backbone has UFBoot 95 or more and is recovered
  by an independent tree. Everything else between components is a polytomy.
- Roots are midpoint roots, a drawing convention. No outgroup supports the
  root of any component or of the backbone.
- Trees in different fold groups are not related to each other by anything in
  this analysis. Their placement on the page carries no meaning.
- The order of components around the fan carries no meaning beyond the three
  supported backbone clades.

## To decide

- Whether Artem's "graft via muscle5" meant profile to profile merging of the
  component alignments (hidden `profalign` and `-guidetreein`). That would
  give one alignment and one ML tree, but between hydrolase components the
  merged columns would sit near 10 percent identity. The structural backbone
  here is the more defensible version of the same idea, and it already shows
  how little is resolvable at that depth.
- Whether to show the star backbone as is, or to drop the backbone and present
  the hydrolase components as a forest too. The data supports either reading.
- Whether to keep the clade of 7 components that has UFBoot 99 but fails both
  cross checks.
- A structure for the PA4-degrading enzyme and IS11 (both over 400 residues)
  would settle the one unconfirmed assignment.

## Outputs

In `out/`:

- `figure_tree.{svg,pdf,png}`, the fan with side trees (step 7)
- `figure_pie.{svg,pdf,png}`, the pie sliced by fold (step 8)
- `figure_pie_table.{svg,pdf,png}`, the fold table for the pie
- `fold_table.tsv`, the table drawn on the pie: fold, CATH accession,
  components, sequences, records, homology verdict
- `cath_domains.xlsx`, `cath_folds.tsv`, `cath_domains.tsv`, the CATH domain
  tables (step 9)
- `trees/gNN.nwk`, one rooted tree per fold group, and `forest.nwk` with all
  of them, one per line
- `tree_tips.tsv`, one row per tip: fold group, component, order in the tree,
  records it stands for, substrate, shape, homology verdict
- `homology_check.tsv`, `backbone_summary.tsv`, `backbone_branches.tsv`,
  `stems.tsv`

In `work/`, regenerable: per component alignments, ensembles and IQ-TREE
output, structures and their download cache, the pairwise homology scores, and
the backbone alignments.
