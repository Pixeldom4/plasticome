# Panel a-i and a-ii, plasticome v2

Regenerates the sequence similarity network (panel a-i) and the family count
bar chart (panel a-ii) from the v2 union set, 609 sequences, replacing the
212 PAZy sequences behind the published version of the figure.

Input: `../plasticome.v2 - plasticome-v1-v260701-union-annotated (1).csv`

## Run

```
conda activate plasticome
PY=/opt/homebrew/Caskroom/miniconda/base/envs/plasticome/bin/python

# one time, all-vs-all alignment (already in work/allvall.tsv)
cd work
diamond makedb --in all.faa -d all
diamond blastp -q all.faa -d all -o allvall.tsv --very-sensitive --masking 0 \
  --comp-based-stats 0 -e 1e-3 --max-target-seqs 0 --outfmt 6 \
  qseqid sseqid pident length qlen slen bitscore evalue

cd ..
$PY build_ssn.py       # panel a-i
$PY build_aii.py       # panel a-ii
```

`build_ssn.py` takes `--component-id`, `--draw-id`, `--coverage`, `--topk`,
`--seed` and `--tag`, so alternative thresholds land in separate output files.

## Two thresholds, two jobs

The identity cutoff that defines a cluster and the cutoff that decides whether
to ink an edge are not the same number, and conflating them is what makes an
SSN unreadable at this size.

**Components, 30 percent identity over 50 percent coverage.** This is the
family cluster definition from panel a-iii. It is also where the data says the
boundary is. Scoring each connected component against the dominant CATH domain
family of its members:

| identity | components | singletons | largest | purity vs domain family |
|---------:|-----------:|-----------:|--------:|------------------------:|
| 20 | 30 | 10 | 416 | 94.3% |
| 25 | 33 | 11 | 409 | 94.3% |
| **30** | **56** | **21** | **329** | **99.8%** |
| 40 | 101 | 41 | 291 | 99.8% |
| 50 | 154 | 83 | 285 | 99.8% |

Purity jumps to 99.8 percent at 30 and holds there. Below 30 the alpha/beta
hydrolase component starts absorbing other families. So 30 is the loosest
threshold at which a component still means one family.

**Drawn edges, 40 percent identity, thinned to a backbone.** The published
panel's width legend runs 40 to 100, so 40 is the floor for inking an edge. A
flat 40 percent cutoff on 609 sequences is 37,097 edges, which is a solid grey
mass. The drawn set is the maximum spanning forest plus each node's 3 strongest
partners, giving 1,091 edges, or E/N 1.79. For calibration, the 212 PAZy
sequences in the published panel read at roughly E/N 1 to 1.5. Every drawn edge
is a real alignment at 40.4 percent identity or better; nothing is invented.

## Layout

Positions come from the component graph at 30 percent, not from the drawn graph
at 40 percent. A sequence whose only similarity to the rest of its family sits
in the 30 to 40 band has no drawn edge at all, and laid out on the drawn graph
it becomes a loose point that force layout throws clear of its own family blob.

The layout graph per component is the drawn backbone plus a maximum spanning
tree of that component's 30 percent graph. The backbone supplies the
substructure, the spanning tree guarantees the component stays in one piece.
Using the full 30 percent subgraph instead produces a featureless disc, because
inside a family that graph is nearly complete.

ForceAtlas2, linlog, strong gravity 4.0, with attraction weighted by
`((pident - 30) / 70) ** 3`, so a 90 percent pair pulls about 30 times harder
than a 40 percent pair. That condenses near identical groups into visible
subclusters while the component as a whole stays near the compact disc radius
`sqrt(n / pi)`. Every component is then rescaled so its median nearest
neighbour spacing is 1 unit, which lets one marker size and one blob radius
work across components of very different size.

## The halo means one component

The grey outline is the component boundary, so two nodes whose halos touch are
in the same component and a break in the halo means a different component. That
only reads correctly if each component's halo is a single connected region.

A plain union of discs of radius `BLOB_PAD` around each node does not
guarantee it: the region splits wherever two members sit further apart than
`2 * BLOB_PAD`. On this data that broke 5 of the 35 multi-node components, worst
of all component 2 (Amidase, 35 sequences) which rendered as 5 separate blobs
and read as 5 components.

`draw_blob()` therefore thresholds two distance fields at once: `BLOB_PAD`
around every node, and a narrower `bridge_pad` (0.42 of it) along a Euclidean
minimum spanning tree of the component's node positions. The spanning tree is
connected by construction, so the union is one region, and because its edges are
the shortest available the necks stay short rather than inflating the blob. A
convex hull was tried first and encloses large empty regions on sparse
components.

Components of *different* families are kept apart by the packing gap
(`2 * BLOB_PAD + 0.9`), so halos of different components can never touch. The
check in the Reproducibility section below verifies both properties.

## Sections

Components are grouped by catalytic family into 8 sections, each headed with
the family name and its component and sequence counts, largest family first:

| section | components | seqs |
|---|---:|---:|
| A/B hydrolase | 36 | 466 |
| L-amino peptidase D-ALA esterase | 1 | 53 |
| Amidase | 1 | 35 |
| Peptidase | 6 | 19 |
| Oxygenase-like | 1 | 15 |
| Other | 7 | 9 |
| Cytochrome C-like | 2 | 6 |
| Trypsin-like | 2 | 6 |

Packing runs at two levels, components within a section and sections on the
page, both through `shelf_pack()`. That is a bottom-left skyline packer rather
than a true shelf packer: a shelf gives every row the height of its tallest box,
so the 329 node component leaves a band of dead space beside everything shorter.
The skyline drops each box at the lowest point it fits and tucks the small
sections alongside the tall one.

`--inner-aspect` and `--outer-aspect` control the two levels. The defaults, 1.00
and 1.10, were chosen by sweeping both and taking the best fill at a roughly
square canvas (70 x 75 units, 44% filled). The run prints those numbers, so a
change can be checked. The canvas is then sized to the content, since the axes
are equal-aspect and a fixed figure size letterboxes the drawing.

## Editing the SVG in Inkscape

`figure_ai_*.svg` is grouped two levels deep, so either a whole family or a
single component can be dragged without pulling its halo, edges, nodes and seed
callouts apart:

```
A/B hydrolase (36 components, 466 seq)      <- section-a-b-hydrolase
   00 A/B hydrolase (cluster, n=329)        <- component-00
   03 A/B hydrolase (cluster, n=29)
   ...
Amidase (1 component, 35 seq)               <- section-amidase
   02 Amidase (cluster, n=35)
```

Section groups also hold the section heading text. Components are listed
largest first, singletons last. Click once to select the section, enter it to
select a component, enter again to reach individual marks.

All labels are real `<text>` elements, so headings, counts and seed names can be
selected and retyped in Inkscape. Matplotlib's default, `svg.fonttype="path"`,
traces text as vector outlines instead, which cannot be edited; both scripts set
`svg.fonttype="none"`. Panel a-i has 50 text elements and a-ii has 29. The
`<use>` references that remain in a-i are the node markers, which are meant to
stay as symbols.

The trade is that the font is referenced by name, not embedded. Editing needs
DejaVu Sans, which matplotlib ships, otherwise Inkscape substitutes and the
metrics shift a little. The PNG and PDF are unaffected and still carry real
glyphs.

Matplotlib writes each artist as its own top level element, so the grouping is a
post-processing pass, `regroup_svg()`, which moves the artists tagged with a
component's gid into one `<g>` with an `inkscape:label`.

Regrouping reorders elements, which would change paint order, so two things keep
it safe and the regrouped file pixel identical to the ungrouped one:

- components are packed with a gap of `2 * BLOB_PAD + 0.7`, so no two grey
  halos can overlap
- label placement scores candidate positions against `foreign_blob_hits()` and
  refuses to sit a label on another component's halo

Verified by rendering both with `rsvg-convert` and diffing: no differing pixels.
If you change `BLOB_PAD` or the packing gap, recheck that.

## Reproducibility

Repeated runs produce byte identical outputs. Getting there needed three fixes,
all of them cases where a set or a tie was resolved in Python's per process
string hash order:

- connected components are sorted by size then smallest member id, not by set
  iteration order
- the drawn backbone is sorted before the graph is built, since edge insertion
  order perturbs the force layout
- `max_spanning()` adds a canonical `1e-9` rank to each weight before running
  Kruskal. Many pairs share an identical pident, and the tied spanning trees
  have equal total weight but different edges, which moved every component

The SVG needed two more. Matplotlib salts element ids (clip paths, marker defs)
from a random uuid per figure unless `svg.hashsalt` is set, and it embeds a
`dc:date`, so two identical renders still produced different files. The salt is
fixed and the date is dropped via `metadata={"Date": None}`.

Verify with `md5 out/ssn_nodes_c30_d40_k3.tsv out/figure_ai_c30_d40_k3.svg`
across two runs. Both are stable.

## Checks worth re-running after a layout change

- **halo connectivity**: rasterise each component's distance field and count
  filled islands; every multi-node component must give exactly 1. Currently
  0 of 35 fail.
- **grouping is lossless**: render the SVG with and without `regroup_svg()`
  via `rsvg-convert` and diff the two PNGs. Must be pixel identical, which
  proves the reordering never changed what covers what. It holds because
  halos of different components cannot touch and labels avoid foreign halos.
- **reproducibility**: two runs, same md5.

If you change `BLOB_PAD`, `GAP` or either aspect, re-run all three.

## Encoding

- node colour: substrate
- node shape: diamond bacterial, hexagon eukaryotic, circle metagenomic
- edge width: aa identity percent, 40 to 100
- grey halo and its label: connected component at 30 percent, named for the
  dominant CATH domain family, labelled once per family on its largest component
- ringed node with a leader line: one of the 13 curated seed enzymes from
  `source-data/plasticome-curated-seeds/plasticome-curated-seeds.tsv`, matched
  into the union set by accession

## Palette

Sampled from the published panel, then re-stepped on the pairs that failed
`scripts/validate_palette.js` from the dataviz skill:

| substrate | published | here | why |
|---|---|---|---|
| PBAT | `#694da6` | `#9a6fd4` | was ΔE 2.0 from PET under deuteranopia, and PET is 56% of nodes |
| PHA | `#c9e7b7` | `#8fbf5f` | below the chroma floor, read as grey |
| PE | `#e82d88` | `#e8559e` | below the normal vision floor against PET |
| PUR | `#fdd871` | `#f0c95e` | lightness band |

PET, PA, PLA, PBS, NR and Combi. are unchanged from the published panel. PC is
new in v2 and had no published colour.

Two checks do not pass and are accepted. PUR sits above the lightness band:
darkening it collides with PHA under protanopia, since yellow and green are
close there, so it stays light and every node carries a thin dark stroke as
relief. PA against PE is ΔE 6.2 under deuteranopia, inside the 6 to 8 band that
is permitted with secondary encoding; PE is 2 nodes and both are direct
labelled. Eleven substrate categories is more than colour alone separates
cleanly, and folding the rare ones into Other was rejected because each named
plastic carries meaning here.

## Known consequence of the v2 annotation

Four of the curated seeds render black as `Combi.` rather than PET blue, because
v2 gives them a second substrate:

| seed | v2 substrate |
|---|---|
| ISPETase | PEF; PET |
| LCC | PC; PEF; PET; PUR |
| TfCut_2 | PEF; PET; PUR |
| HiC | PBAT; PEF; PET |

Any sequence with more than one substrate maps to `Combi.`, which is 43 of 609.
If PEF should be its own colour, or if these should show their primary
substrate, change `substrate_of()` in `build_ssn.py` and both panels follow.

## Panel a-ii binning

`build_aii.py` bins the network's family groups into the five published bars.
L-amino peptidase D-ALA esterase/amidase is counted under Amidase, which is
where most of its 53 members' PA activity belongs. Trypsin-like is counted
under Peptidase and Cytochrome C-like under Oxygenase. Change `FAMILY_BINS` to
move them.

| family | count |
|---|---:|
| A/B hydrolase | 466 |
| Amidase | 87 |
| Peptidase | 25 |
| Oxygenase | 21 |
| Other | 10 |
| total | 609 |

## Outputs

In `out/`:

- `figure_ai_c30_d40_k3.{svg,pdf,png}` and `figure_aii.{svg,pdf,png}`
- `ssn_c30_d40_k3.graphml`, node and edge attributes carried, plus the
  computed x and y so Cytoscape opens on this layout
- `ssn_nodes_c30_d40_k3.tsv`, one row per sequence with substrate group,
  family group, component, component size, drawn degree, colour, shape, x, y
- `ssn_edges_c30_d40_k3.tsv`, every pair at 40 percent or above with a `drawn`
  flag marking the backbone
- `figure_aii_counts.tsv`, the family by substrate table behind the bars
