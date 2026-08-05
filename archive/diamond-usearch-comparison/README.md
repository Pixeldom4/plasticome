# DIAMOND vs USEARCH — step 2/3 intermediates (2026-08-05)

Frozen copy of the clustering and alignment scratch for all four engine combinations, so
the comparison can be re-examined without re-running anything. Deliverables and summaries
stay in the live run directories; only the intermediates are archived here.

Every combination started from the same `01-union.tsv` (610 rows, md5
`beb17580d062b7ffd4a9fc1c669adb82`) at 90% identity, so the engine is the only variable.

## Folders

Named `clust-<engine>.align-<engine>` — the clusterer is step 2 (which sequences become
centroids), the aligner is step 3 (how those centroids group into components). Spelling out
both roles avoids the usual ambiguity about which half of an `a-x-b` name is which.

| folder | clustering | aligner | clusters | components | source run |
|---|---|---|---:|---:|---|
| `clust-usearch.align-usearch/` | USEARCH | USEARCH | 413 | 47 | `runs/2026-08-05.usearch/` |
| `clust-usearch.align-diamond/` | USEARCH | DIAMOND | 413 | 47 | `runs/2026-08-05.diamond-aligner-only/` |
| `clust-diamond.align-usearch/` | DIAMOND | USEARCH | 404 | 47 | `runs/2026-08-05.diamond-usearch-align/` |
| `clust-diamond.align-diamond/` | DIAMOND | DIAMOND | 404 | 48 | `runs/2026-08-05.diamond/` |

**The two off-diagonal runs never built their own cluster table** — that is what holds the
clustering factor fixed while the aligner changes. Their `02-clusters.intermediates/` is
therefore copied from the run that did build it: `clust-usearch.align-diamond/` takes it
from `runs/2026-08-05.usearch/`, and `clust-diamond.align-usearch/` from
`runs/2026-08-05.diamond/`. The cluster tables are byte-identical in each case
(`757aa257e2cbae4abe7462b1d9b06d02` and `5ebd7c772eeadd165efbce5af49e74bb`), so the
intermediates do describe how that folder's node set was produced.

## What is inside

`02-clusters.intermediates/` — the three clustering phases. Contents differ by engine:

- DIAMOND: `00-substrate.faa` / `.dmnd` and `00-edges.tsv` (the single all-vs-all edge graph
  every phase reuses), then `01-seed-*` and `03-denovo-*` edge subgraphs and
  greedy-vertex-cover output. Phase 2 leaves no file: recruitment is read straight off the
  edge graph in Python.
- USEARCH: `01-seed.*`, `02-closedref.*`, `03-denovo.*` — a FASTA, a `.uc` record file and a
  log per phase, since each phase re-searches instead of sharing a substrate.
- Both: `seeds.fasta`, `remaining.fasta`, `notmatched.fasta` (the engine-independent record
  of what each phase was handed) and `provenance.json` (engine, version, identity, coverage
  rule, greedy ordering, counts).

`03-alignment.intermediates/` — the node set and the partition. `combined.fasta` and
`combined_nodes.tsv` are the md5-unique node roster; `combined_pairs.tsv` holds the raw HSPs;
`component_{assignment,edges,members}_*.csv` are the partition; `stats_*.json` records the
aligner, its flags and every count. `combined_rev.fasta` exists in all four folders because
step 1 always writes it, but only the USEARCH arm searches it — `-allpairs_local` is
asymmetric in query/target, whereas DIAMOND reports both directions of a pair inside its one
run and so writes no `combined_rev_pairs.tsv`.

## Reading the four cells

The edge criterion is identical everywhere — ≥30% amino-acid identity **and** e-value < 1e-5,
applied as a post-filter, never as a search parameter — so the cells differ only in what the
engines fed it.

- **Aligner effect** (rows 1→2, 3→4): 42 of 47 components identical either way; six to seven
  nodes cross a boundary on near-threshold edges. DIAMOND passes ~120 more edges on both node
  sets (+121 at 413 nodes, +120 at 404) while reporting *fewer* raw aligned pairs.
- **Clusterer effect** (rows 1→3, 2→4): at the USEARCH aligner the component partition of all
  295 shared centroids is unchanged; at the DIAMOND aligner two nodes move. The clusterer's
  real effect is the centroid set — 404 vs 413, and which row of a duplicate pair labels it.
- Component sizes 2–5 are `45, 19, 16, 14` in all four runs. Only the largest component moves:
  233 → 229 → 224 → 218.

Full analysis: `runs/2026-08-05.diamond/REPORT_diamond_vs_usearch.md`, regenerable with
`runs/2026-08-05.diamond/compare_vs_usearch.py`.
