# Step 0 graphs for Cytoscape

All graphs share **one node table** (`nodes.csv`, 213 nodes) so they overlay. Import that
first, then any edge table.

## Import into Cytoscape
1. **File → Import → Network from File →** an `edges_*.csv`. `source` = Source Node,
   `target` = Target Node; the rest (`pct_id`, `evalue`, `bits`, `presence`) are edge attributes.
2. **File → Import → Table from File →** `nodes.csv`. Key column = `node_id`.
3. Style: color/group nodes by `v1_component` (the paper's 42 labels), `cath` (fold), or a
   recomputed `comp_*` column. `node_id` = `n<rowid>`; `is_redundant`/`md5_group` mark the 8
   duplicate-sequence nodes.

## Graphs (all: ≥30% aa-identity AND e<1e-5, undirected single-linkage)

| edge file | aligner / command | edges | components |
|---|---|---:|---:|
| `edges_allpairs_u11.csv` | **usearch v11.0.667 `-allpairs_local`** — paper-faithful | 3,170 | **42** ✅ |
| `edges_usearch_local_u11.csv` | usearch v11.0.667 `-usearch_local` | 3,173 | 40 |
| `edges_usearch_local_u12.csv` | usearch v12.0-beta `-usearch_local` | 3,155 | 46 |
| `edges_usearch_local_v11_vs_v12.csv` | union of the two `usearch_local` graphs, `presence` = `both`/`v11_only`/`v12_only` | 3,179 | — |

- **`edges_allpairs_u11.csv` is the reproduction of the Logan paper's PAZy network (42 comps).**
  Use this as the reference graph.
- **`edges_usearch_local_v11_vs_v12.csv` is the version-difference view.** Color edges by
  `presence`: **both = 3,149**, **v11_only = 24**, **v12_only = 6**. The difference is small and
  is an e-value **normalization** effect (v12's e-values run 10^4.13 higher), *not* different
  alignments — identity and bit-scores match (outline §2R).

## HSP-collapse rule (read this — it explains the component counts)
usearch reports **multiple local HSPs per pair**. How you collapse them to one edge changes the
count, so the rule is stated per graph:

- **`allpairs_local` uses Rule A** (represent each pair by its highest-identity HSP; keep if that
  HSP has ≥30% id AND e<1e-5). Only 209/5,658 pairs are multi-HSP here, and Rule A reproduces the
  paper's **42** (any-HSP / best-bits give 41 — a one-edge boundary wobble, same magnitude as the
  2026-vs-2024 PAZy snapshot difference).
- **`usearch_local` graphs use Rule C** (edge if **any** HSP has ≥30% id AND e<1e-5 — the literal
  reading of "retain significant alignments"). `usearch_local` emits ~17,600 multi-HSP pairs, and
  Rule A there is unfair: its max-identity HSP is often a short spurious high-id patch that fails
  the e-value, dropping strong pairs (e.g. n1↔n2 at 47.7% id / e=5e-62). Rule C is immune to this
  and gives a fair v11-vs-v12 comparison.

Bottom line: the **42-component reproduction (allpairs) is the stable result**; `usearch_local`
component counts are HSP-rule- and version-sensitive (40–48 across rules/versions) precisely
because heuristic local search is not the paper's exhaustive all-vs-all.
