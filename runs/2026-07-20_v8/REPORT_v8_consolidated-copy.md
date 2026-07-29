---
title: PETadex v8 — Component Assignment (curated 484-node PAZy TSV) — Results
date: 2026-07-20
input: petadex-alignment/data/cleaned_pazy_final_v8-align.tsv (snapshot v8/cleaned_pazy_final.tsv, md5 1d13e83691cc767c7ba9635c9bd2ed60)
v1_overlay: petadex-alignment/v8/plasticome_v1.csv (213 rows / 205 unique md5)
aligner: usearch v11.0.667 -allpairs_local (Docker linux/amd64 on arm64 macOS)
partition_mode: flat de novo (no v1 anchor)
tags: [petadex, plasticome, components, graph, results, v8]
---

> **Consolidation note.** This is the scientific narrative for the run reproduced by the
> consolidated `analysis/` folder (see `README.md`). The `analysis/` pipeline is byte-for-byte
> identical to the v8 run described here — same 484 nodes, same 46 components, same C-labels on
> all 484 nodes. Deliverable filenames below correspond to `analysis/outputs/*_petadex_<date>.csv`.

# PETadex v8 — Results (flat 484-node partition, curated/re-sequenced PAZy)

Re-runs the [[PETadex v2 — Component Assignment Workflow]] Steps 0–3 on the **curated** PAZy TSV
(`data/cleaned_pazy_final_v8-align.tsv`, snapshotted into `v8/`). Same paper-faithful method as
v4–v7: usearch v11.0.667 `-allpairs_local`, permissive search + `≥30% aaid AND e<1e−5`
post-filter, single-linkage connected components, both FASTA orientations. Flat de-novo, no v1
anchor. The input now carries a pre-existing `component_id` column (prior assignment); it is
**preserved for comparison only** (`prior_component_id`) and never seeds or gates the partition.

## What changed vs v7

| | v7 (`cleaned_pazy_final.tsv`, 17:21) | v8 (curated, 14:38) |
|---|---|---|
| Source rows | 491 (1 md5-dup pair) | **484** |
| Nodes (md5-unique) | 490 | **484** (−6) |
| md5-duplicate rows | 1 pair (140=141) | **0** (every row a distinct sequence) |
| Search space N (residues) | 171,118 | **169,346** |

**All node changes are curation re-sequencings, not biological drops.** Nine v7 nodes carry a new
sequence in v8 for the *same* enzyme/accession — including three headline fixes the v7 report had
flagged:

- **`porcine pancreatic lipase`** (id 453): the 112-aa fragment flagged as a G3 short-node in v7 is
  now the **full 448 aa** (still its own singleton, C045 — the correction did not spuriously merge it).
- **`EstB`** (id 15): 218 → **231 aa** (lands in C003, its prior component).
- **`GuaPA`** (id 444): re-sequenced to **324 aa** (C043, its prior component).
- Plus re-sequencings of `Thc_Cut1`, `PHL-7`, `Chath_Est1`, `TfH`, `LCP_STRK3`, `ESM136`.

The node count falls 490 → 484 because six of these corrected sequences now coincide with
sequences already in the set (de-duplicated by the fix) while three become newly distinct. **Zero
sequences were dropped from the biology.**

## The partition is v7, restricted to the shared sequences

Md5-keyed comparison against v7: the two partitions are **identical on all 481 shared sequences**
(43 vs 43 restricted clusters, exact match). Every corrected sequence slots into the component it
belonged to before. So v8 reproduces v7's clustering exactly on everything they share and simply
carries the curation edits.

**Against the input's own prior `component_id`:** a **perfect bijective relabeling** — 0 prior
components split across >1 v8 component, 0 v8 components absorbing >1 prior component. 471/484 rows
carry the identical label; the 13 differences are a single **C031 ↔ C032 ID swap** (canonical
order shifted by the re-sequencing), with membership fully preserved.

## Headline counts

| quantity | value |
|---|---|
| Source rows | 484 |
| Nodes (md5-unique) | **484** |
| Search space N (residues) | **169,346** |
| Aligned pairs (both orientations, best HSP) | 44,659 |
| **Edges passing** ≥30% aaid AND e<1e−5 @ 484-db | **33,101** |
| **Components** | **46** |
| Singletons | 23 |
| Non-singleton components | 23 |
| Largest component | **283** (58% of all nodes) |

E-values thresholded as reported at the 484-node db (decision 2: no down-scaling to the paper's
213-space), on the max-bits + both-orientation best HSP per pair (B4 convention, adopted v3 §7).

## Largest components (v1 overlay is descriptive only)

| rank | comp | n | w/ v1 | canonical (id) | v1 labels | CATH (overlay) |
|---|---|---|---|---|---|---|
| 1 | C001 | 283 | 74 | 1 | 1;20;35 | 3.40.50.1820 (α/β-hydrolase) |
| 2 | C010 | 48 | 6 | 47 | 9 | 3.60.70.12 |
| 3 | C008 | 25 | 17 | 44 | 10 | 3.90.1300.10 |
| 4 | C006 | 18 | 16 | 36 | 2 | 3.40.50.1820 |
| 5 | C004 | 15 | 11 | 30 | 23 | 3.40.50.1820 |
| 6 | C032 | 10 | 5 | 105 | 16 | NA |
| 7 | C011 | 8 | 5 | 52 | 21 | 3.40.50.1820 |
| 8 | C021 | 7 | 6 | 70 | 25 | 3.40.50.1820 |
| 9 | C007 | 5 | 5 | 39 | 4 | 3.40.50.1820 |
| 10 | C014 | 5 | 3 | 56 | 26 | 3.40.50.1820 |

**The α/β-hydrolase superfamily still dominates.** As in v2–v7, one component (C001) absorbs
283/484 nodes — the expected single-linkage merge of the broad α/β-hydrolase fold; its v1 overlay
spans labels 1, 20, 35 and it carries CATH 3.40.50.1820. C001 shed 4 nodes vs v7's 287 (the
re-sequencings that de-duplicated onto existing members).

## v1-overlay coherence (report, don't enforce)

- **191 of 484** nodes carry a v1 overlay.
- **Zero** v1 labels are split across >1 v8 component — every v1 family maps to exactly one v8
  component, identical to v5/v6/v7.

## Fragment / substring caveats (G3)

Per decision 3 no fragments were filtered. The fragment surface **shrank** vs v7:

- The old 112-aa `porcine pancreatic lipase` fragment is gone (corrected to 448 aa). Only **one**
  node is now <100 aa: the **28-aa** `Arylacylamidase` (id 50), which again attaches into C008
  (rank 3) via a short high-identity patch — the same documented G3 short-patch behaviour as
  v4–v7 (edge real per the criterion, rests on ~28 aa).
- The next-shortest full nodes (176 aa `PLAase III` → C013; 186 aa nodes → C001/C006) sit in their
  expected components; no substring-contained node became a spurious singleton.

## Validation gates

| gate | result |
|---|---|
| md5 uniqueness into graph (484, none twice) | ✅ 484 unique node_ids / 484 unique md5 |
| Every graphed node in exactly one component | ✅ 484/484 assigned; 0 blank-seq rows |
| Partition equivalence to v7 (md5-keyed clustering) | ✅ identical on all 481 shared sequences |
| Bijection to input's prior `component_id` | ✅ 0 splits / 0 merges (13 rows = C031↔C032 relabel) |
| **Step 0 engine sanity** (213 v1 nodes → 42 comps) | ✅ **42 components / 3,178 edges** (matches v4–v7 & paper's 42) |

## Deliverables (`petadex-alignment/v8/`)

- `component_edges_v8_2026-07-20.csv` — 33,101 passing edges.
- `component_assignment_v8_2026-07-20.csv` — 484 rows, one per node (`component_id`, `size_rank`,
  identifier(s), accession(s), pazy_id, v1/cath overlay).
- `component_members_v8_2026-07-20.csv` — 46 rows, one per component (canonical label, size,
  v1-label spread, CATH folds, member identifiers).
- `cleaned_pazy_final_v8_components_2026-07-20.csv` — all 484 source rows in original order with
  their v8 component joined, plus `prior_component_id` and `node_status`.
- `combined{,_rev}.fasta`, `combined{,_rev}_pairs.tsv` — node set + raw usearch HSPs.
- `combined_nodes.tsv`, `combined_stats.json`, `stats_v8_2026-07-20.json` — rosters + provenance.
- `cleaned_pazy_final.tsv` — frozen input snapshot (md5 `1d13e83691cc767c7ba9635c9bd2ed60`).
- `step0_check/` — the 213 v1-node engine-sanity inputs, pairs, and partition (→ 42).
- `step0_v8.py`, `step1_nodes_v8.py`, `step23_graph_v8.py`, `annotate_source_v8.py`, `run_v8.sh`
  — pipeline (`step23_graph_v5.py` retained for the accession-labeled Step 0 sanity roster).

## What this run does not claim

- No validation against the paper's 42 (flat mode; the 42 lives only in the Step 0 engine check).
  No join/merge/novel typing. No domain models (Step 4 deferred — HMMER absent).
- The **B4** max-bits convention and its accepted costs (near-30% boundary edges hidden behind a
  higher-bits sub-30% HSP; a few tie-sensitive edges) still apply, carried forward from v3–v7.
- The pre-existing `component_id` in the input is treated as a prior/overlay only; the v8 partition
  is recomputed de-novo and happens to reproduce it exactly (bijective relabel).
