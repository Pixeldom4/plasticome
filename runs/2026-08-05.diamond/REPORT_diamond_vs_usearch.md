# DIAMOND vs USEARCH — 2026-08-05

Same `01-union.tsv` (610 rows), same code revision on every side; the only variable
is the engine.

| | clustering | alignment | run |
|---|---|---|---|
| **U** | usearch | usearch | `2026-08-05.usearch` |
| **A** | usearch (U's table) | diamond | `2026-08-05.diamond-aligner-only` |
| **D** | diamond | diamond | `2026-08-05.diamond` |

**U reproduces the pre-change pipeline exactly** — 413 clusters, 28 recruited, 47
components, 22,415 passing edges, largest component 233 — so it is a valid stand-in
for the previous result. (`2026-08-05.singletons-fixed` predates the seed-promotion
change and books 623 records under `S###` labels, so it is not label-comparable at
step 2; its step-3 numbers are identical to U's.)

## Step 2 — clustering

| | usearch | diamond | Δ |
|---|---:|---:|---:|
| clusters | 413 | 404 | -9 |
| reference (seed) clusters | 13 | 13 | +0 |
| de-novo clusters | 400 | 391 | -9 |
| singletons | 275 | 268 | -7 |
| largest cluster | 21 | 14 | -7 |
| records placed | 610 | 610 | +0 |

- size distribution, usearch: 1×275, 2×112, 3×15, 4×6, 5×3, 6×1, 21×1
- size distribution, diamond: 1×268, 2×107, 3×16, 4×6, 5×3, 6×1, 10×1, 11×1, 14×1

**Centroids.** 295 of usearch's 413 centroids are also diamond centroids (293 of them heading a byte-identical cluster). 118 usearch centroids are demoted to members, and diamond promotes 109 sequences usearch had kept as members.

**Records.** 563 of 610 records sit in a cluster with identical membership on both sides (the cluster may have swapped which member is its centroid).

**Partition.** 339 of 610 records keep the same centroid; 271 move. Over all 185,745 unordered pairs of records the two partitions agree on 185,546 (99.893%) — they disagree about co-clustering for 199 pairs, which is the honest size of the change: most 'moves' are a cluster keeping its membership while swapping which member is called the centroid.

**Seed clusters** (`cluster_id` 1–13, the curated anchors). Phase 2 is where the
coverage rule bites — diamond recruits 17 where usearch recruits 28:

| # | seed centroid | enzyme | usearch | diamond | released by diamond |
|---:|---|---|---:|---:|---|
| 1 | `U0001|WP_054022242.1` | ISPETase | 2 | 2 | — |
| 2 | `U0016|AEV21261.1` | LCC | 3 | 3 | — |
| 3 | `U0024|CBY05530` | TfCut_2 | 21 | 10 | `U0019|CAH17553.1`, `U0020|CAH17554.1`, `U0022|WP_011291330.1`, `U0023|AAZ54920.1`, `U0025|AET05798.1`, `U0026|AET05799.1`, `U0134|WP_104613137.1`, `U0480|Q6A0I4.1` (+3 more) |
| 4 | `U0018|WP_015787089.1` | Cut190 | 2 | 2 | — |
| 5 | `U0043|A0A075B5G4` | HiC | 2 | 2 | — |
| 6 | `U0239|SDZ16714.1` | Mipa-P | 1 | 1 | — |
| 7 | `U0245|WP_116180173.1` | Kubu-P | 1 | 1 | — |
| 8 | `U0051|WP_012476897.1` | NylA | 4 | 4 | — |
| 9 | `U0052|CAA24927.1` | NylB | 2 | 2 | — |
| 10 | `U0053|WP_012476895.1` | NylC | 6 | 6 | — |
| 11 | `U0085|WP_011614907.1` | PhaZ1 | 2 | 2 | — |
| 12 | `U0059|WP_037320891.1` | PLAase I | 1 | 1 | — |
| 13 | `U0048|WP_011061486.1` | PueA | 3 | 3 | — |

Seed clusters hold 50 records under usearch and 39 under diamond: 11 released into de-novo clusters. All 13 seeds stay centroids under both.

**De-novo clusters that change most** (shared centroids only):

- `U0419|MTK05241.1` 1 → 2  — left: —

## Step 3 — alignment

| | U usearch/usearch | A usearch/diamond | D diamond/diamond |
|---|---:|---:|---:|
| nodes aligned | 413 | 413 | 404 |
| raw alignments | 62,771 | 58,255 | 54,094 |
| aligned pairs | 30,859 | 29,343 | 27,259 |
| edges passing | 22,415 | 22,536 | 20,736 |
| components | 47 | 47 | 48 |
| single-node components | 26 | 25 | 25 |
| largest component | 233 | 229 | 218 |

- top-10 component sizes, U: [233, 45, 19, 16, 14, 10, 5, 5, 5, 5]
- top-10 component sizes, A: [229, 45, 19, 16, 14, 10, 5, 5, 5, 5]
- top-10 component sizes, D: [218, 45, 19, 16, 14, 10, 5, 5, 5, 5]

### Aligner isolated — U vs A (identical 413-centroid input)

Edges passing ≥30% aaid AND e<1e-5: 22,415 → 22,536 (+121).

413 centroids are nodes in both runs, in 47 components on the left and 47 on the right. 42 components are identical on both sides.

| from | to | n | nodes |
|---|---|---:|---|
| C001 | C001 (keeps identity) | 228 | — |
| C001 | C039 (keeps identity) | 4 | `U0040|RLI42440.1`, `U0117|noacc`, `U0118|noacc`, `U0122|noacc` |
| C034 | C036 (keeps identity) | 3 | `U0146|ACY56506.1`, `U0181|WP_034768559.1`, `U0576|WP_013956815.1` |
| C043 | C039 | 1 | `U0184|WP_034767800` |
| C037 | C036 | 1 | `U0064|WP_011157072.1` |
| C034 | C033 (keeps identity) | 1 | `U0454|WDT94443` |
| C033 | C001 | 1 | `U0458|noacc` |
| C001 | C040 (keeps identity) | 1 | `U0195|RGD93181.1` |

- **C001** absorbs 2 left-hand components (C001, C033) — 229 nodes.
- **C036** absorbs 2 left-hand components (C034, C037) — 4 nodes.
- **C001** (233 nodes) splits across C001, C039, C040.
- **C034** (4 nodes) splits across C033, C036.

### End to end — U vs D

Over the 295 centroids both runs have (usearch's other 118 are members, not nodes, under diamond, so their components cannot be compared directly):

295 centroids are nodes in both runs, in 33 components on the left and 36 on the right. 30 components are identical on both sides.

| from | to | n | nodes |
|---|---|---:|---|
| C001 | C001 (keeps identity) | 175 | — |
| C034 | C021 (keeps identity) | 2 | `U0146|ACY56506.1`, `U0576|WP_013956815.1` |
| C034 | C030 (keeps identity) | 1 | `U0454|WDT94443` |
| C033 | C028 (keeps identity) | 1 | `U0458|noacc` |
| C001 | C028 | 1 | `U0231|MSQ87502.1` |
| C001 | C024 (keeps identity) | 1 | `U0195|RGD93181.1` |
| C001 | C012 (keeps identity) | 1 | `U0040|RLI42440.1` |

- **C001** (178 nodes) splits across C001, C012, C024, C028.
- **C034** (3 nodes) splits across C021, C030.

## Bottom line

- Step 2 moves: 413 → 404 clusters. The mechanism is phase 2 — diamond's 90% member-coverage rule recruits 17 to the seeds where usearch's identity-only rule recruits 28 — plus degree-ordered rather than length-ordered centroid choice, which reshuffles which member of an unchanged cluster is called the representative.
- The partitions are near-identical in substance: they disagree about co-clustering for 199 of 185,745 pairs (0.107%).
- Step 3 is close to a no-op when isolated: same 47 components, +121 edges, and the giant α/β-hydrolase component shifting by a handful of nodes.
- Export artifact: `04-*.fasta` carries 413 records under usearch and 404 under diamond. These are different centroid sets and only one should feed the nr-search.
