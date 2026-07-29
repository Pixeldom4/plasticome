# DIAMOND all-vs-all alignment — experiment report

**Date:** 2026-07-21
**Folder:** `petadex-alignment-diamond/allvsall/`
**Engine:** DIAMOND v2.2.4 (Linux x86-64 ELF, run via Docker `linux/amd64` on arm64 macOS)
**Goal:** All-vs-all align the `substrings-collapsed` sequence set with settings as
close as possible to the project's usearch pipeline (v4–v8 method), then compare
the resulting single-linkage clustering against the usearch **v10** run on the
same input.

---

## 1. Input

| | |
|---|---|
| Source TSV | `../substrings-collapsed/plasticome-substrings-collapsed.tsv` |
| Sequences | **481** protein sequences (column `aa_sequence`, keyed by `identifier` 1…481) |
| FASTA built | `seqs.faa` (481 records, all ids unique, none empty) |
| DB letters | 168,497 |

FASTA header = the `identifier` column; body = `aa_sequence`. No sequences were
dropped or deduplicated at this stage (the input is already substring-collapsed).

---

## 2. Method — usearch → DIAMOND parameter mapping

The usearch reference method (from `petadex-alignment-usearch/analysis/config.py`
and `run.sh`) is:

```
usearch11 -allpairs_local <fasta> -userout <out> \
          -userfields query+target+id+evalue+bits -acceptall
```
permissive search, then a **post-filter of pct_id ≥ 30% AND e-value < 1e-5**,
single-linkage components, run in **both FASTA orientations**.

DIAMOND command used:

```bash
diamond blastp \
  --query seqs.faa --db seqs.dmnd \
  --out allpairs.tsv \
  --outfmt 6 qseqid sseqid pident evalue bitscore \
  --ultra-sensitive \
  --max-target-seqs 0 \
  --evalue 10 \
  --masking 0 \
  --comp-based-stats 0 \
  --threads 4
```

| usearch setting | DIAMOND setting | Rationale |
|---|---|---|
| `-allpairs_local` | `blastp`, query = db | Query against itself = all-vs-all; both are **local** (Smith-Waterman-style) alignment. |
| exhaustive allpairs | `--ultra-sensitive` | usearch aligns *every* pair exhaustively; DIAMOND is seed-based. `--ultra-sensitive` is its most exhaustive mode — the closest heuristic approximation. At 481 seqs it recovers essentially everything. |
| `-acceptall` | `--evalue 10` | Mirrors "accept everything at search time, filter later." 10 ≫ the 1e-5 gate, so the search never drops a candidate edge. The **post-filter is the single source of truth.** |
| reports all pairs | `--max-target-seqs 0` | **Critical.** DIAMOND defaults to top-25 hits/query, which would silently truncate a dense all-vs-all. `0` = unlimited. |
| `-userfields query+target+id+evalue+bits` | `--outfmt 6 qseqid sseqid pident evalue bitscore` | Same five columns, same order. |
| (no low-complexity masking) | `--masking 0` | DIAMOND 2.x masks by default (tantan); usearch allpairs does not. Disabling keeps identity/coverage comparable. |
| (plain SW scoring) | `--comp-based-stats 0` | DIAMOND defaults to composition-based e-value adjustment; usearch does not. Off = plainer model, closer to usearch. |
| default BLOSUM62, gap 11/1 | DIAMOND defaults | Both default to BLOSUM62, gapopen 11 / gapextend 1 — left as-is. |
| **run both FASTA orientations** | **single run** | The orientation-doubling is a usearch quirk: its allpairs is order-dependent and computes each pair once. DIAMOND's query-vs-db already evaluates A→B *and* B→A, so one run covers both orientations. |

**Post-filter + clustering** (`cluster.py`): keep an edge iff `pct_id ≥ 30.0 AND
evalue < 1e-5`, drop self-hits, then single-linkage (union-find) components.
Filtering is done post-hoc, never as a search parameter — identical to the
usearch pipeline.

### Reproduce

```bash
./run.sh            # from petadex-alignment-diamond/allvsall/
```
Steps: TSV→FASTA → `diamond makedb` → `diamond blastp` (above) → `cluster.py`.

---

## 3. Results

| Quantity | Value |
|---|---|
| Raw hits emitted | 83,248 |
| Edges passing filter (id ≥ 30% & e < 1e-5, self-hits dropped) | 32,622 unique undirected |
| Nodes in ≥ 1 edge | 458 |
| **Single-linkage components** | **46** |
| Singletons | 23 |
| Largest 8 component sizes | `281, 48, 25, 18, 15, 10, 8, 7` |

Outputs:
- `allpairs.tsv` — raw all-vs-all hits (5 columns)
- `component_assignment.csv` — `identifier, component, component_size`
- `seqs.faa`, `seqs.dmnd` — FASTA + DIAMOND database

---

## 4. Comparison vs usearch v10

Reference: `petadex-alignment-usearch/analysis/outputs/v10/component_assignment_v10_2026-07-21.csv`
(same 481 identifiers, joined on `identifier`).

| Metric | Value |
|---|---|
| Components | **46 (DIAMOND) vs 46 (usearch)** |
| Adjusted Rand Index | **0.9841** |
| Size spectrum | `[281,48,25,18,15,10,8,7…]` vs `[280,48,25,18,15,10,8,7…]` — identical except the top |
| Pairs co-clustered in both | 40,621 |
| Pairs together in DIAMOND only | 559 |
| Pairs together in usearch only | 282 |
| Sequences at the disagreement boundary | **4 of 481 (0.8%)** |
| Components matching 1:1 exactly | 42 of 46 |

> Note: a naive "differing co-membership" count returns 286 nodes, but that is a
> statistical artifact — because the giant cluster's roster changes by ±2 nodes,
> every one of its ~280 members trivially gets a different partner-set. The honest
> measure is the **4 boundary sequences** below; ARI (0.984) already accounts for
> this correctly.

### The 4 boundary sequences — every real difference

| id | enzyme | usearch | DIAMOND | direction |
|---|---|---|---|---|
| 164 | Ces19_14 (*C. thermoamylovorans*, len 240) | singleton `C038` | joins big cluster | DIAMOND **merges** |
| 434 | GuaPA (len 324) | singleton `C043` | joins big cluster | DIAMOND **merges** |
| 175 | HG-3 (*Clostridiales* AM23-16LB, len 270) | in big cluster `C001` | singleton | DIAMOND **splits off** |
| 430 | TA26 (*Thermoleophilum album*, len 321) | in `C034` (n=4) | singleton | DIAMOND **splits off** |

Net effect on the giant component: usearch 280 → DIAMOND 281 (+164, +434, −175).
`C034` loses 430 (4 → 3). All 42 other components (~477 sequences) are identical.

---

## 5. Interpretation & caveats

- **DIAMOND reproduces the usearch clustering to within 4 borderline sequences**
  (ARI 0.984, matching component counts). This is the expected residual from
  swapping alignment engines, not a methodology discrepancy.
- All 4 disagreements sit on the giant cluster's rim or in a small component —
  i.e. the *low-confidence* edges straddling the 30% / 1e-5 threshold, precisely
  where the engines' different alignment heuristics and e-value models tip an edge
  in or out.
- **E-values are not directly comparable across engines.** usearch and DIAMOND use
  different Karlin–Altschul parameterizations and effective DB sizes, so the same
  alignment gets a different e-value from each. The ≥ 30% identity term is the more
  transferable criterion; the 1e-5 gate is the softer of the two.
- Component labels differ cosmetically (DIAMOND `1,2,3…` vs usearch `C001,C034…`);
  components were matched by membership overlap, not by label.

## Files
```
allvsall/
├── run.sh                                  # end-to-end driver
├── cluster.py                              # post-filter + single-linkage
├── seqs.faa / seqs.dmnd                    # input FASTA + DIAMOND DB
├── allpairs.tsv                            # raw all-vs-all hits
├── component_assignment.csv                # 46-component clustering
└── REPORT_diamond_2026-07-21.md            # this file
```
