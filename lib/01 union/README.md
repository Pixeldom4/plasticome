# Plasticome v1.1 – v260701 Union

Merges the plasticome **v1.1** and **v260701** enzyme tables into one
accession-indexed table. v260701 is the canonical base; v1.1 is folded in.

## Inputs

| Role | File | Rows | Key columns |
|------|------|------|-------------|
| **v260701** | `source-data/plasticome.v260701/cleaned_pazy-260701_retrieving_from_accession.tsv` | 473 | `accession`, `enzyme_name`, `pazy_id`, `aa_sequence` |
| **v1.1** | `source-data/plasticome.v1.1/plasticome.v1.1.tsv` | 212 | `retrieved` (accession), `enzyme_name`, `component`, `aa_sequence`, + curated fields |

v260701 is the newer canonical table (PAZy snapshot `260701` = 2026-07-01); v1.1
is the cleaned form of the older curated table.

## Output

Six columns, one row per input record (see [Duplicates](#duplicates-are-intentional)).
**609 rows** on the current default inputs (the v260701 singletons table); 607 on
the earlier `cleaned_pazy-260701_retrieving_from_accession.tsv`.

| Column | Value |
|--------|-------|
| `plasticome_id` | unique 1..N, ordered on `pazy_id` (primary) then `component` (secondary) |
| `enzyme_name` | v260701 for v260701/both rows; v1.1 for v1.1-only rows |
| `accession` | v260701 `accession`; v1.1-only rows use v1.1 `retrieved` (blank if empty or a placeholder) |
| `pazy_id` | v260701 `pazy_id`; blank for v1.1-only rows |
| `aa_sequence` | see [Sequences](#sequences) — **607/607 populated** |
| `source` | `v260701`, `both`, or `v1.1` |

All other v1.1 curated columns (`substrate`, `EC`, `kingdom`, `reference`,
`pazy_accession`, `pdb_accession`, `uniprot_accession`, …) are dropped here but
remain in `plasticome.v1.1.tsv`.

### Ordering

`pazy_id` exists only on v260701-derived rows and `component` only on
v1.1-derived rows, so the documented sort degenerates cleanly into: rows with a
`pazy_id` first in numeric `pazy_id` order (ids 1–468), then the rest in numeric
`component` order, ties broken by input order. The 5 v260701 rows with no
`pazy_id` have no `component` either and land at the end.

## Merge logic

**Join key:** a v260701 row and a v1.1 row are the same enzyme iff v260701
`accession` equals v1.1 `retrieved` (only that column; empty and `-` ignored).
Matching is 1:1.

Building on v260701 as the base:

- v260701 row **matched** by a v1.1 `retrieved` → `both` (v1.1 row consumed)
- v260701 row **unmatched** → `v260701`
- v1.1 row **unmatched** → appended as `v1.1`

| `source` | Rows | Meaning |
|----------|------|---------|
| `v260701` | 395 | v260701 accession not in v1.1 |
| `both` | 78 | v260701 accession matched a v1.1 `retrieved` |
| `v1.1` | 134 | v1.1-only, appended |
| **Total** | **607** | all 473 v260701 rows + 134 v1.1-only |

### Duplicates are intentional

The same enzyme can appear on **two** rows — once as `v260701`/`both` and once
as `v1.1` — when its two records don't share an accession to join on. The
Erickson enzymes (`Enzyme 1xx`) and `jmPE13` are the main cases: their v1.1
accessions are blank or placeholders, so they don't merge. (`jmPE13` is in both
tables with no accession on either side, so it appears twice.) This is **by
design** — we do not deduplicate at this stage. Deduplication happens later, at
the **clustering** step. For the same reason, accessions appearing twice within
v260701 are kept as-is.

## Sequences

> **The rule:** sequences are always the ones retrieved from accession, never
> carried over from v1.1 — except where manually assigned during v1 cleaning.

| Row kind | Sequence source | Rows |
|---|---|---|
| `v260701` / `both` | v260701's `aa_sequence` (already retrieved from accession) | 469 |
| `both`, manually assigned | v1.1's — the manual assignment outranks the database record | 4 |
| `v1.1`, manually assigned | v1.1's | 18 |
| `v1.1`, everything else | retrieved from accession by Stage 2 | 116 |

**Manual assignments (22 rows).** Selected by `is_manual()` in
[`build_union.py`](build_union.py): the Erickson primary/only enzymes, named
`Enzyme <n>` / `Enzyme <n> like`, whose sequences were set to Erickson
supplementary table D1 verbatim; plus `jmPE13`, carried from its paper
supplement (no database record exists). Rows where Erickson is a *secondary*
reference keep their primary name (`Est1; Enzyme 708`, `MtCut; Enzyme 606`,
`RgCut-II`) and are **not** manual — hence a `startswith` test, not a substring
match.

Four of the 22 merge into v260701, and for those the manual sequence wins over
the database record — this is the whole point of the exception:

| Enzyme | accession | v1.1 (kept) | v260701 (discarded) |
|---|---|---|---|
| Enzyme 503 like | EGD44994.1 | 267 | 294 |
| Enzyme 611 | WP_093412886.1 | 269 | 293 |
| Enzyme 701 | WP_104613137.1 | 270 | 301 |
| Enzyme 702 | ADM47605.1 | 271 | 262 |

The deltas are the signal sequences / His-tags that the v1 cleaning step
deliberately normalized away.

**Stage 2 — accession lookup ([`fetch_sequences.py`](fetch_sequences.py)).** The
116 non-manual v1.1-only rows are left **blank** by Stage 1 specifically so this
script fills them from their accession. Cascade, first hit wins: **NCBI**
protein (65) → **UniProt** (+45) → **RCSB/PDB** (+2) → **UniParc** (+4). The
UniParc step recovers entries deleted from UniProtKB but archived in UniParc.
**0 unresolved.** Per-accession provenance is written to
`fetched_sequences.json` next to the union table.

> Stage 2 makes live HTTP calls. Set `NCBI_API_KEY` to raise the NCBI rate limit
> (3→10 req/s); it also runs without one. Results are stable (keyed on fixed
> accessions). It only ever writes into blank cells, so it is idempotent and
> cannot clobber a manual assignment.

⚠️ **Do not "fix" Stage 1 by carrying v1.1's sequence forward for all v1.1-only
rows.** `plasticome.v1.1.tsv` has `aa_sequence` populated on all 212 rows, so a
blanket carry-forward silently turns Stage 2 into a no-op and voids the
never-carried-over guarantee. (This was the behaviour before 2026-07-30.)

## Reproducing

```bash
python3 lib/union/build_union.py \
    source-data/plasticome.v1.1/plasticome.v1.1.tsv \
    source-data/plasticome.v260701/cleaned_pazy-260701_retrieving_from_accession.tsv \
    -o "runs/<run>/01 plasticome_v1.v260701-union.tsv"          # Stage 1, offline

python3 lib/union/fetch_sequences.py "runs/<run>/01 plasticome_v1.v260701-union.tsv"
                                                                 # Stage 2, networked
```

`build_union.py` rewrites the table from scratch each run and is deterministic.
`fetch_sequences.py` fills blanks in place and is idempotent.
