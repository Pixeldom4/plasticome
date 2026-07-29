# Plasticome v1–v2 Union

Merges the plasticome **v1** and **v2** enzyme tables into one accession-indexed
table: [`plasticome.v1-v2-union.tsv`](plasticome.v1-v2-union.tsv). v2 is the
canonical base; v1 is folded in.

## Inputs

| Role | File | Rows | Key columns |
|------|------|------|-------------|
| **v2** | `cleaned_pazy-260701_retrieving_from_accession.tsv` | 473 | `accession`, `enzyme_name`, `pazy_id`, `aa_sequence` |
| **v1** | `plasticome.v1.1.tsv` | 213 | `retrieved` (accession), `enzyme_name`, `aa_sequence`, + curated fields |

v2 is the newer canonical table (PAZy snapshot `260701` = 2026-07-01); v1 is the
older curated table.

## Output

Five columns; one row per input record (see [Duplicates](#duplicates-are-intentional)).

| Column | Value |
|--------|-------|
| `enzyme_name` | v2 for v2/merged rows; v1 for v1-only rows |
| `accession` | v2 `accession`; v1-only rows use v1 `retrieved` (blank if empty) |
| `pazy_id` | v2 `pazy_id`; blank for v1-only rows |
| `aa_sequence` | see [Sequences](#sequences) — **611/611 populated** |
| `source` | `v2`, `merged`, or `v1` |

All other v1 curated columns (`substrate`, `EC`, `kingdom`, `reference`, `pazy_accession`,
`pdb_accession`, `uniprot_accession`, …) are dropped here but remain in `plasticome.v1.1.tsv`.

## Merge logic

**Join key:** a v2 row and a v1 row are the same enzyme iff v2 `accession` equals v1
`retrieved` (only that column; empty and `-` ignored). Matching is 1:1.

Building on v2 as the base:

- v2 row **matched** by a v1 `retrieved` → `merged` (kept in v2 form; v1 row consumed)
- v2 row **unmatched** → `v2`
- v1 row **unmatched** → appended as `v1`

| `source` | Rows | Meaning |
|----------|------|---------|
| `v2` | 398 | v2 accession not in v1 |
| `merged` | 75 | v2 accession matched a v1 `retrieved` |
| `v1` | 138 | v1-only, appended |
| **Total** | **611** | all 473 v2 rows + 138 v1-only |

### Duplicates are intentional

The same enzyme can appear on **two** rows — once as `v2`/`merged` and once as `v1` —
when its v2 and v1 records don't share an accession to join on. The Erickson enzymes
(`Enzyme 1xx`) and `jmPE13` are the main cases: their v1 accessions are blank or
placeholders, so they don't merge. This is **by design** — we do not deduplicate at
this stage. Deduplication happens later, at the **clustering** step of the pipeline.
(For the same reason, three accessions that appear twice within v2 —
`BAB86909`, `WP_252811421.1`, `ANG60415` — are simply kept as-is.)

## Sequences

Every row carries an `aa_sequence`, filled in two stages:

**Stage 1 — merge ([`build_union.py`](build_union.py)).** v2/merged rows take v2's
sequence. v1-only rows carry v1's sequence when v1 has one — **18 rows**, listed
below. For the 12 with no real accession, v1 is the *only* possible source (they can
never be database-fetched); `jmPE13`'s sequence was curated by hand into
`plasticome.v1.1.tsv`.

| Enzyme | v1 accession | Len |
|--------|--------------|-----|
| Enzyme 101 like | *(none)* | 310 |
| Enzyme 102 like | *(none)* | 305 |
| Enzyme 202 | 7QJM | 388 |
| Enzyme 204 like | *(none)* | 250 |
| Enzyme 211 like | *(none)* | 272 |
| Enzyme 214 like | *(none)* | 248 |
| Enzyme 301 like | *(none)* | 219 |
| Enzyme 305 like | *(none)* | 240 |
| Enzyme 307 like | *(none)* | 186 |
| Enzyme 403 | *(none)* | 292 |
| Enzyme 405 | WP_082414832.1 | 283 |
| Enzyme 406 | *(none)* | 313 |
| Enzyme 407 | WP_117215036.1 | 401 |
| Enzyme 409 | *(none)* | 277 |
| Enzyme 412 | *(none)* | 291 |
| Enzyme 607 | WP_107095481.1 | 271 |
| Enzyme 711 | WP_083947829.1 | 269 |
| jmPE13 | *(placeholder)* | 273 |

The first 17 are the **Erickson** functional-metagenomics enzymes.

**Stage 2 — accession lookup ([`fetch_sequences.py`](fetch_sequences.py)).** The other
**120 v1-only rows** were blank in v1 but each has a real accession. They are filled by
a cascade, first hit wins: **NCBI** protein (68) → **UniProt** (+46) → **RCSB/PDB**
(+2) → **UniParc** (+4). The UniParc step recovers four entries deleted from UniProtKB
but archived in UniParc (`A0A0P0ZE81`, `A0A291HVH1`, `E5BBQ2`, `P94146`). **0 unresolved.**
Per-accession provenance is written to [`fetched_sequences.json`](fetched_sequences.json).

> Stage 2 makes live HTTP calls. Set `NCBI_API_KEY` to raise the NCBI rate limit
> (3→10 req/s); it also runs without one. Results are stable (keyed on fixed accessions).

## Reproducing

```bash
python3 build_union.py      # Stage 1: merge, carry v1 sequences forward (offline)
python3 fetch_sequences.py  # Stage 2: fill remaining rows by accession (networked)
```

`build_union.py` rewrites the table from scratch each run. `fetch_sequences.py` fills
blanks in place and is idempotent (skips rows that already have a sequence).
