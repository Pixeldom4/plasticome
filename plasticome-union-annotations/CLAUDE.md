# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project state

This directory is for generating figures from the **PETadex plasticome v2** dataset (plastic-degrading enzymes). It lives inside the `plasticome` pipeline repo (see the root `README.md`), but it's standalone. It holds loose Python scripts run from this directory, with no build, test or lint setup. `pipeline_run.bash` doesn't call anything here.

**Name collision:** this directory's `annotate_union.py` is not `lib/05 annotate/annotate_union.py`. Step 5 adds cluster and component assignments to the union; this script adds substrate, organism and PDB metadata. The inputs are the same union, though: `plasticome-v1-v260701-union.csv` has the same 609 rows as `runs/2026-08-06.final-usearch.3/01-union.tsv`, with the same `plasticome_id`, accession, `pazy_id`, source and sequence.

Environment: Homebrew Python 3.14 with `pandas` and `matplotlib`. **`openpyxl` is not installed**, so any `pd.read_excel` call (including in `annotate_union.py`) fails until you run `pip3 install openpyxl`.

## Pipeline

```
python3 annotate_union.py   # union.csv + workbook  →  plasticome-v1-v260701-union.annotated.corrected.csv
```

`annotate_union.py` joins each union row back to its source table(s) according to `source`:

- **v1.1 rows:** it matches the accession against any of `pazy_accession`/`uniprot_accession`/`pdb_accession`/`BLAST_link`, with version suffixes stripped. It falls back to the md5 of the sequence.
- **v260701 rows:** it matches on `pazy_id`.

It adds `substrate, organism, superkingdom, metagenomic, uniprot_accession, pdb_accession, has_structure` before `aa_sequence`. Where the two sources disagree, PAZy wins. Substrates and IDs are unioned and joined with `"; "`. Superkingdom is normalised through the `SUPERKINGDOM` map. The output keeps the same two-row header convention (description row, then header row). The script prints fill rates and any rows it couldn't match; check those after changes.

PDB typos and obsolete IDs in the sources are corrected through the `PDB_FIXES` map. Each replacement was verified by aligning the RCSB sequence to the row's `aa_sequence`, so add new corrections there rather than editing the output. Rows with an empty `pazy_id` fall back to matching PAZy on `accession`.

Output history:

- `plasticome-v1-v260701-union.annotated.csv`: the earlier output, before the corrections.
- `plasticome.v2 - plasticome-v1-v260701-union-annotated.csv`: a Google Sheets export of that earlier output. It is identical apart from line endings (LF instead of CRLF).
- `*.annotated.corrected.csv`: the current output.

## Data files

### `plasticome.v2.xlsx`: the full curation workbook

Sheets appear roughly in pipeline order:

- `plasticome.v1`, `erickson-changes`, `plasticome.v1.1`, `plasticome.v1.1.no-seq`: the legacy v1 plasticome.
- `pazy-260701_pull`, `pazy-260701_manual-curation`, `cleaned_pazy-260701`, `cleaned_pazy-260701-singletons-*`: a fresh PAZy database pull (2026-07-01) and its manual cleanup.
- `singletons-validation.1` / `.2`: a review of connected components, including singletons.
- `plasticome-v1-v260701-union`: v1.1 and the cleaned PAZy pull merged and collapsed on accession, with no deduplication.
- `plasticome-union.non-redundant`: the union collapsed to one row per unique amino-acid sequence, using an md5 hash.
- `plasticome_centroids_final`: the union clustered at 90% AA identity, with one row per centroid and its CATH assignment(s).
- `cath_domain_annotations`, `plasticome.v2.sequence_annotati…`: per-CATH-fold / per-domain annotations, including domain names and catalytic residues.

**Parsing gotcha:** in several derived sheets (the union, non-redundant and centroids sheets), and in every CSV here, **row 1 is a free-text description and the real header is on row 2**. Use `header=1` in `pd.read_excel` / `pd.read_csv`; the raw source sheets (`plasticome.v1.1`, `pazy-260701_pull`) use a normal row-1 header. Read with `dtype=str, keep_default_na=False` so the `NA`/`-` sentinels stay as strings. Several sheets also have trailing empty columns, and `plasticome_centroids_final` has a duplicated `accession` column. Sheet names are truncated to Excel's 31-character limit.

### `plasticome-v1-v260701-union.csv`

A CSV export of the `plasticome-v1-v260701-union` sheet (~609 enzymes). Columns: `plasticome_id, enzyme_name, accession, doi, pazy_id, source, catalytic-cath, accessory-cath, domain, centroid, centroid_id, aa_sequence`.

- `source` ∈ {`v1.1`, `v260701`, `both`}: which dataset(s) the enzyme came from.
- `centroid` is `yes`/`no`. `centroid_id` links non-centroid members to their 90%-identity cluster representative (`NA` when the row is itself a centroid).
- `catalytic-cath` holds CATH superfamily IDs; `3.40.50.1820` (α/β hydrolase) dominates. `accessory-cath` uses `-` for none. Missing values are written as the literal strings `NA` or `-`, not left empty.

### `record_structure_map.tsv`

A 73 MB TSV with about 576k rows and a normal header. Columns: `rec_id, query_petadex_id, component_id, matched_accession, uniprot_ac, s3_key`. It maps each homolog record, keyed by plasticome ID (`PL###`) and component (`C###`), to an AlphaFold DB model stored under `s3://petadex-protein-structures/af_db_v260701/`. Nothing here reads it yet.
