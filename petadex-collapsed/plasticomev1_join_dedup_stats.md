# PlastiComeV1 Join — Deduplication & Column Stats

**File:** `plasticomev1_join_retrievingfromaccession.tsv`
**Generated:** 2026-07-21

## Summary

**Duplication is defined solely by amino acid sequence.** By that criterion the
dataset has **0 duplicates** — all 501 sequences are unique, so no records are
redundant and a dedup pass removes nothing. Repeats in `accession`,
`enzyme_name`, or `organism` do **not** count as duplicates.

| Metric | Value |
|---|---|
| Data records | 501 |
| Columns | `identifier`, `enzyme_name`, `pazy_id`, `accession`, `organism`, `aa_sequence` |
| Identifier range | 1–501 |

> **Note:** The last line has no trailing newline, so a naive `wc -l` reports
> 500 data rows. Counting with `awk` (which counts the unterminated final
> line) gives the correct **501**.

## Uniqueness by column

| Column | Unique values | Notes |
|---|---|---|
| `aa_sequence` | 501 / 501 | Fully unique — **0 duplicates** (the only dedup criterion) |
| `accession` | 472 / 501 | Repeats do **not** count as duplicates (see below) |
| `enzyme_name` | 490 / 501 | Repeats do **not** count as duplicates |

## Accession breakdown

The 472-vs-501 gap in the `accession` column is **not** duplicate accessions:

- **29 rows have an empty `accession` field.**
- Only **1 accession appears on two rows:** `ALS54749.1`
  - id 19 — `Chath_Est1` — *Hungatella hathewayi*
  - id 490 — `Chath_Est1` — *Clostridium hathewayi DSM-13479*
  - The two organism names are synonyms, but the **two amino acid sequences
    differ** — so these are **not duplicates** under the sequence-only
    criterion, and both are correctly retained.

The `accession` field is compound (semicolon-joined, e.g.
`WP_054022242.1;A0A0K8P6T7;5XFY`); the counts above compare the full string.

## Enzyme-name repeats (top)

| Count | enzyme_name |
|---|---|
| 4 | `PURase;Amidase` |
| 2 | `jmPE14`, `PhaZAfa`, `PhaZ5`, `Lcp2`, `Lcp1`, … |

These share a name but correspond to distinct sequences/organisms.

## Takeaways

- ✅ **0 duplicates** by the sequence-only criterion — all 501 sequences unique;
  dedup is a no-op.
- ℹ️ 29 records lack an accession (does not affect duplicate status).
- ℹ️ `ALS54749.1` (id 19 / id 490) shares an accession + enzyme name but has two
  **distinct sequences** → not a duplicate; both retained.

## Reproduce

```bash
# record count (handles missing final newline)
awk 'NR>1' plasticomev1_join_retrievingfromaccession.tsv | awk 'END{print NR}'
# unique / duplicate amino acid sequences
awk -F'\t' 'NR>1{print $6}' plasticomev1_join_retrievingfromaccession.tsv | sort -u | wc -l
awk -F'\t' 'NR>1{print $6}' plasticomev1_join_retrievingfromaccession.tsv | sort | uniq -d | wc -l
# empty accessions
awk -F'\t' 'NR>1 && $4==""' plasticomev1_join_retrievingfromaccession.tsv | wc -l
```

---

# Change Log — Collapse of numbered enzymes

**Date:** 2026-07-21
**Input:** `plasticomev1_join_retrievingfromaccession.tsv` (501 records)
**Output:** `plasticomev1_collapsed.tsv` (489 records) — original left untouched
**Script:** `collapse_numbered_enzymes.py`

## Rationale

Some enzymes appear as **two separate rows** for the same catalytic entity:

1. a **bare-number** record whose `enzyme_name` is (or contains) the number
   itself — e.g. `101`, `202`, `RgCut-II;403`; and
2. an **`Enzyme N [like]`** record — e.g. `Enzyme 101 like`, `Enzyme 202`.

These were collapsed into a single row per number.

## Rules applied

- **Key = the enzyme number** (`Enzyme 101 like` → `101`, `Enzyme 202` → `202`).
- **The bare-number record is the truth.** The merged row keeps the
  bare-number record's `identifier` **and its `aa_sequence`**.
- **Amino acid sequences are NOT concatenated.** Only the truth sequence is
  kept; the `Enzyme N` record's sequence is **dropped**. (The two sequences
  genuinely differ — see the dropped lengths below.)
- **All other metadata** (`enzyme_name`, `pazy_id`, `accession`, `organism`)
  is merged with `;`, order-preserving, de-duplicated, empties dropped.

## Scope — what collapsed and what did not

- ✅ **12 numbers collapsed** (truth row + `Enzyme N` row → one row).
- ⏸️ **Not collapsed — number already in the same row** as its name
  (`Enzyme 405;405`, `Enzyme 503 like;503`, `606`, `611`, `701`, `702`, `711`):
  nothing to merge.
- ⏸️ **Not collapsed — no bare-number truth record exists**: `504` (id 83),
  `607` (id 464), `708` (id 234), and `406`/`409` (the stray `406`/`409` tokens
  live *inside* other `Enzyme N` names, they are not standalone records).
- ⚠️ **Borderline, left alone:** id 55 `607-Nsp` looks related to `Enzyme 607`
  (same organism) but its name is not *only* a number, so it was not treated as
  a truth record.

## The 12 collapses

| # | truth id (kept) | merged-in id (removed) | merged enzyme_name | kept seq len | dropped seq len |
|---|---|---|---|---|---|
| 101 | 35 | 491 | `101;Enzyme 101 like` | 310 | 613 |
| 102 | 36 | 486 | `102;Enzyme 102 like` | 305 | 565 |
| 202 | 37 | 465 | `202;Enzyme 202` | 388 | 380 |
| 204 | 38 | 482 | `204;Enzyme 204 like` | 250 | 258 |
| 211 | 39 | 481 | `211;Enzyme 211 like` | 272 | 263 |
| 214 | 40 | 466 | `214;Enzyme 214 like` | 248 | 456 |
| 301 | 41 | 489 | `301;Enzyme 301 like` | 219 | 495 |
| 305 | 42 | 496 | `305;Enzyme 305 like` | 240 | 335 |
| 307 | 44 | 487 | `307;Enzyme 307 like` | 186 | 257 |
| 403 | 45 | 49 | `RgCut-II;403;Enzyme 403;409` | 312 | 269 |
| 407 | 48 | 463 | `407;Enzyme 407` | 436 | 434 |
| 412 | 50 | 47 | `412;Enzyme 412;406` | 311 | 305 |

Notes on merged metadata:
- **403** (id 45) also absorbs pazy `141;145`, accessions
  `RLU00646.1;A0A3L8BW54;RLT92980.1;A0A3L8BDT3`, organisms
  `Ketobacter sp.;Rhizobacter gummiphilus` (the shared `Ketobacter sp.` was
  de-duplicated). The stray `409` token from `Enzyme 403;409` is carried along
  in the name field.
- **407** (id 48): organisms kept as `Allorhizocola rhizosphaerae;Allorhizocoloa
  rhizosphaerae` — note the **spelling variant** ("Allorhizocoloa") in the
  source; not auto-corrected.
- **412** (id 50): merges pazy `146;143`, accessions `WP_158643351.1;ODU60407.1`.
  The stray `406` token is carried in the name field.

## Result / integrity

| Metric | Before | After |
|---|---|---|
| Records | 501 | 489 |
| Unique aa_sequences | 501 | 489 |
| Duplicate aa_sequences | 0 | 0 |

- No amino acid sequence was duplicated or merged; each output record still has
  a unique sequence.
- Source file uses **CRLF** line endings; the loader strips them, so output
  sequences contain no stray `\r`. Output is written with the same CRLF style.

## Reproduce

```bash
python3 collapse_numbered_enzymes.py
# verify
awk 'NR>1' plasticomev1_collapsed.tsv | awk 'END{print "records:",NR}'
awk -F'\t' 'NR>1{print $6}' plasticomev1_collapsed.tsv | sort -u | wc -l   # unique seqs
```
