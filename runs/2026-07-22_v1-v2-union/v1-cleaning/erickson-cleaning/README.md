# Erickson cleaning — `aa_sequence` reconciliation to Source Data D1

**File affected:** `plasticome.v1.1-erickson.tsv`
**Reference:** `../D1-erickson.csv` (Erickson et al., 2022 — Source Data Table D1: "Sequences for the 74 PET hydrolase candidates and the signal peptide-containing variants")
**Date:** 2026-07-22

## What was done

Every `aa_sequence` in `plasticome.v1.1-erickson.tsv` was **replaced with the verbatim "Protein Sequence" cell from D1**, so each row's sequence now matches the Erickson D1 reference exactly.

- **Rows changed: 21 of 21** (all rows).
- **Only the `aa_sequence` column was modified** — verified by diff against the backup. All other columns (including the curation-notes column) are untouched.
- The D1 protein sequences are the **expression constructs**: engineered N-terminal start-Met + mature domain (native signal peptide removed) + C-terminal `LEHHHHHH` His-tag. These were copied **byte-for-byte** (His-tag included), per request "so that it matches D1 exactly."

## How rows were matched — by `enzyme_name` number

Each row was mapped to the D1 entry whose **Enzyme ID equals the number in `enzyme_name`** (e.g. `Enzyme 403` → D1 #403, `Enzyme 101 like` → D1 #101).

| rowid | enzyme_name | → D1 # | old len | new len |
|------:|-------------|-------:|--------:|--------:|
| 103 | Enzyme 101 like | 101 | 613 | 310 |
| 109 | Enzyme 102 like | 102 | 565 | 305 |
| 3 | Enzyme 202 | 202 | 380 | 388 |
| 85 | Enzyme 204 like | 204 | 258 | 250 |
| 37 | Enzyme 211 like | 211 | 263 | 272 |
| 35 | Enzyme 214 like | 214 | 456 | 248 |
| 38 | Enzyme 301 like | 301 | 495 | 219 |
| 133 | Enzyme 305 like | 305 | 335 | 240 |
| 39 | Enzyme 307 like | 307 | 257 | 186 |
| 50 | Enzyme 403 | 403 | 269 | 292 |
| 73 | Enzyme 405 | 405 | 302 | 283 |
| 70 | Enzyme 406 | 406 | 292 | 313 |
| 83 | Enzyme 407 | 407 | 434 | 401 |
| 79 | Enzyme 409 | 409 | 283 | 277 |
| 41 | Enzyme 412 | 412 | 305 | 291 |
| 32 | Enzyme 503 like | 503 | 294 | 267 |
| 81 | Enzyme 607 | 607 | 310 | 271 |
| 78 | Enzyme 611 | 611 | 293 | 269 |
| 80 | Enzyme 701 | 701 | 301 | 270 |
| 8 | Enzyme 702 | 702 | 262 | 271 |
| 74 | Enzyme 711 | 711 | 284 | 269 |

## Why — background from the validation

Before the fix, each `aa_sequence` was the **native full-length NCBI protein** for the row's `retrieved` accession (verified byte-exact against NCBI for the spot-checked accessions). Those native sequences differ from D1 in a systematic way (native signal peptide present; no engineered Met; no His-tag), and validation surfaced two classes of discrepancy:

1. **"like" homologs** (`Enzyme 101/102/204/211/214/301/305/307 like`) — retrieved as BLAST homologs from *different* organisms, so their native sequence is only partially similar (6–99 % of the D1 core) to the D1 enzyme of the same number.
2. **400-series number offset** — four rows carried an accession whose native protein is actually a *different* D1 enzyme than the row's label (403↔409, 406↔504, 409↔412, 412↔406 by sequence identity).

In the input sheet these were flagged by hand: the **`retrieved` accession was blanked** and the notes column set to *"Accession removed bc. sequences didn't match"* on exactly the rows where the native accession sequence did **not** match D1 for that row's enzyme number. Rows whose native sequence *did* match D1 at the same number kept their accession.

The requested fix makes the `aa_sequence` authoritative to **D1 by enzyme number**, consistent with that curation: every row now holds the D1 sequence for its labeled enzyme.

## Change magnitude (before → after)

Full before/after sequences with per-row metrics are in the companion file
**`aa_sequence_before_after.tsv`** (columns: `rowid, enzyme_name, len_before, len_after,
delta_len, pct_similarity, aa_before, aa_after`).

"similarity" is the global sequence-similarity (Python `difflib` ratio) between the old and
new `aa_sequence`; "magnitude" is its complement.

| rowid | enzyme_name | len before | len after | Δlen | similarity | magnitude |
|------:|-------------|-----------:|----------:|-----:|-----------:|-----------|
| 103 | Enzyme 101 like | 613 | 310 | −303 | 51% | 49% changed |
| 109 | Enzyme 102 like | 565 | 305 | −260 | 52% | 49% changed |
| 3 | Enzyme 202 | 380 | 388 | +8 | 99% | 1% changed |
| 85 | Enzyme 204 like | 258 | 250 | −8 | 93% | 7% changed |
| 37 | Enzyme 211 like | 263 | 272 | +9 | 98% | 2% changed |
| 35 | Enzyme 214 like | 456 | 248 | −208 | 68% | 32% changed |
| 38 | Enzyme 301 like | 495 | 219 | −276 | 53% | 47% changed |
| 133 | Enzyme 305 like | 335 | 240 | −95 | 61% | 39% changed |
| 39 | Enzyme 307 like | 257 | 186 | −71 | 75% | 25% changed |
| 50 | Enzyme 403 | 269 | 292 | +23 | 20% | 80% changed |
| 73 | Enzyme 405 | 302 | 283 | −19 | 94% | 6% changed |
| 70 | Enzyme 406 | 292 | 313 | +21 | 39% | 61% changed |
| 83 | Enzyme 407 | 434 | 401 | −33 | 94% | 6% changed |
| 79 | Enzyme 409 | 283 | 277 | −6 | 82% | 18% changed |
| 41 | Enzyme 412 | 305 | 291 | −14 | 14% | 86% changed |
| 32 | Enzyme 503 like | 294 | 267 | −27 | 92% | 8% changed |
| 81 | Enzyme 607 | 310 | 271 | −39 | 91% | 9% changed |
| 78 | Enzyme 611 | 293 | 269 | −24 | 93% | 7% changed |
| 80 | Enzyme 701 | 301 | 270 | −31 | 92% | 8% changed |
| 8 | Enzyme 702 | 262 | 271 | +9 | 98% | 2% changed |
| 74 | Enzyme 711 | 284 | 269 | −15 | 94% | 6% changed |

The changes fall into three tiers:

- **Small (≤ 9 %, 11 rows)** — 202, 204, 211, 405, 407, 503, 607, 611, 701, 702, 711.
  Only the flanks moved: the native N-terminal signal peptide was trimmed to an engineered
  start-Met and `LEHHHHHH` was appended. Same protein / catalytic core. (Rows 702 and 211
  actually *gained* residues because their native form was already mature — only the Met +
  His-tag were added.)
- **Moderate (18–49 %, 6 rows)** — 101, 102, 214, 301, 305, 307 ("like" homologs).
  Genuinely different proteins: the previous native sequence was a BLAST homolog from another
  organism, now replaced by the actual D1 reference enzyme.
- **Large (60–86 %, 3 rows, plus 409 at 18 %)** — 403, 406, 412 (and 409). Near-complete
  replacement, because the old accession's protein was a *different enzyme entirely*. The
  sequences trade places in a cross-pattern (e.g. old #412 ≈ new #406, old #406 ≈ new #504),
  which is exactly the mislabeling this fix corrected.

## Verification

- All 21 sequences are **byte-identical** to their D1 source cell.
- All 21 end in the `...LEHHHHHH` His-tag.
- File structure intact: **22 columns** on every row; notes column preserved.
- Diff vs. backup confirms **`aa_sequence` is the only column that changed**.

## Backups

| file | contents |
|------|----------|
| `plasticome.v1.1-erickson.tsv.bak2` | the 21-row sheet **before** this replacement (native sequences) |
| `plasticome.v1.1-erickson.tsv.bak` | an earlier 14-row version (pre-dates the "complete set" update) |

## Known caveat — metadata not reconciled

The **sequences** now match D1, but some **labels/accessions still describe a different protein**:

- The eight `... like` rows keep the "like" qualifier though their sequence is now the exact D1 reference enzyme.
- On rows where the accession was **kept** (e.g. `Enzyme 503 like` → `EGD44994.1`, and the 405/407/607/611/701/702/711 group), the `retrieved`/`uniprot`/`pazy_accession` fields point at the native homolog, whose sequence is close to but not identical to the D1 construct now stored.

If a fully self-consistent table is needed, the `enzyme_name` "like" qualifiers and the residual accession fields should be reviewed against the new sequences.
