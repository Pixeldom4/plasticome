# Diff: `cleaned_pazy_final (2).csv` (A) vs `cleaned_pazy_final-tmp.csv` (B)

_Date: 2026-07-21_

| File | Rows | Columns |
|------|------|---------|
| **A** — `cleaned_pazy_final (2).csv` | 484 | 8 |
| **B** — `cleaned_pazy_final-tmp.csv` | 481 | 9 |

**B is the newer, re-curated version.** Diffing naively on the `accession` field reports ~150 rows changed on each side, but that is misleading — `accession` is a semicolon-joined list that was heavily expanded in B. Matching on `aa_sequence` (the true biological key) shows the datasets are **almost identical**, with B being an enriched/deduplicated version of A.

## Schema change

- **New column in B: `sequence_length`** — populated for all 481 rows (0 empty).

## Genuine row changes (net −3 rows)

Matching on protein sequence, **479 rows are shared**. The remaining differences are dedup/merges, **not real losses**.

**5 sequences in A absent from B — all collapsed into other rows in B:**

- `607-Nsp` (KPI31299.1) + `Enzyme 607` (WP_107095481.1) → merged into one row `WP_107095481.1;KPI31299.1`
- `406` (ODU60407.1) → folded into the `412` row (`WP_158643351.1;ODU60407.1;…`)
- `409` (RLT92980.1) → folded into the `403/409` row (`RLU00646.1;…;RLT92980.1;…`)
- `202` → sequence was replaced (A's `202` had a placeholder/empty accession)

**2 "new" sequences in B** are just the canonical re-sequenced versions of:

- `202` (now `YNPsite05_CeleraDRAFT_401410`; `7QJM`, _Chloroflexus sp._)
- the merged `607`

## Field enrichment on the 479 shared rows

B mostly adds cross-references and synonyms — nothing structural was rewritten.

| Change | Rows affected |
|--------|---------------|
| `accession` expanded (added UniProt / PDB / GenBank IDs) | 89 |
| `enzyme_name` expanded (added aliases) | 95 |
| `organism` expanded (added taxonomic synonyms) | 108 |
| `sequence_length` newly populated | 479 / 479 |
| `component_id` changed | 13 |
| `cath` changed | 1 |

## Bottom line

B ("tmp") is the **same PAZy dataset** as A, with:

1. a new `sequence_length` column,
2. richer `accession` / `enzyme_name` / `organism` annotations on ~20% of rows, and
3. 5 rows deduplicated down into 2 (484 → 481).

**No enzymes were genuinely dropped** — the "missing" ones were merged into existing entries.

## Open follow-ups

- Produce a row-by-row changelog CSV (every accession/name/organism token added per sequence).
- Inspect the 13 `component_id` changes and the 1 `cath` change specifically.
