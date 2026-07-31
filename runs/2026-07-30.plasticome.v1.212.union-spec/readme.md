# Union rebuilt to the documented spec (2026-07-30)

Regenerated from:

1. `sources/plasticome.v1.1/plasticome.v1.1.tsv` (212 rows)
2. `sources/plasticome.v260701/cleaned_pazy-260701_retrieving_from_accession.tsv` (473 rows)

via `lib/union/build_union.py` (Stage 1) + `lib/union/fetch_sequences.py` (Stage 2).

Supersedes `runs/2026-07-30.plasticome.v1.212/`, which was built by the
pre-2026-07-30 script and diverged from the spec in four ways: v1.1 sequences
were carried forward wholesale (making Stage 2 a no-op), manual Erickson
assignments were overwritten on merge, `plasticome_id` and the pazy_id/component
ordering were absent, and the `source` labels were `v1`/`v2`/`merged`. That run
directory is left untouched; its `02 clusters.tsv` derives from the old union
and is **stale** with respect to this one.

## Result

607 rows — 395 `v260701`, 78 `both`, 134 `v1.1`. 607/607 sequences populated,
0 unresolved. Stage 1 verified byte-identical across two runs; Stage 2 re-run
fills 0 rows.

---

# Layout & naming convention

```
01-plasticome_v1.v260701-union.tsv          607 rows, the union
02-clusters.tsv                             411 rows, one per 90% cluster
02-clusters.intermediates/                    └ USEARCH seed/closed-ref/de-novo scratch
03-alignment.tsv                            411 rows, + component_id
03-alignment.intermediates/                   └ all-pairs alignment + graph scratch
04-plasticome.v1.212.union-spec.fasta       411 records, centroid FASTA
fetched_sequences.json                      Stage 2 per-accession provenance
```

**Preferences this run establishes** (applied 2026-07-30, adopted for future runs):

1. **Numbered, hyphenated step files.** `NN-<name>.<ext>` — a two-digit step
   prefix, then a hyphen, then the name. Spaces in filenames are out: every
   invocation needed quoting, and an unquoted path silently split into two
   argv entries rather than failing loudly.
2. **Sidecars are `<stem>.intermediates/`, not `<file>.work/`.** The extension is
   **replaced**, not appended, so the folder reads as a peer of the deliverable
   (`02-clusters.tsv` ↔ `02-clusters.intermediates/`) instead of as a second
   extension on it. `.work` was opaque about what it held; `.intermediates` says
   it. Nothing downstream reads these folders — they are regenerable scratch, and
   the only parts worth keeping are the `*_stats*.json`, which pin the tool
   version and thresholds that produced `component_id`.
3. **One deliverable per step, sidecar for everything else.** A step writes
   exactly one file at the top level; all scratch goes in its sidecar.

The generator defaults were updated to match, so a re-run reproduces these names
rather than recreating `.work/`: `lib/02 clustering/cluster_reference_seeded.py`
(`--workdir` default), `lib/03 alignment/run_from_clusters.sh` (`OUT` default),
and `lib/03 alignment/README.md`. Sidecars written before 2026-07-30 in other run
directories still use the old `.work/` form and were left untouched.

**Two caveats on the rename.**

* Four provenance files inside `03-alignment.intermediates/` (`combined_stats.json`,
  `stats_clusters_2026-07-30.json`, `comparison_vs_2026-07-23.json`,
  `REPORT_compare_vs_2026-07-23.md`) plus `usearch.log` record the absolute path
  as it was at execution time and so still say `03 alignment.tsv.work`. They were
  **deliberately not rewritten** — they are a record of what ran, not a live
  index. Their `input_md5` and stats remain correct; only the directory string is
  historical.
* The step-file names here are ahead of the generators, which still default to
  writing `03 alignment.tsv` (space, no hyphen). Renaming those defaults would
  change output names for every other run directory, so it was left alone; pass
  `OUT_TSV=` / `-o` explicitly when re-running this one.

## Step 04 — centroid FASTA

`04-plasticome.v1.212.union-spec.fasta`, built by `lib/fasta/clusters_to_fasta.py`
from `03-alignment.tsv`. 411 records, one line per sequence, five-field header:

```
>identifier|accession|genbank_accessions|pazy_id|component
>PL1|WP_054022242.1||1|C001
```

`identifier` is a clean 1..411 re-index (`PL<n>`), **not** carried over from the
union table's `plasticome_id` — the row set is the 411 centroids, not the 607
union rows, so `PL16` here is a different enzyme than `PL16` in
`runs/2026-07-23_plasticome.union.v1.v2/`. `genbank_accessions` is reserved and
emitted empty; the pipe is kept so the field count is fixed at five.

Ordering mirrors `build_union.py:sort_key` — pazy_id numerically, then component
numerically, then amino-acid sequence — with rows lacking a pazy_id sorting last
(ids 391–411). Since pazy_id is unique across the 390 rows that have one, the
component and sequence keys only discriminate within that blank tail; the file is
**not** grouped by component. 25 records have a blank accession
(`>PL391||||C001`), faithful to the source.

---

# Diff vs `runs/2026-07-30.plasticome.v1.212/01 plasticome.v1-v2-union.tsv`

A line-by-line diff is uninformative — the row order changed completely (1 of
607 rows sits at the same file position). Rows below were matched on identity
(`enzyme_name`, `accession`, `pazy_id`, falling back to `enzyme_name` + `pazy_id`).

## Structural

| | old | new |
|---|---|---|
| Columns | `enzyme_name, accession, pazy_id, aa_sequence, source` | `plasticome_id, enzyme_name, accession, pazy_id, aa_sequence, source` |
| Rows | 607 | 607 |
| Source labels | `v2` / `merged` / `v1` | `v260701` / `both` / `v1.1` |
| Order | v260701 rows in input order, v1.1 appended | `plasticome_id` 1–607 on pazy_id then component |

**Row membership is identical** — all 607 matched, 0 old-only, 0 new-only. After
remapping the labels, **0 rows changed source** (395 / 78 / 134 in both). The
merge decided exactly the same thing in both runs; only sequences, ordering, and
labels moved.

## Content — 5 cell-level differences across 607 rows

**`aa_sequence` — 4 rows**, all `both`, all the manual-assignment exception
firing. (`enzyme_name` in the union comes from v260701, which names these rows
`503`/`611`/`701`/`702`; they are `Enzyme 503 like` etc. in v1.1.)

| id | name | v1.1 name | accession | old len | new len | delta |
|---|---|---|---|---|---|---|
| 129 | 503 | Enzyme 503 like | EGD44994.1 | 294 | 267 | −27 |
| 133 | 611 | Enzyme 611 | WP_093412886.1 | 293 | 269 | −24 |
| 134 | 701 | Enzyme 701 | WP_104613137.1 | 301 | 270 | −31 |
| 135 | 702 | Enzyme 702 | ADM47605.1 | 262 | 271 | +9 |

Old took v260701's database record; new keeps v1.1's manually assigned Erickson
supplementary table D1 sequence. The deltas are the signal sequences / His-tags
the v1 cleaning step deliberately normalized away.

**`accession` — 1 row**: id 495, `jmPE13`, `'jmPE13'` → `''`. The placeholder is
no longer emitted as an accession (blank-accession count 39 → 40).

**whitespace — 1 row**: `UMG-SP-3` / WBR49958.1 had a trailing space in old,
trimmed in new. Does not show up as a sequence difference above because the
comparison strips both sides; the residue string is unchanged.

## The notable non-difference

The 116 v1.1-only rows that Stage 2 now retrieves from accession came back
**byte-identical** to the values the old run had carried forward from v1.1 —
they do not appear in the diff at all.

That is the substantive finding: the old blanket carry-forward was not
corrupting sequences, it was merely unverified. What the rewrite buys is
provenance, not corrections — except for the 4 Erickson rows above, where the
spec's manual-assignment exception was genuinely being violated.

Stage 2 cascade for those 116: NCBI 65 → UniProt +45 → RCSB/PDB +2 →
UniParc +4. Per-accession provenance in `fetched_sequences.json`.

## Open questions (not resolved by this run)

- The spec's Step 3 says "15 Erickson primary/only enzymes" but its breakdown
  sums to 20 (11 + 6 + 3), and v1.1 actually has **21** rows named
  `Enzyme <n>`. This run treats all 21 plus `jmPE13` (22 rows) as manual. The
  "removed 12 GenBank accessions" figure does check out exactly.
- The spec says v1 n=213 → 611 rows (398 / 75 / 138); the actual inputs give
  212 → 607 (395 / 78 / 134). One fewer input row does not explain merged
  going 75 → 78, so `plasticome.v1.1.tsv` has changed in more than row count
  since the spec was written.
