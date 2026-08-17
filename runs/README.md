# Runs

One directory per pipeline run, or per contiguous stretch of one. Everything
shared across runs lives outside this tree: code in `lib/`, binaries in `bin/`,
inputs in `source-data/`.

Data files here are gitignored; `*.md`, `*.py`, `*.sh` and `*stats*.json` are not,
so the human-readable record and the provenance survive in git.

---

## Layout of a run

**Adding a run.** Create `runs/<YYYY-MM-DD>.<name>/` and a `readme.md` recording
what changed versus the parent run. Reference `lib/`, `bin/`, `source-data/` and
`archive/cache/` — never copy code into a run.

**Run-internal layout.** Adopted 2026-07-30, replacing the earlier
`inputs/` + `work/` + `results/` split:

```
runs/<run>/
├── 01-<name>.tsv                deliverable of step 1
├── 02-<name>.tsv                deliverable of step 2
├── 02-<name>.intermediates/       └ that step's scratch
├── 03-<name>.tsv
├── 03-<name>.intermediates/
├── 04-<name>.fasta
└── readme.md
```

- **`NN-<name>.<ext>`** — two-digit step prefix, hyphen, name. Spaces in filenames
  are out; they forced quoting on every invocation and split silently when missed.
  Runs before 07-30 use `NN <name>.<ext>` with a space and are left as they are.
- **One deliverable per step.** A step writes exactly one file at the top level.
- **`<stem>.intermediates/` sidecars** — the extension is *replaced*, not appended,
  so the folder reads as a peer of its deliverable rather than a second extension
  on it. The older form was `<deliverable>.work/`; folders written before 07-30
  keep it. Sidecars are regenerable scratch that nothing downstream reads — the
  only parts worth keeping are the `*stats*.json`, which pin the tool version and
  thresholds that produced the result.

---

## Run index

Two lineages. `step0`–`v11` align the **PAZy-only** `cleaned_pazy_final` table;
the 07-23 and later runs align **union-derived** centroids. Chronological order is
*not* causal order across the two — read `parent` before assuming lineage.

| Run | Date | Lineage | What it did |
|---|---|---|---|
| `2026-07-02_step0` | 07-02 | PAZy | engine sanity check — 213 v1 nodes → 42 components |
| `2026-07-09_activity-annotation` | 07-09 | — | DOI list + paper corpus for activity/proof QC |
| `2026-07-09_v2` | 07-09 | PAZy | first component assignment + sensitivity models |
| `2026-07-10_accession-validation` | 07-10 | PAZy | first validation pass (tiers 1/2/S) |
| `2026-07-13_v3` | 07-13 | PAZy | re-run on combined node set |
| `2026-07-13_v4` | 07-13 | PAZy | paper-faithful method fixed (both orientations, max-bits HSP) |
| `2026-07-16_v5` | 07-16 | PAZy | `cleaned_pazy_final.csv`, 489 nodes |
| `2026-07-17_v6` | 07-17 | PAZy | TSV input, stable `PL<n>` node labels |
| `2026-07-17_v7` | 07-17 | PAZy | 490 nodes; flagged short-node/re-sequencing issues |
| `2026-07-20_v8` | 07-20 | PAZy | curated 484-node set — the frozen reference result (46 components) |
| `2026-07-20_collapse-v9-v10` | 07-20→07-21 | PAZy | **the full 2.4→2.6 loop, one run** (see below) |
| `2026-07-22_v1-v2-union` | 07-22 | union | v1 cleaning + v1∪v2 merge → 611 records |
| `2026-07-23_v11` | 07-23 | PAZy | alignment of `plasticome.v1.1` centroids |
| `2026-07-23_curated-seed-validation` | 07-23 | union | 14 putative seeds in → **13 curated seeds out** |
| `2026-07-23_plasticome.union.v1.v2` | 07-23 | union | early run of the union centroid alignment |
| `2026-07-23_singleton-cleaned-union.v1.v260701` | 07-23 | union | **three clustering passes** (see below) |
| `2026-07-23_plasticome.v1.v260701-union` | 07-23 | union | later run of the union centroid alignment |
| `2026-07-30.plasticome.v1.212` | 07-30 | union | pre-spec union rebuild — **superseded**, and its `02 clusters.tsv` is stale |
| `2026-07-30.plasticome.v1.212.union-spec` | 07-30 | union | **current** — union rebuilt to spec, clustered, aligned, FASTA emitted |

The 07-30 runs use `.` rather than `_` after the date. Older runs keep `_`; both
forms are in the tree and neither is parsed by anything.

### Runs that span more than one step folder

**`2026-07-20_collapse-v9-v10`** is one continuous session (07-20 16:14 → 07-21 18:33)
that executed the flowchart's whole 2.4 → 2.5 → 2.6 loop. Timestamps establish the
interleaving:

```
accession-validation/out    →  collapse/  →  v9 align  →  singletons-validated
   →  v10 re-align  →  final-diff  →  identifier_assignment  →  DIAMOND arm
```

It is one run directory with `accession-validation/`, `collapse/`,
`alignment-usearch/{v9,v10}/`, and `alignment-diamond/` inside it.

**`2026-07-23_curated-seed-validation`** exists to test seed non-redundancy. Its
14-seed input is in `inputs/`; the all-vs-all put **BhrPETase (GBD22443.1, PL17)
in component C002 alongside ISPETase (WP_054022242.1, PL1)**, so BhrPETase was
dropped. The resulting **13-seed set is this run's output** and the input to every
downstream clustering run; it now lives at
`source-data/plasticome-curated-seeds/plasticome-curated-seeds.tsv`.

**`2026-07-23_singleton-cleaned-union.v1.v260701`** clustered **three times**,
producing three genuinely different centroid sets — this is preserved, not deduped:

| Pass | Centroids md5 | Aligned at | Components |
|---|---|---|---|
| 1 | `88816d81…` (15:17) | `pass1-alignment/` | 47 |
| 2 | `cc2b3269…` (17:45) | *(no alignment recorded)* | — |
| 3 | `f4367699…` (17:54) | `pass3-alignment/` | 48 |

`pass3-alignment/` also holds the usearch-vs-DIAMOND comparison
(`REPORT_compare_2026-07-23.md`, `node_disagreements.tsv`).
