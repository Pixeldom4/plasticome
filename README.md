# Plasticome

Curation and comparative-sequence pipeline for the plasticome, a catalogue of
plastic-degrading enzymes assembled from PAZy, the Erickson functional-metagenomics
set, and literature curation.

The pipeline merges the v1 and v2 tables into a union, backfills every sequence
from its accession, clusters at 90% amino-acid identity, partitions the centroids
into homology components by all-vs-all alignment, and emits the centroid FASTA. A
parallel branch collapses the union to one row per distinct sequence for accession
resolution. Flowchart:
[archive/figure-design/plasticome-v2-figures/pipeline.mmd](archive/figure-design/plasticome-v2-figures/pipeline.mmd).

---

## Layout

Organized **by run**, not by pipeline step. Each execution owns one directory
under `runs/`; everything shared across runs is hoisted out of the run tree. The
step vertical had drifted, with the same run named differently in different step
folders and handoff artifacts hand-copied between them. The run is the unit that
gets reasoned about, so it is the unit on disk.

```
plasticome/
├── pipeline_run.bash  end-to-end driver, steps 1-6 into one run directory
├── runs/              one directory per run          → runs/README.md
├── lib/               step-generic code, no data     → per-step READMEs below
├── bin/               third-party binaries (usearch11/12, DIAMOND)
├── source-data/       immutable upstream inputs of record
└── archive/           retired material, caches, figures, one-off investigations
```

| Directory | Step | Contents |
|---|---|---|
| `lib/01 union/` | 1 | `build_union.py`, `fetch_sequences.py` |
| `lib/02 clustering/` | 2 | `cluster_reference_seeded.py` + method doc |
| `lib/03 alignment/` | 3 | all-vs-all component assignment, `config.py`, `scripts/`, notebook |
| `lib/04 fasta/` | 4 | TSV→FASTA generators, accession-order utilities |
| `lib/05 annotate/` | 5 | `annotate_union.py`, the union at one row per sequence |
| `lib/06 nr/` | 6 | `build_nr.py`, `crosswalk.py` |
| `lib/common/` | — | `membership.py`, shared by steps 5 and 6 |
| `lib/summary/` | — | `summarize_run.py`, rebuilds `summary.md` from disk |
| `lib/00-archived-steps/` | — | retired: `accession/`, `annotation/`, `collapse/` |

Each `lib/NN */README.md` documents its own step. `source-data/` holds
`plasticome.v1.1/` (212 rows), `plasticome.v260701/` (473 rows) and
`plasticome-curated-seeds/` (13 seeds); these are snapshots, never edited in place.

---

## Running it

```bash
./pipeline_run.bash <run-name>                 # new run under runs/<today>.<run-name>/
./pipeline_run.bash --run-dir runs/<dir>       # resume or re-run an existing one
./pipeline_run.bash <run-name> --to 6          # include steps 5 and 6
./pipeline_run.bash <run-name> --dry-run       # print the commands, run nothing
```

Range is 1 to 6, default 1 to 4. Steps 5 and 6 are opt-in. **Resumable**: a step
whose deliverable exists is skipped unless `--force`. `--help` lists the rest
(`--id`, `--engine`, `--only`, `--from`, `--v1`, `--v2`, `--seeds`).

`--engine` selects the engine for **both** steps 2 and 3, `diamond` (default) or
`usearch`. Not interchangeable: diamond enforces member coverage where usearch
recruits on identity alone, and their edge sets differ. State which produced a run.

**Docker** is required for usearch always, and for diamond only when no native one
is on PATH (`brew install diamond`). `bin/usearch11`/`12` are Linux x86-64 ELF with
no macOS build. Containers mount the repo root at `/d`, so inputs and outputs must
both live under it or they are invisible inside the container.

---

## Steps and deliverables

| # | In | Out | Rows |
|---|---|---|---|
| 1 | v1.1, v260701 | `01-union.tsv` | 609 |
| 2 | 01 + seeds | `02-clusters.tsv` | 411 |
| 3 | 02 | `03-alignment.tsv` | 411 |
| 4 | 03 + 01 | `04-<run>.fasta` / `.tsv` | 411 |
| 5 | 03 + 01 | `05-union-with-components.tsv` | 609 |
| 6 | 01 only | `06-nr.tsv`, `06-nr.fasta` | 493 |

Three deliverable sets: **A** the union (609), **B** the 90% centroids (411, the
node set for the component partition), **C** the 100% non-redundant set (493, for
resolving alternative accessions).

**B and C are parallel branches off step 1, not a chain.** Building B on top of C
would shrink the database steps 2 and 3 search, and a smaller database gives
smaller e-values for the same score, so the fixed `evalue < 1e-5` post-filter would
become more permissive. See `lib/06 nr/README.md`.

Step 4's FASTA header is five positional fields:
`>identifier|accession|genbank_accessions|pazy_id|component`. Field 3 is currently
always empty; the pipe is kept so the field count is fixed.

`06-nr.fasta` uses the same five-field shape, so one `split("|")` parser reads
both, but the two are confusable from a header line alone: `06-nr.fasta` has 493
records against step 4's 411, and its field 5 is empty until `crosswalk.py` fills
it from branch B. See `lib/06 nr/README.md`.

---

## Conventions

**Quote every `lib/` path.** Directory names contain spaces. Unquoted, argparse
reports `unrecognized arguments` rather than a missing file, which misreads as a
flag problem.

**`NN-<name>.<ext>`** at the top of a run directory, one deliverable per step, with
scratch in a `<stem>.intermediates/` sidecar. Runs before 2026-07-30 use a space
instead of the hyphen and are left as they are. See `runs/README.md`.

**Runs are self-contained for data, shared for code.** Copying a data artifact into
a run is expected; copying a script into one is not.

**Identical filename ≠ identical content.** Several artifacts recur under one name
with different bytes across passes. Compare by hash before assuming.

**Provenance files are not re-pathed.** The `*stats*.json`, `*.log` and `REPORT_*`
inside a sidecar record absolute paths as they were at execution time. Later
renames deliberately leave them stale: they record what ran, not a live index.

**Timestamps.** File mtime equals birthtime here, so a copy resets it. Where a date
is embedded in a filename or frontmatter, that date wins.

---

## Known gaps

- `todo/plasticome.v211-manual-removal/` was deleted from the working tree but is
  still tracked; it survives in git history only.
- The step-3 driver still defaults its deliverable to `03 alignment.tsv` (space,
  no hyphen). Pass `OUT_TSV=` to get the hyphenated name; changing the default
  would rename outputs for every older run.
- The pre-2026-07-29 snapshot did **not** cover `analysis/outputs/*` — a nested
  `.gitignore` excluded it. Those files migrated intact and are tracked now, but
  they had no undo path during the move.
- `v9`, `v10` and `v11` have no written report; `runs/2026-07-20_v8/REPORT_v8_consolidated-copy.md`
  is a duplicate of the v8 report. **Which alignment result is authoritative is
  still undecided.**
- Pass 1 and pass 3 of the singleton run disagree (47 vs 48 components) because
  they consumed different centroid sets. Not investigated.
- `runs/2026-07-22_v1-v2-union/clustering-readme.md` documents a 611→411 clustering
  whose artifacts are not in the repo; superseded by the reference-seeded runs.
- Per-run `manifest.yaml` files (parent run, tool versions, input md5s) are **not
  yet written**. Until they exist, lineage lives only in this README.

---

History of reorganisations: [MIGRATIONS.md](MIGRATIONS.md). Run index: [runs/README.md](runs/README.md).
