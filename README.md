# Plasticome

Curation and comparative-sequence pipeline for the plasticome — a catalogue of
plastic-degrading enzymes assembled from PAZy, the Erickson functional-metagenomics
set, and literature curation.

The pipeline compiles a curated enzyme table, resolves and validates every accession,
merges the v1 and v2 tables, clusters the union at 90% amino-acid identity, partitions
the cluster centroids into homology components by all-vs-all alignment, and emits the
centroid FASTA. The canonical flowchart lives in
[archive/figure-design/plasticome-v2-figures/pipeline.mmd](archive/figure-design/plasticome-v2-figures/pipeline.mmd).

---

## Layout

The repository is organized **by run**, not by pipeline step. Each execution of the
pipeline — or of a contiguous stretch of it — owns one directory under `runs/`.
Everything shared across runs is hoisted out of the run tree.

```
plasticome/
├── run_pipeline.sh end-to-end driver — steps 1-4 into one run directory
├── runs/           one directory per pipeline run
├── lib/            step-generic code — no data, no outputs
├── bin/            third-party binaries (usearch11/12, DIAMOND)
├── source-data/    immutable upstream inputs of record
└── archive/        retired material — fetch caches, figures, one-off investigations
```

**Why run-first.** The step vertical had drifted: the same run carried different
names in different step folders, handoff artifacts were hand-copied between steps,
and cross-step paths had already broken in an earlier move. Reports were already
written per-run, never per-step. The run is the unit that gets reasoned about, so
it is now the unit on disk.

### `lib/` — step-generic code

Live pipeline steps are numbered in execution order; the number is part of the
directory name, so the pipeline reads top-to-bottom in a directory listing.

| Directory | Step | Contents |
|---|---|---|
| `lib/01 union/` | 1 | `build_union.py` (v1 ∪ v260701 merge), `fetch_sequences.py` (sequence backfill) |
| `lib/02 clustering/` | 2 | `cluster_reference_seeded.py` + the `reference-seed-clustering.md` method doc |
| `lib/03 alignment/` | 3 | all-vs-all component assignment — `config.py`, `run.sh`, `run_from_clusters.sh`, `scripts/`, `PETadex_alignment.ipynb` |
| `lib/04 fasta/` | 4 | TSV→FASTA generators and accession-order utilities — see `lib/04 fasta/README.md` |
| `lib/05 annotate/` | 5 | `annotate_union.py` — the union at one row per sequence, with its cluster and component |
| `lib/06 nr/` | 6 | `build_nr.py` (100% non-redundant set) + `crosswalk.py` (the B × C crosswalk) |
| `lib/common/` | — | `membership.py`, the positional cluster-membership reconstruction steps 5 and 6 share |
| `lib/summary/` | — | `summarize_run.py`, which rebuilds `summary.md` from whatever is on disk |
| `lib/00-archived-steps/` | — | retired steps: `accession/` (8-tier validation), `annotation/` (paper fetching), `collapse/` (substring/numbered-enzyme scans) |

Directory names contain spaces, so **every path must be quoted** in the shell. An
unquoted path splits into two argv entries and argparse rejects it with an
`unrecognized arguments` error rather than a missing-file error, which reads as a
flag problem and is easy to misdiagnose.

`lib/03 alignment/config.py` is run-agnostic: it resolves the repo root and reads
`PLASTICOME_RUN` to select which run's data to operate on.

```bash
PLASTICOME_RUN=runs/2026-07-23_v11 python3 "lib/03 alignment/scripts/step23_graph.py" ...
```

### `source-data/` — inputs of record

`plasticome.v1.1/` (the 212-row v1 table), `plasticome.v260701/` (the 473-row PAZy
pull), and `plasticome-curated-seeds/` (the 13 curated seeds). These are upstream
snapshots — never regenerated, never edited in place.

### `bin/` and `archive/`

**DIAMOND** is the pipeline's engine for both the clustering (step 2) and the
all-vs-all alignment (step 3). A native `diamond` on PATH (`brew install diamond`)
is used when there is one, and needs no Docker; `bin/diamond` is the bundled Linux
x86-64 ELF fallback, run through Docker `linux/amd64` (`bin/diamond.sh` wraps that
call for interactive use).

`usearch11` / `usearch12` are Linux x86-64 ELF binaries with no macOS distribution,
so they **only run through Docker `linux/amd64`** — Docker Desktop must be running
before any step that uses them: `--engine usearch`, `lib/03 alignment/run.sh`, the
notebook, and the Step-0 sanity fixture.

`archive/cache/` holds the NCBI / BLAST / supplementary-file caches and the
~26-paper corpus; all of it is keyed on stable identifiers, so deleting it only
costs re-fetch time. `archive/figure-design/` holds the pipeline figures and
sequence logos; `archive/usearch12-bugreport/` is the usearch v11-vs-v12
investigation, which is not a pipeline run.

---

## Running the pipeline

Use the driver — it does all four steps, in order, into one run directory:

```bash
./run_pipeline.sh <run-name>                    # new run under runs/<today>.<run-name>/
./run_pipeline.sh --run-dir runs/<dir>          # resume or re-run an existing one
./run_pipeline.sh <run-name> --from 3           # only steps 3-4
./run_pipeline.sh <run-name> --dry-run          # print the commands, run nothing
```

It is **resumable**: a step whose deliverable already exists is skipped unless
`--force`, so an interrupted run picks up where it stopped and re-running only the
last step is cheap. It also adopts a deliverable that a previous run named
differently, rather than building a second one beside it. `--help` lists the rest
(`--id`, `--engine`, `--v1`, `--v2`, `--seeds`, `--only`, `--to`).

`--engine` selects the engine for **both** steps 2 and 3 — `diamond` (default) or
`usearch`, the engine used through 2026-08-04. They are not interchangeable outputs:
diamond enforces member coverage where usearch recruited on identity alone, and the
two aligners' edge sets differ slightly, so state which one produced a run.

The driver passes `V1=` explicitly to step 3, which sidesteps the stale-path bug
described under *Known gaps* below. Invoking the steps by hand does not, and loses
the v1 overlay silently.

### The steps by hand

Equivalent to the above, if you need to vary something the driver does not expose.
Each step writes exactly one deliverable at the top of the run directory;
everything else lands in that deliverable's sidecar.

```bash
RUN=runs/2026-07-30.plasticome.v1.212.union-spec

# 1. union — 212 + 473 rows in, 607 out, then backfill sequences from accession
python3 "lib/01 union/build_union.py" \
    source-data/plasticome.v1.1/plasticome.v1.1.tsv \
    source-data/plasticome.v260701/cleaned_pazy-260701_retrieving_from_accession.tsv \
    -o "$RUN/01-plasticome_v1.v260701-union.tsv"
python3 "lib/01 union/fetch_sequences.py" "$RUN/01-plasticome_v1.v260701-union.tsv"

# 2. reference-seeded clustering at 90% — 607 rows → clusters  [diamond]
python3 "lib/02 clustering/cluster_reference_seeded.py" \
    "$RUN/01-plasticome_v1.v260701-union.tsv" \
    source-data/plasticome-curated-seeds/plasticome-curated-seeds.tsv \
    -o "$RUN/02-clusters.tsv"

# 3. all-vs-all alignment of the centroids → adds component_id  [diamond]
OUT_TSV="$RUN/03-alignment.tsv" "lib/03 alignment/run_from_clusters.sh" "$RUN/02-clusters.tsv"

# 4. centroid FASTA
python3 "lib/04 fasta/clusters_to_fasta.py" "$RUN/03-alignment.tsv" \
    -o "$RUN/04-plasticome.v1.212.union-spec.fasta"
```

Steps 2 and 3 shell out to DIAMOND, which runs natively when one is on PATH and
otherwise falls back to `bin/diamond` under Docker — as does all of `--engine
usearch`, which fails immediately if Docker Desktop is not running. Step 1 hits
NCBI/UniProt; export `NCBI_API_KEY` first to lift the rate limit from 3 to 10
requests/second.

**Row counts through the pipeline**, as a smoke test — if these drift, something
upstream changed:

| Step | Out | Note |
|---|---|---|
| 1 union | 607 rows | 395 `v260701`, 78 `both`, 134 `v1.1`; 607/607 sequences populated |
| 2 clustering | 411 clusters | 13 seed clusters (ids 1–13) + 398 de novo |
| 3 alignment | 411 rows | + `component_id`; 46 components, largest 232, 25 singletons |
| 4 FASTA | 411 records | 5-field pipe header |

Step 2's internal funnel is worth knowing when a count looks wrong: the 13 curated
seeds enter as centroids, **22** union rows are identical to a seed and fold
straight into its cluster (marked `dup` rather than a percent identity), leaving
**585**; closed-reference search recruits **28** of those to seeds (4.8%); the
remaining **557** cluster de novo into **398**. So 13 + 398 = 411 clusters, and
the member column totals 620 = 607 union rows + the 13 seed centroids.

### FASTA header format

```
>identifier|accession|genbank_accessions|pazy_id|component
>PL1|WP_054022242.1||1|C001
```

`identifier` is a clean 1..N re-index over the run's own rows, **not** carried
across runs — `PL16` in a centroid FASTA is a different enzyme than `PL16` in the
607-row union table. `genbank_accessions` is reserved and currently emitted empty;
the pipe is kept so the field count is fixed at five. Records are ordered by
pazy_id numerically, then component, then sequence, with rows lacking a pazy_id
sorted last — so the file is **not** grouped by component.

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

---

## Conventions

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

**Runs are self-contained for data, shared for code.** Copying a data artifact
into a run is correct and expected; copying a script into a run is not.

**Provenance files are not re-pathed.** The `*stats*.json`, `*.log` and `REPORT_*`
files inside a sidecar record absolute paths as they were at execution time. When a
directory is later renamed, those strings are deliberately left stale — they are a
record of what ran, not a live index. Their md5s and statistics remain valid.

**Timestamps.** File mtime equals birthtime throughout this repo, so a copy resets
it. Where a date is embedded in a filename or report frontmatter, that date wins.
Three run folders were re-stamped by a bulk move on 2026-07-13 11:18 and carry
their true dates only in their filenames: `step0` (07-02) and `v2` (07-09).

**Identical filename ≠ identical content.** Several artifacts recur under one name
with different bytes across passes. Always compare by hash before assuming two
copies are the same file.

**Docker.** Containerized invocations mount the **repo root** at `/d`, so tool paths
inside the container are `/d/bin/usearch11` and data paths are `/d/$REL/...`. The
clustering step is the one exception: it mounts `bin/` at `/b` and only the run's
intermediates directory at `/w`, since that is all it reads or writes. Relocated
drivers discover the repo root by walking up to the directory containing both `lib/`
and `bin/`, so they work from any depth. A consequence: **inputs and outputs must
both live under the repo root**, or they are invisible inside the container. None of
this applies when DIAMOND runs natively, which is the normal case.

---

## Migration record — 2026-07-29

Reorganized from a step vertical (10 top-level step folders) to the run vertical
above. A full `git` snapshot was taken first (commit `e0c9cfc`); every distinct
blob from that snapshot is still present, verified by hash comparison.

**Moves.** `petadex-alignment/`, `petadex-clustering/`, `petadex-collapsed/`,
`fasta-generation/`, `accession-validation/`, `activity-annotation/`,
`plasticome-v1-cleaning/`, `v1-v2-union/` were dissolved into `runs/`, `lib/`,
`bin/`, `sources/` and `cache/`.

**Path repairs (13 files).** Two hardcoded absolute paths in `run_clustering.sh`,
six relative tool references in the alignment drivers, two in the usearch12 bug
report, plus `config.py`, `run.sh` and `collapse/config.yaml`.

**Resolved during migration.**

- `plasticome.union.v1.v2` is the **early run** of `plasticome.v1.v260701-union`, not a separate dataset.
- The `singleton-cleaned-union` folders named `.v1.v2` and `.v1.v260701` were one run under two names; canonical name is `.v1.v260701`.
- `usearch/v9/` (input) and `analysis/outputs/v9/` (output) were one run.
- The curated-seed directories (`curated_seed`, `curated_seed_putative`, `curated-seed-validation`) were one run plus duplicates and an empty stub.
- `petadex-collapsed` was one run, and it spans v9 and v10.
- `accession-validation/results/` and `out/` were two different runs (07-10 and 07-21).
- `plasticome.v1.1` is a source datafile, not a run.

**Deduplicated.** Two byte-identical copies of `usearch11` and of
`plasticome.v1.csv`. Two files sharing the basename `tsv_to_fasta.py` were *not*
duplicates — the PAZy/PL-header-specific one was preserved as
`tsv_to_fasta_pl_headers.py`, then retired on 2026-07-30 once both input schemas
were settled (see `lib/04 fasta/README.md`).

## Migration record — 2026-08-17

| Was | Is |
|---|---|
| `lib/fasta/` | `lib/04 fasta/` |

Every live step directory is now numbered in execution order, so `lib/` reads
top-to-bottom as the pipeline runs. `lib/fasta/` was the only live step still
named by what it does rather than when it happens, which put step 4 between
`03 alignment` and `05 annotate` alphabetically but not visibly.

The new directory name contains a space, like the others, so the examples in
`lib/04 fasta/README.md` are now quoted. `fix_accession_order.py` still resolves
its input as a bare filename and so still has to be run from inside that
directory.

Added `lib/05 annotate/` (step 5, which previously had no code at all -- the
`05-union-with-components.tsv` on `2026-08-06.final-usearch.2` was built outside
the repo), `lib/06 nr/` (step 6), and `lib/common/membership.py`, shared by both.
The driver's step range is now 1-6; the default is still 1-4, so a plain run
produces exactly what it always did.

Dated migration prose below is left as written. Those entries describe paths as
they were at the time, and rewriting them would make the record inaccurate rather
than current.

## Migration record — 2026-07-30

A second reorganization, on top of the run vertical. **Not yet committed**; the
whole of it is uncommitted working-tree state.

| Was | Now |
|---|---|
| `sources/` | `source-data/` |
| `cache/` | `archive/cache/` |
| `figure-design/` | `archive/figure-design/` |
| `experiments/usearch12-bugreport/` | `archive/usearch12-bugreport/` |
| `lib/union/` | `lib/01 union/` |
| `lib/clustering/` | `lib/02 clustering/` |
| `lib/alignment/` | `lib/03 alignment/` |
| `lib/accession/`, `lib/annotation/`, `lib/collapse/` | `lib/00-archived-steps/` |
| `todo/` | *(removed from the working tree)* |

Added `run_pipeline.sh` (the end-to-end driver), `lib/fasta/clusters_to_fasta.py`
(step 4), and the two `2026-07-30.*` runs.

**Path repairs.** The move left five live references pointing at the old
`sources/` tree. All are fixed: `lib/03 alignment/config.py` (`SOURCES`, and `LIB`,
which still said `lib/alignment`), `run_from_clusters.sh` (the `V1` default), and
the usage examples in `lib/01 union/README.md`, `lib/fasta/README.md` and
`build_union.py`. Every path in `config.py` now resolves.

`source-data/plasticome.v1.csv` and `source-data/pazy_pull_2026-06-30.csv` were
restored from commit `8e243d1`, which is where they were last present — the move
had dropped them entirely. They are byte-identical to the committed blobs.
`config.V1_CSV` points at the first, and it loads 205 overlay entries, so the
`--v1` defaults of `step0_sanity.py` and `step1_nodes.py` work again.

The v1 overlay is **descriptive only** — it fills `v1_component` and `cath`, and
the partition never reads them. That is why its absence produced a warning rather
than an error, and why it could go unnoticed. Two different tables are in use as
the overlay: `config.V1_CSV` is `plasticome.v1.csv` (213 rows, used by the frozen
v4–v8 PAZy lineage) while `run_from_clusters.sh` defaults to `plasticome.v1.1.tsv`
(212 rows). Both were left as they were — unifying them changes what the frozen
runs produce and is a scientific call, not a path repair.

### Known gaps

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
