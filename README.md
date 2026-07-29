# Plasticome

Curation and comparative-sequence pipeline for the plasticome — a catalogue of
plastic-degrading enzymes assembled from PAZy, the Erickson functional-metagenomics
set, and literature curation.

The pipeline compiles a curated enzyme table, resolves and validates every accession,
merges the v1 and v2 tables, clusters the union at 90% amino-acid identity, and
partitions the cluster centroids into homology components by all-vs-all alignment.
The canonical flowchart lives in
[figure-design/plasticome-v2-figures/pipeline.mmd](figure-design/plasticome-v2-figures/pipeline.mmd).

---

## Layout

The repository is organized **by run**, not by pipeline step. Each execution of the
pipeline — or of a contiguous stretch of it — owns one directory under `runs/`.
Everything shared across runs is hoisted out of the run tree.

```
plasticome/
├── runs/           one directory per pipeline run, named <date>_<run-name>
├── lib/            step-generic code — no data, no outputs
├── bin/            third-party binaries (usearch11/12, DIAMOND)
├── sources/        immutable upstream inputs of record
├── cache/          run-independent fetch caches and the paper corpus
├── experiments/    investigations that are not pipeline runs
├── figure-design/  pipeline figures and sequence logos (not part of a run)
└── todo/           pending work, not yet run
```

**Why run-first.** The step vertical had drifted: the same run carried different
names in different step folders, handoff artifacts were hand-copied between steps,
and cross-step paths had already broken in an earlier move. Reports were already
written per-run, never per-step. The run is the unit that gets reasoned about, so
it is now the unit on disk.

### `lib/` — step-generic code

| Directory | Contents |
|---|---|
| `lib/accession/` | 8-tier accession validation (`tier*.py`, `ncbi.py`, `blast.py`, `supp.py`, `run_all.py`) |
| `lib/alignment/` | consolidated PETadex component pipeline — `config.py`, `run.sh`, `scripts/`, `PETadex_alignment.ipynb` |
| `lib/annotation/` | paper fetching for activity annotation |
| `lib/clustering/` | generalized reference-seeded clustering method |
| `lib/collapse/` | substring / numbered-enzyme collapse scans |
| `lib/fasta/` | TSV→FASTA generators and accession-order utilities |
| `lib/union/` | v1∪v2 merge and sequence backfill |

`lib/alignment/config.py` is run-agnostic: it resolves the repo root and reads
`PLASTICOME_RUN` to select which run's data to operate on.

```bash
PLASTICOME_RUN=runs/2026-07-23_v11 python3 lib/alignment/scripts/step23_graph.py ...
```

### `sources/` — inputs of record

`plasticome.v1.1/` (the v1 table and its FASTA), `plasticome.v1.csv`, and
`pazy_pull_2026-06-30.csv`. These are upstream snapshots — never regenerated,
never edited in place.

### `bin/` and `cache/`

`usearch11` / `usearch12` are Linux x86-64 ELF binaries and **only run through
Docker `linux/amd64`** on this arm64 macOS host. DIAMOND likewise, via
`bin/diamond.sh`. `cache/` holds the NCBI / BLAST / supplementary-file caches and
the ~26-paper corpus; all of it is keyed on stable identifiers and shared by every
run, so deleting it only costs re-fetch time.

---

## Run index

Two lineages. `step0`–`v11` align the **PAZy-only** `cleaned_pazy_final` table;
the 07-23 runs align **union-derived** centroids. Chronological order is *not*
causal order across the two — read `parent` before assuming lineage.

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
dropped. The resulting **13-seed set in `results/` is this run's output** and the
input to every downstream clustering run.

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

**Adding a run.** Create `runs/<YYYY-MM-DD>_<name>/`, put the run's own inputs in
`inputs/`, intermediates in `work/`, deliverables in `results/`, and a `README.md`
recording what changed versus the parent run. Reference `lib/`, `bin/`, `sources/`
and `cache/` — never copy code into a run.

**Runs are self-contained for data, shared for code.** Copying a data artifact
into a run is correct and expected; copying a script into a run is not.

**Timestamps.** File mtime equals birthtime throughout this repo, so a copy resets
it. Where a date is embedded in a filename or report frontmatter, that date wins.
Three run folders were re-stamped by a bulk move on 2026-07-13 11:18 and carry
their true dates only in their filenames: `step0` (07-02), `v2` (07-09), and
`sources/pazy_pull_2026-06-30.csv` (06-30).

**Identical filename ≠ identical content.** Several artifacts recur under one name
with different bytes across passes. Always compare by hash before assuming two
copies are the same file.

**Docker.** Every usearch/DIAMOND invocation mounts the **repo root** at `/d`, so
tool paths inside the container are `/d/bin/usearch11` and data paths are
`/d/$REL/...`. Relocated drivers discover the repo root by walking up to the
directory containing both `lib/` and `bin/`, so they work from any depth.

---

## Migration record — 2026-07-29

Reorganized from a step vertical (10 top-level step folders) to the run vertical
above. A full `git` snapshot was taken first (commit `e0c9cfc`); every distinct
blob from that snapshot is still present, verified by hash comparison.

**Moves.** `petadex-alignment/`, `petadex-clustering/`, `petadex-collapsed/`,
`fasta-generation/`, `accession-validation/`, `activity-annotation/`,
`plasticome-v1-cleaning/`, `v1-v2-union/` were dissolved into `runs/`, `lib/`,
`bin/`, `sources/` and `cache/`. `figure-design/` and `todo/` were left untouched.

**Path repairs (13 files).** Two hardcoded absolute paths in `run_clustering.sh`,
six relative tool references in the alignment drivers, two in the usearch12 bug
report, plus `lib/alignment/config.py`, `lib/alignment/run.sh` and
`lib/collapse/config.yaml`. All 15 shell drivers pass `bash -n`; `config.py`
resolves every path.

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
duplicates — the PAZy/PL-header-specific one is preserved as
`lib/fasta/tsv_to_fasta_pl_headers.py`.

**Deleted.** Three empty directories; `__pycache__`, `.DS_Store`, and a `.gitkeep`.

### Known gaps

- The pre-migration snapshot did **not** cover `analysis/outputs/*` — a nested `.gitignore` excluded it. Those files migrated intact and are tracked now, but they had no undo path during the move.
- `v9`, `v10` and `v11` have no written report; `runs/2026-07-20_v8/REPORT_v8_consolidated-copy.md` is a duplicate of the v8 report carried by the old consolidated `analysis/` folder. **Which alignment result is authoritative is still undecided.**
- Pass 1 and pass 3 of the singleton run disagree (47 vs 48 components) because they consumed different centroid sets. Not investigated.
- `runs/2026-07-22_v1-v2-union/clustering-readme.md` documents a 611→411 clustering whose artifacts are not in the repo; superseded by the reference-seeded runs.
- Per-run `manifest.yaml` files (parent run, tool versions, input md5s) are **not yet written**. Until they exist, lineage lives only in this README.
