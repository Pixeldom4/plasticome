# Migration record

Dated reorganisations of this repository, newest first. Paths inside each entry
describe the tree **as it was at that date**; they are deliberately not rewritten
by later moves, because the entry is a record of what happened rather than a live
index.

For the current layout see [README.md](README.md).

---

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
