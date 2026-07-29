# PETadex-alignment — consolidated pipeline

Self-contained, single-folder version of the PETadex component-assignment analysis
(the v8 run). One central notebook orchestrates the whole pipeline from the frozen
input and then analyses the result; the same steps are available as standalone
scripts and a `run.sh` driver.

Reproduces the v8 result **exactly**: 484 md5-unique nodes → **46 components**,
33,101 passing edges, largest component 283 (the α/β-hydrolase superfamily). The
partition is byte-for-byte identical to `../v8/` (same clustering, same C-labels on
all 484 nodes).

## Layout

```
analysis/
├── PETadex_alignment.ipynb   ← central notebook: orchestrate (Steps 0–3) + analyse
├── run.sh                    ← shell twin of the notebook (same steps/outputs)
├── config.py                 ← all paths + scientific constants + usearch(Docker) helper
├── requirements.txt          ← notebook deps (scripts are pure stdlib)
├── scripts/
│   ├── step0_sanity.py       ← rebuild 213 v1 nodes (engine gate → 42 comps / 3178 edges)
│   ├── step1_nodes.py        ← curated PAZy TSV → 484 md5-unique nodes + FASTA
│   ├── step23_graph.py       ← usearch pairs → filtered graph → connected components
│   └── annotate_source.py    ← join components back onto every source row
├── data/
│   ├── cleaned_pazy_final.tsv (frozen input, md5 1d13e83691cc767c7ba9635c9bd2ed60)
│   └── plasticome_v1.csv      (v1 overlay — descriptive only)
├── bin/usearch11             ← Linux ELF aligner (run only via Docker linux/amd64)
└── outputs/                  ← all generated artifacts (gitignored)
```

## Prerequisites

- **Docker** running (`docker ps` must work) — usearch is a Linux ELF binary run
  through `debian:bookworm-slim` on `linux/amd64`. There is no macOS/arm64 native path.
- The **`plasticome`** conda env for the notebook (jupyter, pandas, matplotlib, networkx —
  see `requirements.txt`). A kernel named *Python 3 (plasticome)* is registered.

## Run it

**Notebook (recommended):** open `PETadex_alignment.ipynb` with the *Python 3 (plasticome)*
kernel and Run All. Every step executes from the raw input; assertions guard the known-good
counts (42 / 484 / 46); the Results section rebuilds the report tables and figures.

**Shell / headless:**
```bash
./run.sh                       # uses the plasticome env by default
PY=/path/to/python ./run.sh    # or point at any Python 3.9+
```

**Individual steps** (each is an independent CLI, `--help` for options):
```bash
python scripts/step1_nodes.py --tsv data/cleaned_pazy_final.tsv --outprefix outputs/combined
python scripts/step23_graph.py --prefix outputs/combined --outdir outputs --tag petadex --date 2026-07-21
```

## Method (paper-faithful, carried forward v4–v8; do not change without a version bump)

- Node identity = md5 of `normalize()` (letters only, uppercased) — the paper's component definition.
- Aligner: `usearch v11.0.667 -allpairs_local -acceptall`, **both FASTA orientations**.
- HSP selection: best per unordered pair by **max bits** (ties → lower e-value) — the B4 convention.
- Edge criterion (post-filter only, never a search parameter): **≥30% aaid AND e-value < 1e-5**,
  thresholded as reported at this run's db (no down-scaling).
- Partition: undirected **single-linkage** connected components over all nodes.
- v1 `component`/`cath` and the input's prior `component_id` are **overlays only** — never seed
  or gate the partition.

## Outputs (`outputs/`)

| file | what |
|---|---|
| `component_assignment_petadex_<date>.csv` | one row per node → component |
| `component_members_petadex_<date>.csv`    | one row per component (size, v1 spread, CATH, members) |
| `component_edges_petadex_<date>.csv`      | passing edges |
| `cleaned_pazy_final_components_<date>.csv`| every source row + joined component + `prior_component_id` |
| `combined{,_rev}.fasta`, `combined{,_rev}_pairs.tsv` | node set + raw usearch HSPs |
| `stats_petadex_<date>.json`, `combined_stats.json`   | provenance + counts |
| `step0_check/`                            | the 213-node engine-sanity inputs/pairs/partition (→ 42) |

The older `../v2` … `../v8` iteration folders are left untouched; this folder is the
consolidated, canonical entry point.
