# File tree

Map of the repo. Paths are relative to the repo root. Generated and cache folders
(`__pycache__`, `node_modules`, `.DS_Store`, `.git`) are left out. The 2026-11-03
reorganisation that produced this layout is in [MIGRATIONS.md](MIGRATIONS.md).

## Quick links

To find a figure, start at `figures/CURRENT/`. Each entry there is a symlink to the
folder its script writes into, so it always holds the latest render.

| What | Path |
|---|---|
| Every current figure | `figures/CURRENT/` |
| Mechanism-group tree | `figures/CURRENT/tree-mechanism-groups/` |
| Fold-group tree | `figures/CURRENT/tree-fold-groups/` |
| Tree scripts (steps 01 to 18) | `figures/tree/` |
| Heatmap scripts | `figures/heatmap/` |
| SSN and MST scripts | `figures/ssn/` |
| Hand-annotated structures | `figures/structures/annotations/` |
| Feedback rounds | `figures/feedback/` |
| Curation workbook | `annotation/union/plasticome.v2.xlsx` |
| Pipeline entry point | `pipeline_run.bash` |

## Root

```
plasticome/
├── README.md
├── pipeline_run.bash        pipeline driver, steps 1 to 6
├── pipeline_run.config
├── lib/                     pipeline step code
├── bin/                     diamond, usearch11, usearch12 (Linux ELF, run through Docker)
├── source-data/             upstream inputs of record
├── runs/                    one folder per pipeline run, see runs/README.md
├── annotation/
├── figures/
├── docs/
└── archive/
```

`lib/`, `bin/`, `source-data/` and `runs/` stay at the root because the pipeline
takes the repo root as two levels above `lib/`, and mounts it into Docker.

## figures/

```
figures/
├── CURRENT/                          symlinks, one per figure (see CURRENT/README.md)
├── plasticome.v2 - ...-annotated (1).csv   table the tree and SSN read (corrected PDB IDs)
├── petadex-logan-workflow.jpg
├── tree/                             forest of trees, ring figure, structures
│   ├── 01_components.py ... 18_ring_tree.py
│   ├── common.py  folds.py  catalytic_struct.py  ring_scene.py
│   ├── templates/  iqtree-models/  molstar/  reference_pdb/
│   ├── figure-pie-hand-edited/
│   ├── dennis-build/
│   ├── archive/
│   └── runs/
│       ├── all-sequences/{intermediate,results}/
│       │                        results/figures/{mechanism-groups,fold-groups,plain}/
│       ├── centroids/{intermediate,results}/
│       ├── mechanisms/               mechanism cards
│       └── mechanisms-mechanism/     mechanism cards, grouped colours for step 18
├── heatmap/                          01_identity ... 07_heatmap_v2, out/, out_centroids/, out_v2/
├── ssn/                              build_ssn, build_mst, build_aii
├── activity-probe/
├── palette/                          written by tree/16_palette_reference.py
├── structures/
│   ├── annotations/                  hand-annotated folds (good_folds/)
│   ├── fold-representatives/         folding bundle from tree step 11, structures from step 12
│   └── fold-representatives.tar.gz
├── bundles/
│   └── centroid-tree/                self-contained cluster bundle (keeps the old layout inside)
└── feedback/
    └── 2026-10-27/                   REPORT.md, previews/
```

## annotation/

```
annotation/
├── union/                   annotate_union.py, plasticome.v2.xlsx, annotated union CSVs, CLAUDE.md
└── wiki/                    01 plasticome, 02 CATH Assignment, 03 Component Assignment
```

## docs/ and archive/

```
docs/
├── FILE_TREE.md             this file
├── MIGRATIONS.md            dated reorganisations
└── aws-changes.md

archive/
├── cache/                   NCBI, BLAST, paper caches the pipeline reuses
├── diamond-usearch-comparison/
├── figure-design/
├── usearch12-bugreport/
├── thomas/                  files shared with Thomas
├── v1-doi-assignment/       DOI resolution for v1.1
└── plasticome_centroids_final.2026-08-11.tsv
```
