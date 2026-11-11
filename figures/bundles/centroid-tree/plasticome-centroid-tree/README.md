# plasticome-centroid-tree: bundle for the cluster

Finishes the 90 percent centroid tree of the plasticome figure. Packed on
2026-10-01 from `plasticome/logan-plasticome-figure/tree`.

## Contents

```
PROMPT.md            the prompt for Claude Code on the server
run_centroids.sh     runs the remaining steps, then verify.py
verify.py            checks the result
environment.yml      conda environment (linux-64), tool versions pinned
logan-plasticome-figure/
   plasticome.v2 - ... (1).csv      the union table, input
   ssn/                             SSN components and the all-vs-all hits, input
   tree/                            the pipeline, README, and its state:
      work/  out/                   finished all sequences run (read only here)
      work_centroids/               centroid run so far: 55 of 56 components
source-data/plasticome-curated-seeds/   the 13 seed enzymes, input
```

The directory layout mirrors the project on the laptop, because the scripts
find their inputs by relative path. Do not flatten it.

## On the server

```
tar xzf plasticome-centroid-tree.tar.gz
cd plasticome-centroid-tree
claude          # then give it PROMPT.md, for example: "Read PROMPT.md and do what it says."
```

Or by hand:

```
conda env create -f environment.yml && conda activate plasticome-tree
PLASTICOME_TREE_THREADS=32 ./run_centroids.sh
```

## Bringing the result back

On the server, from the bundle root:

```
tar czf centroid_results.tar.gz logan-plasticome-figure/tree/out_centroids logan-plasticome-figure/tree/work_centroids
```

On the laptop, unpack it over the project so the paths line up:

```
cd /Users/Pixel/Documents/projects/plasticome
tar xzf centroid_results.tar.gz
```

The figures are then in `logan-plasticome-figure/tree/out_centroids/`.
