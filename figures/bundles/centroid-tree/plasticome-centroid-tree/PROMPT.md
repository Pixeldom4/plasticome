# Task: finish the plasticome 90 percent centroid tree on this cluster

This folder was packed on a laptop where the last heavy step was too slow to
run. Your job is to run that step and the short steps after it here, check the
result, and pack the outputs so they can be copied back. You are not being
asked to change the analysis.

## What this is

A phylogenetics pipeline for plastic degrading enzymes (the "plasticome").
The sequences fall into 56 similarity network components. Each component gets
its own alignment ensemble and maximum likelihood tree, components with the
same fold are joined on a structure based backbone, and the result is drawn as
figures. `logan-plasticome-figure/tree/README.md` explains every step and why
it is done that way. Read its "Run" and "Two sets" sections before starting.

The pipeline has already been run in full for all 492 sequences (`tree/work/`,
`tree/out/`). It is now being rerun for the 411 sequences that are 90 percent
centroids (`tree/work_centroids/`, `tree/out_centroids/`). For that run, 55 of
56 component trees are finished and shipped. Only component 0 is missing: 225
sequences, 16 Muscle5 alignments, then an IQ-TREE search with 1000 ultrafast
bootstraps, then 100 FastTree trees. On a 10 core laptop that is about an
hour. Everything after it takes minutes.

## What to do

1. Create the environment from `environment.yml` (conda or mamba, channels
   conda-forge and bioconda) and activate it. If an exact version pin is not
   available for this platform, relax that one pin, and tell me which tool
   and which version you got.
2. Find out how jobs are meant to be run here. If there is a scheduler
   (Slurm, PBS, SGE, LSF), do not run the heavy step on a login node: submit
   `run_centroids.sh` as a batch job on one node. A sensible request is 16 to
   64 cores, 16 GB of memory, 4 hours. Set `PLASTICOME_TREE_THREADS` to the
   number of cores actually allocated, not the node total. If you cannot tell
   how jobs should be submitted here, ask me before submitting anything.
3. Run `./run_centroids.sh`. It sets `PLASTICOME_TREE_SET=centroids`, runs
   steps 2, 4, 5, 6, 7, 8 and 9, and then `verify.py`.
4. When it finishes, run `python verify.py` yourself and read its output. All
   checks must pass.
5. Pack the results from the bundle root:
   `tar czf centroid_results.tar.gz logan-plasticome-figure/tree/out_centroids logan-plasticome-figure/tree/work_centroids`
6. Report back, as described under "What to tell me".

## Things that matter

- **Do not change scientific settings.** No edits to models, bootstrap counts,
  thresholds, the number of alignments, or seeds, and do not skip a step to
  save time. The centroid run has to be comparable with the all sequences run
  that was made with these exact settings. If a step fails, fix the cause
  (environment, paths, resources). If the only way forward would change the
  method, stop and tell me instead.
- **Do not touch `tree/work/` or `tree/out/`.** They hold the finished all
  sequences run. The centroid run reads from `work/` (it reuses a structure
  cache there) and writes only to `work_centroids/` and `out_centroids/`.
- **Do not rerun step 1 or step 3.** Their outputs are shipped. Step 3 looks
  up structures online; every lookup it needs is already cached, but compute
  nodes often have no internet and there is no reason to risk it.
- **No network access is needed** for any step in `run_centroids.sh`.
- **Threads.** Muscle5 and the FastTree bootstraps scale with
  `PLASTICOME_TREE_THREADS`. IQ-TREE does not gain from many threads on an
  alignment this short (about 150 columns), so it is capped at 8 through
  `PLASTICOME_TREE_IQ_THREADS`. Leave that cap unless you measure a gain.
- **The run is resumable.** If the job is killed, submit it again. A finished
  alignment ensemble is reused, an interrupted IQ-TREE search continues from
  its checkpoint, and finished components are skipped. Do not delete
  `work_centroids/components/c00/` to "start clean" unless a file in it is
  corrupt.
- **Paths have a space in one file name** (the input CSV). Quote paths.
- The scripts find binaries on PATH when `PLASTICOME_TREE_BIN` does not point
  at them, which `run_centroids.sh` arranges. Python and the tools must come
  from the same activated environment.

## How to tell it worked

`verify.py` checks the essentials: 56 components, component 0 with 225 tips
and 16 alignments, 100 bootstrap trees, 411 tips in the final trees, one tree
per fold group, and all figure files present. It also prints a short summary
of component 0 and the backbone.

For orientation, the all sequences run gave these numbers for component 0
(284 tips): 148 of 1598 alignment columns kept, median UFBoot 96, median
ensemble support 52. The centroid numbers will differ, but if they are wildly
different (for example fewer than 50 columns kept, or ensemble support near
zero everywhere) something is wrong and I want to know.

Small numeric differences from a laptop run are expected. IQ-TREE results are
not bit identical across platforms and thread counts.

## What to tell me

- Whether every check in `verify.py` passed, and its summary lines.
- Tool versions actually installed, if any pin was relaxed.
- How it was run (scheduler, cores, wall time for step 2).
- Anything you had to change to make it run, with the reason.
- The path to `centroid_results.tar.gz`.
