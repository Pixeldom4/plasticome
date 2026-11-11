# CURRENT

Each entry here is a symlink to the folder a figure script writes into, so it
always shows the latest render in every format (svg, pdf, png, html). Nothing in
this folder is a copy. Open a figure from here, edit the scripts where they live.

| Link | Points to | Written by |
|---|---|---|
| `tree-mechanism-groups/` | `tree/runs/all-sequences/results/figures/mechanism-groups/` | `tree/18_ring_tree.py` |
| `tree-fold-groups/` | `tree/runs/all-sequences/results/figures/fold-groups/` | `tree/18_ring_tree.py` |
| `tree-plain/` | `tree/runs/all-sequences/results/figures/plain/` | `tree/18_ring_tree.py` |
| `tree-centroids-plain/` | `tree/runs/centroids/results/figures/plain/` | `tree/18_ring_tree.py`, `PLASTICOME_TREE_SET=centroids` |
| `mechanism-cards/` | `tree/runs/mechanisms/` | `tree/17_mechanism_cards.py` |
| `mechanism-cards-grouped/` | `tree/runs/mechanisms-mechanism/` | `tree/17_mechanism_cards.py`, called from step 18 |
| `pie-hand-edited/` | `tree/figure-pie-hand-edited/` | by hand, then steps 13 and 14 |
| `heatmap-v2/` | `heatmap/out_v2/` | `heatmap/07_heatmap_v2.py` |
| `heatmap-centroids/` | `heatmap/out_centroids/` | `heatmap/02` to `05`, centroid set |
| `heatmap-all/` | `heatmap/out/` | `heatmap/02` to `05`, all sequences |
| `ssn/` | `ssn/out/` | `ssn/build_ssn.py`, `build_mst.py`, `build_aii.py` |
| `activity-probe/` | `activity-probe/out/` | `activity-probe/02_probe.py` |
| `palette/` | `palette/` | `tree/16_palette_reference.py` |

When a script starts writing somewhere new, add a link here instead of copying
files in.
