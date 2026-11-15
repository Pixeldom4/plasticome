#!/usr/bin/env python3
"""
The working vertical mechanism tree (promoted from the demo on 2026-11-14),
sized for print as a page-wide figure.

Run 17b_mechanism_colour_renders.py --spread first (the structure renders).
This runs 18_ring_tree.py with:

  views       each structure centred on its catalytic residues and turned so they
              do not overlap on screen (17b --spread renders, via --pictures)
  structures  each protein scaled until it fills its column or its row, captions
              aligned along the row, the number beside the title, no structure ID,
              the catalytic residues named on the protein with a thin line to each
              (--fill-cells --fit-content --no-ids --label-residues); no close ups
  order       groups from the most tips to the fewest, numbered in that order
              (--order-by-size)
  tree side   only the numbered circles beside the tree, spaced apart with a
              tick back to their bar when moved (--no-group-names); no tip dots
              (--no-tips); plastic columns in their own colours (one grey,
              --mono-plastics, was tried and reverted 2026-11-14)
  height      the tree's rows shortened until the structure rows have no spare
              height (TIGHT, --row-h): the structures keep their size and the
              figure gets shorter, since the tree is read by group, not by branch
  print       FONT throughout; text and branch widths set so that at PRINT_MM
              wide the text is TEXT_PT and the branches at least LINE_PT
              (--text-pt, --line-scale), re-rendered until both hold

Writes runs/all-sequences/results/figures/mechanism-groups/tree-vertical/tree-vertical.{png,svg,pdf,html},
or with --insets the parallel demo (a round close up of every active site at the
top right of its column) to .../mechanism-groups/demo/.
Check that SVG for Inkscape edits (grep -c inkscape) before running: this overwrites it.

Usage (plasticome env): python 18b_vertical_demo.py [--insets]
"""

import argparse
import os
import re
import subprocess
import sys

TREE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(TREE, "runs", "all-sequences", "results", "figures", "mechanism-groups")
PICTURES = os.path.join(TREE, "runs", "mechanisms-mechanism-grey-spread")   # 17b --spread

PRINT_MM = 180                 # double-column width the figure is sized for (journal not chosen)
TEXT_PT, LINE_PT = 7.0, 0.25   # smallest text and branch width at that width
FONT = "Arial"
TIGHT = True                   # shorten the tree's rows until the structure rows have no spare height


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--insets", action="store_true",
                    help="the parallel demo: the same figure with a round close up of each "
                         "active site, written to demo/ instead of tree-vertical/")
    args = ap.parse_args()
    OUT = os.path.join(FIGS, "demo" if args.insets else "tree-vertical")
    os.makedirs(OUT, exist_ok=True)
    extra = ["--insets"] if args.insets else []
    vfs, ls, row_h = 37.5, 1.0, 0.03
    for _ in range(8):
        out = subprocess.run([sys.executable, "18_ring_tree.py", "--vertical", "--type-bars",
                              "--grouping", "mechanism", "--fill-cells", "--fit-content",
                              "--no-group-names", "--no-tips", "--no-ids", "--label-residues", "--order-by-size",
                              "--font", FONT, *extra, "--row-h", f"{row_h:.5f}",
                              "--text-pt", f"{vfs:.2f}", "--line-scale", f"{ls:.3f}",
                              "--pictures", PICTURES, "--stem", os.path.join(OUT, "tree-vertical")],
                             cwd=TREE, check=True, capture_output=True, text=True).stdout
        w_in, h_in = map(float, re.search(r"canvas ([\d.]+) x ([\d.]+) in", out).groups())
        need, have, tree_h = map(float, re.search(
            r"fit rows: need ([\d.]+) in, have ([\d.]+) in, tree ([\d.]+) in", out).groups())
        scale = PRINT_MM / 25.4 / w_in
        print(f"  {vfs:.1f} pt text, lines x{ls:.2f}, rows {row_h:.4f} in: canvas {w_in:.1f} x "
              f"{h_in:.1f} in, prints {PRINT_MM} x {h_in * scale * 25.4:.0f} mm at {vfs * scale:.1f} pt, "
              f"structure rows need {need:.1f} of {have:.1f} in")
        new_vfs, new_ls = TEXT_PT / scale, max(1.0, LINE_PT / (0.55 * scale))
        # the tree shorter by the height the two structure rows do not use
        new_row_h = max(0.01, row_h * (tree_h - 2 * (have - need)) / tree_h) if TIGHT else row_h
        if abs(new_vfs - vfs) < 0.5 and abs(new_ls - ls) < 0.05 and abs(new_row_h - row_h) < 0.0005:
            break
        vfs, ls, row_h = new_vfs, new_ls, new_row_h
    print("".join(l + "\n" for l in out.splitlines()
                  if "canvas" in l or "layout" in l or "inset" in l), end="")


if __name__ == "__main__":
    main()
