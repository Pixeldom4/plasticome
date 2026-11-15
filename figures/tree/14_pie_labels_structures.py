#!/usr/bin/env python3
"""
Step 14: the pie with the structures beside the fold labels.

A second arrangement of the same material as step 13, which puts the
structures in a strip along the bottom. Here they run up the margins instead,
in a column beside the fold labels and in the same order, so each structure
sits next to the label of the fold it stands for. Both files are kept; this
one does not replace step 13's.

Why a column beside the labels rather than one structure per label. The fold
labels are 23 points apart at the top of the stack and span 560 points in all,
while ten legible structures need about 900. There is no way to sit each one
on its own label without either shrinking them to the height of a line of type
or pulling the labels apart, and the labels were placed by hand. So the
structures take the margin: the left column holds the nine folds whose labels
are on the left, top to bottom in label order, and the right column holds the
one fold labelled on the right. Each is as close to its own label as the
others allow, by the same `spread` the pie uses for the labels themselves.

What ties a structure to its label is the tint, which is the slice's own, the
order, which is the label order, and the slice id under the structure. No
leader lines: a line from a margin column to a right aligned label has to
cross the labels below it, which costs more than it explains.

The original is copied to a backup and never edited in place. The structures
go in their own Inkscape layer, `structures`, so the columns move or delete as
one object.

Writes figure-pie-hand-edited/figure_pie_labelled_structures.{svg,png}.
"""

import importlib
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

from common import OUT, read_tsv

pie = importlib.import_module("08_draw_pie")        # for spread()

HERE = os.path.dirname(os.path.abspath(__file__))
PIE = os.path.join(HERE, "figure-pie-hand-edited", "figure_pie.svg")
BACKUP = os.path.join(HERE, "figure-pie-hand-edited", "backups",
                      "figure_pie_backup_before_structures.svg")
NEW = os.path.join(HERE, "figure-pie-hand-edited", "figure_pie_labelled_structures.svg")
SDIR = os.path.join(HERE, "..", "structures", "fold-representatives", "structures")

MARGIN = 100.0         # points of canvas added at each side
BOX = 68.0             # the square a structure is drawn in
CAPTION = 19.0         # room under a box for its two caption lines
GAP = 7.0              # between one item and the next in a column
EDGE = 15.0            # smallest gap from a column to the canvas edge
# Nine items at BOX + CAPTION + GAP do not fit the pie's own height, so the
# canvas grows by this much at the bottom. The pie stays where it is.
HEIGHT_ADD = 70.0
BG = "#fcfcfb"
FONT = ("font-family: 'DejaVu Sans', 'Bitstream Vera Sans', 'Computer Modern "
        "Sans Serif', 'Lucida Grande', 'Verdana', 'Geneva', 'Lucid', 'Arial', "
        "'Helvetica', 'Avant Garde', sans-serif")
SVG = "http://www.w3.org/2000/svg"
INK = "http://www.inkscape.org/namespaces/inkscape"


def shrink(svg, prefix):
    """Namespace the ids and drop pointless coordinate precision."""
    for i in sorted(set(re.findall(r'id="([^"]+)"', svg)), key=len, reverse=True):
        svg = svg.replace(f'id="{i}"', f'id="{prefix}{i}"')
        svg = svg.replace(f'url(#{i})', f'url(#{prefix}{i})')
        svg = svg.replace(f'href="#{i}"', f'href="#{prefix}{i}"')
    return re.sub(r'(-?\d+\.\d{3,})',
                  lambda m: f"{float(m.group(1)):.2f}".rstrip("0").rstrip("."), svg)


def inner(svg):
    vb = re.search(r'viewBox="([^"]+)"', svg).group(1)
    body = svg[svg.index(">", svg.index("<svg")) + 1:svg.rindex("</svg>")]
    return re.sub(r"<metadata>.*?</metadata>", "", body, flags=re.S), vb


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def label_anchors(path):
    """slice -> (y of its fold label, which side it is on)."""
    root = ET.parse(path).getroot()
    out = {}
    for g in root.iter(f"{{{SVG}}}g"):
        gid = g.get("id") or ""
        if not re.fullmatch(r"label-g\d\d", gid):
            continue
        t = g.find(f".//{{{SVG}}}text")
        side = "left" if "text-anchor: end" in (t.get("style") or "") else "right"
        out[gid.replace("label-", "")] = (float(t.get("y")), side)
    return out


def column(slices, anchors, height, item):
    """y of each box top, as near its label as the column allows."""
    order = sorted(slices, key=lambda g: anchors[g][0])
    # spread() is written for matplotlib axes, where y runs up and a list read
    # top to bottom is decreasing. SVG y runs down, so it goes in negated.
    ys = [-u for u in pie.spread([-anchors[g][0] for g in order],
                                 [item] * len(order))]
    top = [y - item / 2 for y in ys]
    # spread keeps the centre of mass, which can run off the canvas when the
    # labels are bunched. Slide the whole column back inside, keeping its order
    # and its spacing.
    if top and top[0] < EDGE:
        top = [y + (EDGE - top[0]) for y in top]
    over = (top[-1] + item) - (height - EDGE) if top else 0
    if over > 0:
        top = [y - over for y in top]
    return list(zip(order, top))


def main():
    reps = {r["slice"]: r for r in read_tsv(os.path.join(OUT, "fold_representatives.tsv"))}
    src = {r["slice"]: r for r in read_tsv(os.path.join(SDIR, "sources.tsv"))}
    svgs = {f.split("_")[0]: os.path.join(SDIR, "svg", f)
            for f in os.listdir(os.path.join(SDIR, "svg")) if f.endswith(".svg")}

    pie_svg = open(PIE).read()
    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', pie_svg)
    W, H = float(m.group(1)), float(m.group(2))
    if not os.path.exists(BACKUP):
        shutil.copyfile(PIE, BACKUP)
        print(f"backup: {os.path.relpath(BACKUP, HERE)}")
    else:
        print(f"backup already there, left alone: {os.path.relpath(BACKUP, HERE)}")

    anchors = label_anchors(PIE)
    item = BOX + CAPTION + GAP
    H2 = H + HEIGHT_ADD
    placed = []
    for side, x in (("left", (MARGIN - BOX) / 2),
                    ("right", MARGIN + W + (MARGIN - BOX) / 2)):
        members = [g for g in reps if anchors[g][1] == side]
        need = len(members) * item + 2 * EDGE
        if need > H2:
            raise SystemExit(f"{side} column needs {need:.0f} pt of {H2:.0f}; "
                             f"lower BOX or raise HEIGHT_ADD")
        for g, y in column(members, anchors, H2, item):
            placed.append((g, x, y, side))

    W2 = W + 2 * MARGIN
    note = {"rcsb": "", "afdb": "AlphaFold", "esmfold": "ESMFold"}
    parts = [f'<rect x="0" y="0" width="{W2}" height="{H2}" fill="{BG}" />',
             f'<g transform="translate({MARGIN},0)">',
             pie_svg[pie_svg.index(">", pie_svg.index("<svg")) + 1:pie_svg.rindex("</svg>")],
             '</g>',
             '<g id="layer-structures" inkscape:groupmode="layer" '
             'inkscape:label="structures">']
    for g, x, y, side in placed:
        rep, s = reps[g], src[g]
        body, vb = inner(open(svgs[g]).read())
        parts.append(f'<g id="struct-{g}" inkscape:label="{esc(rep["representative"])}">')
        parts.append(f'<svg x="{x:.2f}" y="{y:.2f}" width="{BOX}" height="{BOX}" '
                     f'viewBox="{vb}" overflow="visible">'
                     + shrink(body, f"{g}-") + '</svg>')
        cx = x + BOX / 2
        tail = note[s["source"]] or rep["rep_pdb"]
        for dy, size, weight, colour, text in (
                (BOX + 9.5, 8.0, "bold", "#1a1a1a", f'{g}  {rep["representative"]}'),
                (BOX + 18.5, 6.8, "normal", "#6b6b6b",
                 f'{rep["catalytic_cath"]}  ({tail})')):
            parts.append(
                f'<text x="{cx:.2f}" y="{y + dy:.2f}" style="font-size: {size}px; '
                f'{FONT}; font-weight: {weight}; text-anchor: middle; fill: {colour}">'
                f'{esc(text)}</text>')
        parts.append('</g>')
    parts.append('</g>')

    head = pie_svg[:pie_svg.index(">", pie_svg.index("<svg")) + 1]
    head = re.sub(r'width="[\d.]+pt"', f'width="{W2:.6f}pt"', head)
    head = re.sub(r'height="[\d.]+pt"', f'height="{H2:.6f}pt"', head)
    head = re.sub(r'viewBox="0 0 [\d.]+ [\d.]+"', f'viewBox="0 0 {W2} {H2}"', head)
    if "xmlns:inkscape" not in head:
        head = head.replace("<svg ", f'<svg xmlns:inkscape="{INK}" ', 1)

    with open(NEW, "w") as fh:
        fh.write(head + "\n" + "\n".join(parts) + "\n</svg>\n")
    png = NEW[:-4] + ".png"
    subprocess.run(["rsvg-convert", "-w", str(int(W2 * 4)), "-o", png, NEW], check=True)

    print(f"wrote {os.path.relpath(NEW, HERE)}  "
          f"{os.path.getsize(NEW) / 1e6:.1f} MB, canvas {W2:.0f} x {H2:.0f} pt")
    print(f"wrote {os.path.relpath(png, HERE)}")
    for g, x, y, side in sorted(placed, key=lambda p: (p[3], p[2])):
        print(f"  {g}  {side:<5} box top {y:6.1f}  label y {anchors[g][0]:6.1f}")


if __name__ == "__main__":
    main()
