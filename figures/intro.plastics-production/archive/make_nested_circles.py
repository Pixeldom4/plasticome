"""Nested proportional circles of carbon mass (Gt C, 2026): life vs plastics.

Life groups are nested circles touching the baseline at one point. Plastics
sits beside the nest on the same baseline. Its area is all plastic ever
produced, split into three slices: landfill, environment or in use; recycled;
incinerated (notebook colours). Circle and slice areas are proportional to
carbon mass. Numbers from plastics_numbers.py. Units are mm.

Previous version (two slices, y axis): backup/make_nested_circles.2026-11-04.py

Run: python make_nested_circles.py  ->  1a.plastics-production.circles.svg
"""
import math
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # plastics_numbers.py lives one level up
import plastics_numbers as N

OUT = Path(__file__).with_name("1a.plastics-production.circles.svg")

# name, Gt C, fill, stroke
LIFE = [
    ("Plants",   N.LIFE_C["Plants"],   "#7de8a1", "#24ca5b"),
    ("Bacteria", N.LIFE_C["Bacteria"], "#c59bed", "#a15ee3"),
    ("Fungi",    N.LIFE_C["Fungi"],    "#f0c5a6", "#e59f6c"),
    ("Archaea",  N.LIFE_C["Archaea"],  "#ecde97", "#ddc23c"),
    ("Protists", N.LIFE_C["Protists"], "#7cd9de", "#33c0c8"),
    ("Animals",  N.LIFE_C["Animals"],  "#f2b3b4", "#e1595b"),
]
P_STROKE = "#1b5998"
# name, Gt C, fill. One blue hue at three well-separated lightness levels
# (CIE L* 71, 37, 90) so the slices stay distinct in greyscale print.
# Incinerated is palest: that plastic is gone.
SLICES = [
    ("Landfill, Environment or In-Use", N.REMAINING_C, "#7db3e8"),
    ("Recycled", N.RECYCLED_C, "#1b5998"),
    ("Incinerated", N.INCINERATED_C, "#d2e3f5"),
]

FONT = "'Ubuntu Sans Mono', monospace"
HEAD_FONT = "'Ubuntu Sans Mono'"   # panel letter + title, same string as the Inkscape panels
HEAD_SIZE, HEAD_GAP = 7.0, 9.2     # letter and title share one size
INK, MUTED = "#000000", "#555555"

R_MAX = 60.0                        # radius of the plants circle, mm
K = R_MAX / math.sqrt(LIFE[0][1])   # mm per sqrt(Gt C)
rad = lambda v: K * math.sqrt(v)

LABEL_W = 47.0                      # room for the label column on the left
X0, Y_TITLE = 12.0, 12.0
Y0 = Y_TITLE + 7.0                  # top of the plants circle
CX = X0 + LABEL_W + 6.0 + R_MAX     # shared centre line of the nest
BASE = Y0 + 2 * R_MAX               # shared baseline
SW, GAP = 0.5, 6.0

parts = []
add = parts.append
fmt = lambda v: (f"{v:.1f}".rstrip("0").rstrip(".")) + " Gt C"


def darken(hex_, f=0.72):
    """Stroke colour pulled toward black so it reads as text on white."""
    r, g, b = (int(hex_[i:i + 2], 16) for i in (1, 3, 5))
    return "#%02x%02x%02x" % tuple(round(c * f) for c in (r, g, b))


# header
add(f'<text x="{X0}" y="{Y_TITLE}" font-family="{HEAD_FONT}" font-size="{HEAD_SIZE}" '
    f'font-weight="bold" fill="{INK}">A.</text>')
add(f'<text x="{X0+HEAD_GAP}" y="{Y_TITLE}" font-family="{HEAD_FONT}" font-size="{HEAD_SIZE}" '
    f'font-weight="bold" fill="{INK}">The Scale of Global Plastics Production</text>')

# life circles, largest first
for name, v, fill, stroke in LIFE:
    r = rad(v)
    add(f'<circle cx="{CX:.3f}" cy="{BASE-r:.3f}" r="{r:.3f}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{SW}"/>')

# plastics circle: area = all ever produced, slices clockwise from 12 o'clock
rp = rad(N.PRODUCED_C)
pcx, pcy = CX + R_MAX + GAP + rp, BASE - rp
pt = lambda a: (pcx + rp * math.cos(a), pcy + rp * math.sin(a))
a0 = -math.pi / 2
for name, v, fill in SLICES:
    a1 = a0 + 2 * math.pi * v / N.PRODUCED_C
    (x0, y0), (x1, y1) = pt(a0), pt(a1)
    big = 1 if a1 - a0 > math.pi else 0
    add(f'<path d="M{pcx:.3f},{pcy:.3f} L{x0:.3f},{y0:.3f} A{rp:.3f},{rp:.3f} 0 {big} 1 '
        f'{x1:.3f},{y1:.3f} Z" fill="{fill}" stroke="#ffffff" stroke-width="0.35" '
        f'stroke-linejoin="round"/>')
    a0 = a1
add(f'<circle cx="{pcx:.3f}" cy="{pcy:.3f}" r="{rp:.3f}" fill="none" '
    f'stroke="{P_STROKE}" stroke-width="{SW*1.6}"/>')

# label column: one leader per circle top, label in the group colour
AX = X0 + LABEL_W
ticks = [(n, fmt(v), s, CX, 2 * rad(v)) for n, v, _, s in LIFE]
ticks.append(("Plastics", fmt(N.PRODUCED_C), P_STROKE, pcx, 2 * rp))
ticks.sort(key=lambda t: t[4])
MIN_GAP, FS = 4.8, 3.4              # label baselines no closer than MIN_GAP mm
ly = []
for t in ticks:
    ly.append(max(t[4], (ly[-1] + MIN_GAP) if ly else 4.0))
for (name, val, stroke, tx, h), lh in zip(ticks, ly):
    col = darken(stroke)
    ty = BASE - h
    # solid guide at the true height of the circle top
    add(f'<line x1="{AX:.2f}" y1="{ty:.2f}" x2="{tx:.2f}" y2="{ty:.2f}" stroke="{stroke}" '
        f'stroke-width="0.3"/>')
    add(f'<circle cx="{tx:.2f}" cy="{ty:.2f}" r="0.5" fill="{col}"/>')
    # pointer to the label: solid if the label sits at the true height,
    # dashed if it was nudged to avoid a neighbour
    lx = AX - 5.0
    dash = "" if abs(lh - h) < 0.05 else ' stroke-dasharray="0.6 0.45"'
    add(f'<line x1="{AX:.2f}" y1="{ty:.2f}" x2="{lx+0.8:.2f}" y2="{BASE-lh:.2f}" '
        f'stroke="{col}" stroke-width="0.3"{dash}/>')
    add(f'<text x="{lx:.2f}" y="{BASE-lh+FS*0.35:.2f}" text-anchor="end" font-family="{FONT}" '
        f'fill="{col}"><tspan font-size="{FS}" font-weight="bold">{name}</tspan>'
        f'<tspan font-size="{FS*0.85:.2f}" dx="1.2">{val}</tspan></text>')

# slice legend under the plastics circle, right-aligned, dark ink; colour on the swatches
bx = pcx + rp
by = BASE + 6.0
for i, (name, v, fill) in enumerate(SLICES):
    y = by + i * 4.6
    add(f'<rect x="{bx-2.6:.2f}" y="{y-2.5:.2f}" width="2.6" height="2.6" rx="0.4" fill="{fill}" '
        f'stroke="{P_STROKE}" stroke-width="0.3"/>')
    add(f'<text x="{bx-4.2:.2f}" y="{y:.2f}" text-anchor="end" font-family="{FONT}" fill="{INK}">'
        f'<tspan font-size="3.2">{name}</tspan>'
        f'<tspan font-size="2.9" fill="{MUTED}" dx="1.2">{fmt(v)}</tspan></text>')
legend_bottom = by + (len(SLICES) - 1) * 4.6

# key
kr = rad(1)
ky = BASE + 6.0
add(f'<circle cx="{X0+kr:.3f}" cy="{ky+kr:.3f}" r="{kr:.3f}" fill="#aeaeae" stroke="#989898" '
    f'stroke-width="{SW}"/>')
add(f'<text x="{X0+2*kr+2:.2f}" y="{ky+kr+1.2:.2f}" font-family="{FONT}" font-size="3.4" '
    f'fill="{INK}">= 1 billion tonnes of carbon</text>')

W = bx + 8
H = max(legend_bottom, ky + 2 * kr) + 8
OUT.write_text(f'<?xml version="1.0" encoding="UTF-8"?>\n'
               f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}mm" height="{H:.1f}mm" '
               f'viewBox="0 0 {W:.1f} {H:.1f}">\n<rect width="100%" height="100%" fill="#ffffff"/>\n'
               + "\n".join(parts) + "\n</svg>\n")
print(OUT, f"{W:.1f} x {H:.1f} mm")
for name, v, *_ in LIFE:
    print(f"  {name:9s} {v:7.3f} Gt C  radius {rad(v):6.2f} mm")
for name, v, _ in SLICES:
    print(f"  plastics, {name}: {v:.3f} Gt C")
