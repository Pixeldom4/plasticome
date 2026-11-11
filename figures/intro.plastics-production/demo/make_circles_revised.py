"""Demo 3: the nested circles figure with the suggested changes.

1. plastics circle in three slices: landfill/environment/in use and recycled
   are filled (still plastic); incinerated is hatched with a dashed edge (now CO2)
2. callout: plastic that still exists vs all animals, with an animals circle at scale
3. plastics axis label reads "produced", since the circle is all plastic ever made
4. no axis line or "0"; the labels are a column with leaders
5. subtitle with year, unit and sources
6. legend text in dark ink, wider label spacing, less space above the title

Run: python demo/make_circles_revised.py  ->  demo/demo-3.circles-revised.svg
"""
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import plastics_numbers as N  # noqa: E402

OUT = HERE / "demo-3.circles-revised.svg"

LIFE = [
    ("Plants",   N.LIFE_C["Plants"],   "#7de8a1", "#24ca5b"),
    ("Bacteria", N.LIFE_C["Bacteria"], "#c59bed", "#a15ee3"),
    ("Fungi",    N.LIFE_C["Fungi"],    "#f0c5a6", "#e59f6c"),
    ("Archaea",  N.LIFE_C["Archaea"],  "#ecde97", "#ddc23c"),
    ("Protists", N.LIFE_C["Protists"], "#7cd9de", "#33c0c8"),
    ("Animals",  N.LIFE_C["Animals"],  "#f2b3b4", "#e1595b"),
]
P_STROKE = "#1b5998"
EXISTS_C = N.PRODUCED_C - N.INCINERATED_C          # still plastic: not incinerated
# name, Gt C, fill ("hatch" = gone), still plastic?
SLICES = [
    ("landfill, environment or in use", N.REMAINING_C, "#7db3e8", True),
    ("recycled", N.RECYCLED_C, "#1b5998", True),
    ("incinerated, now CO₂", N.INCINERATED_C, "hatch", False),
]

FONT = "'Ubuntu Sans Mono', monospace"
HEAD_FONT = "'Ubuntu Sans Mono'"
HEAD_SIZE, HEAD_GAP = 7.0, 9.2
INK, MUTED = "#000000", "#555555"
SUB = ("Carbon mass in 2026. Plastics: Geyer 2017, Pottinger 2024. "
       "Life: Bar-On et al. 2018.")

R_MAX = 60.0
K = R_MAX / math.sqrt(LIFE[0][1])
rad = lambda v: K * math.sqrt(v)

LABEL_W = 47.0
X0, Y_TITLE = 12.0, 12.0
Y0 = Y_TITLE + 12.0                     # top of the plants circle
CX = X0 + LABEL_W + 6.0 + R_MAX
BASE = Y0 + 2 * R_MAX
SW, GAP = 0.5, 6.0

parts = []
add = parts.append
fmt = lambda v: (f"{v:.1f}".rstrip("0").rstrip(".")) + " Gt C"


def darken(hex_, f=0.72):
    r, g, b = (int(hex_[i:i + 2], 16) for i in (1, 3, 5))
    return "#%02x%02x%02x" % tuple(round(c * f) for c in (r, g, b))


add('<defs><pattern id="hatch" patternUnits="userSpaceOnUse" width="1.4" height="1.4" '
    'patternTransform="rotate(45)"><rect width="1.4" height="1.4" fill="#ffffff"/>'
    f'<line x1="0" y1="0" x2="0" y2="1.4" stroke="{P_STROKE}" stroke-width="0.3" '
    'stroke-opacity="0.55"/></pattern></defs>')

# header
add(f'<text x="{X0}" y="{Y_TITLE}" font-family="{HEAD_FONT}" font-size="{HEAD_SIZE}" '
    f'font-weight="bold" fill="{INK}">A.</text>')
add(f'<text x="{X0+HEAD_GAP}" y="{Y_TITLE}" font-family="{HEAD_FONT}" font-size="{HEAD_SIZE}" '
    f'font-weight="bold" fill="{INK}">The Scale of Global Plastics Production</text>')
add(f'<text x="{X0+HEAD_GAP}" y="{Y_TITLE+5.2}" font-family="{FONT}" font-size="3.1" '
    f'fill="{MUTED}">{SUB}</text>')

# life circles
for name, v, fill, stroke in LIFE:
    r = rad(v)
    add(f'<circle cx="{CX:.3f}" cy="{BASE-r:.3f}" r="{r:.3f}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{SW}"/>')

# plastics circle: area = all ever produced, slices clockwise from 12 o'clock
rp = rad(N.PRODUCED_C)
pcx, pcy = CX + R_MAX + GAP + rp, BASE - rad(N.PRODUCED_C)
pt = lambda a: (pcx + rp * math.cos(a), pcy + rp * math.sin(a))
a0 = -math.pi / 2
arcs = []
for name, v, fill, exists in SLICES:
    a1 = a0 + 2 * math.pi * v / N.PRODUCED_C
    (x0, y0), (x1, y1) = pt(a0), pt(a1)
    big = 1 if a1 - a0 > math.pi else 0
    f = "url(#hatch)" if fill == "hatch" else fill
    add(f'<path d="M{pcx:.3f},{pcy:.3f} L{x0:.3f},{y0:.3f} A{rp:.3f},{rp:.3f} 0 {big} 1 '
        f'{x1:.3f},{y1:.3f} Z" fill="{f}" stroke="#ffffff" stroke-width="0.35" '
        f'stroke-linejoin="round"/>')
    arcs.append((a0, a1, big, exists))
    a0 = a1
# outer edge: solid where the plastic still exists, dashed where it was burned
for a0, a1, big, exists in arcs:
    (x0, y0), (x1, y1) = pt(a0), pt(a1)
    dash = "" if exists else ' stroke-dasharray="1.1 0.7"'
    add(f'<path d="M{x0:.3f},{y0:.3f} A{rp:.3f},{rp:.3f} 0 {big} 1 {x1:.3f},{y1:.3f}" '
        f'fill="none" stroke="{P_STROKE}" stroke-width="{SW*1.6}"{dash}/>')

# label column (no axis line): one leader per circle top
AX = X0 + LABEL_W
ticks = [(n, fmt(v), s, CX, 2 * rad(v)) for n, v, _, s in LIFE]
ticks.append(("Plastics produced", fmt(N.PRODUCED_C), P_STROKE, pcx, 2 * rp))
ticks.sort(key=lambda t: t[4])
MIN_GAP, FS = 4.8, 3.4
ly = []
for t in ticks:
    ly.append(max(t[4], (ly[-1] + MIN_GAP) if ly else 4.0))
for (name, val, stroke, tx, h), lh in zip(ticks, ly):
    col = darken(stroke)
    ty = BASE - h
    add(f'<line x1="{AX:.2f}" y1="{ty:.2f}" x2="{tx:.2f}" y2="{ty:.2f}" stroke="{stroke}" '
        f'stroke-width="0.3" stroke-dasharray="1 0.8"/>')
    add(f'<circle cx="{tx:.2f}" cy="{ty:.2f}" r="0.5" fill="{col}"/>')
    lx = AX - 3.4
    add(f'<polyline points="{AX:.2f},{ty:.2f} {AX-1.6:.2f},{ty:.2f} {lx+0.6:.2f},{BASE-lh:.2f}" '
        f'fill="none" stroke="{col}" stroke-width="0.35"/>')
    add(f'<text x="{lx:.2f}" y="{BASE-lh+FS*0.35:.2f}" text-anchor="end" font-family="{FONT}" '
        f'fill="{col}"><tspan font-size="{FS}" font-weight="bold">{name}</tspan>'
        f'<tspan font-size="{FS*0.85:.2f}" dx="1.2">{val}</tspan></text>')

# legend under the plastics circle, right-aligned, dark ink; colour on the swatches
bx = pcx + rp
by = BASE + 6.0
for i, (name, v, fill, _) in enumerate(SLICES):
    y = by + i * 4.6
    f = "url(#hatch)" if fill == "hatch" else fill
    dash = ' stroke-dasharray="0.6 0.4"' if fill == "hatch" else ""
    add(f'<rect x="{bx-2.6:.2f}" y="{y-2.5:.2f}" width="2.6" height="2.6" rx="0.4" fill="{f}" '
        f'stroke="{P_STROKE}" stroke-width="0.3"{dash}/>')
    add(f'<text x="{bx-4.2:.2f}" y="{y:.2f}" text-anchor="end" font-family="{FONT}" fill="{INK}">'
        f'<tspan font-size="3.2">{name}</tspan>'
        f'<tspan font-size="2.9" fill="{MUTED}" dx="1.2">{fmt(v)}</tspan></text>')

# callout: what still exists vs all animals, animals circle at the same scale
cy = by + 3 * 4.6 + 3.6
add(f'<line x1="{bx-78:.2f}" y1="{cy-3.6:.2f}" x2="{bx:.2f}" y2="{cy-3.6:.2f}" '
    f'stroke="#c8c8c8" stroke-width="0.25"/>')
ra = rad(N.LIFE_C["Animals"])
acx, acy = bx - ra, cy + ra + 0.4
add(f'<circle cx="{acx:.2f}" cy="{acy:.2f}" r="{ra:.2f}" fill="#f2b3b4" stroke="#e1595b" '
    f'stroke-width="{SW}"/>')
tx = bx - 2 * ra - 2.4
add(f'<text x="{tx:.2f}" y="{acy-0.6:.2f}" text-anchor="end" font-family="{FONT}" fill="{INK}">'
    f'<tspan font-size="3.4" font-weight="bold">Still plastic: {fmt(EXISTS_C)}</tspan></text>')
add(f'<text x="{tx:.2f}" y="{acy+3.6:.2f}" text-anchor="end" font-family="{FONT}" '
    f'font-size="3.1" fill="{darken(P_STROKE)}"><tspan font-weight="bold">'
    f'{EXISTS_C / N.LIFE_C["Animals"]:.1f}×</tspan> the carbon in all animals (circle)</text>')

# key
kr = rad(1)
ky = BASE + 6.0
add(f'<circle cx="{X0+kr:.3f}" cy="{ky+kr:.3f}" r="{kr:.3f}" fill="#aeaeae" stroke="#989898" '
    f'stroke-width="{SW}"/>')
add(f'<text x="{X0+2*kr+2:.2f}" y="{ky+kr+1.2:.2f}" font-family="{FONT}" font-size="3.4" '
    f'fill="{INK}">= 1 billion tonnes of carbon (area)</text>')

W = bx + 8
H = acy + ra + 8
OUT.write_text(f'<?xml version="1.0" encoding="UTF-8"?>\n'
               f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}mm" height="{H:.1f}mm" '
               f'viewBox="0 0 {W:.1f} {H:.1f}">\n<rect width="100%" height="100%" fill="#ffffff"/>\n'
               + "\n".join(parts) + "\n</svg>\n")
print(OUT, f"{W:.1f} x {H:.1f} mm")
