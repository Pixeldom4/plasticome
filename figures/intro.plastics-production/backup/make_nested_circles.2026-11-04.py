"""Nested proportional circles of carbon mass (Gt C): life vs plastics.

Same values and colours as 1a.plastics-production.svg. Life groups are nested
circles that all touch the baseline at one point; plastics sits beside the
nest on the same baseline, split into "landfill, environment or in use" and
"recycled or incinerated" slices (notebook numbers, plastics_numbers.py). Circle and slice areas are proportional to mass. Units are mm.

Run: python make_nested_circles.py  ->  1a.plastics-production.circles.svg
"""
import math
from pathlib import Path

OUT = Path(__file__).with_name("1a.plastics-production.circles.svg")

# name, Gt, fill, stroke
import plastics_numbers as N

LIFE = [
    ("Plants",   N.LIFE_C["Plants"],   "#7de8a1", "#24ca5b"),
    ("Bacteria", N.LIFE_C["Bacteria"], "#c59bed", "#a15ee3"),
    ("Fungi",    N.LIFE_C["Fungi"],    "#f0c5a6", "#e59f6c"),
    ("Archaea",  N.LIFE_C["Archaea"],  "#ecde97", "#ddc23c"),
    ("Protists", N.LIFE_C["Protists"], "#7cd9de", "#33c0c8"),
    ("Animals",  N.LIFE_C["Animals"],  "#f2b3b4", "#e1595b"),
]
# Plastics, notebook numbers for 2026 in Gt of carbon (see plastics_numbers.py).
# The circle is all plastic carbon ever produced; slices split it into what
# remains (landfill, environment, in use) and what was recycled or incinerated.
PLASTIC = ("Plastics", N.PRODUCED_C, "#7db3e8", "#1b5998")
# name, Gt C, fill (colours from the notebook's production plot)
PLASTIC_PARTS = [
    ("landfill, environment or in use", N.REMAINING_C, "#7db3e8"),
    ("recycled or incinerated", N.RECYCLED_C + N.INCINERATED_C, "#1b5998"),
]
KEY = (1, "#aeaeae", "#989898")

FONT = "'Ubuntu Sans Mono', monospace"
HEAD_FONT = "'Ubuntu Sans Mono'"   # panel letter + title, same string as the Inkscape panels
HEAD_SIZE = 7.0          # letter and title share one size
HEAD_GAP = 9.2           # title starts this far right of the letter, mm
INK = "#000000"
MUTED = "#555555"

R_MAX = 60.0                        # radius of the plants circle, mm
K = R_MAX / math.sqrt(LIFE[0][1])   # mm per sqrt(Gt)
rad = lambda v: K * math.sqrt(v)

LABEL_W = 40.0                      # room for the axis labels on the left
X0, Y0 = 12.0, 32.0                 # page origin of the panel
CX = X0 + LABEL_W + 6.0 + R_MAX     # shared centre line of the nest (6 mm clear of the axis)
BASE = Y0 + 2 * R_MAX               # shared baseline
SW = 0.5
GAP = 6.0                           # space between plants and plastics

parts = []
add = parts.append
fmt = lambda v: (f"{v:.1f}".rstrip("0").rstrip(".")) + " Gt C"


def label(x, y, name, value, anchor="start", size=4.2):
    add(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-family="{FONT}" fill="{INK}">'
        f'<tspan font-size="{size}" font-weight="bold">{name}</tspan>'
        f'<tspan font-size="{size*0.82:.2f}" fill="{MUTED}" dx="1.2">{value}</tspan></text>')


add(f'<text x="{X0:.2f}" y="{Y0-4.6:.2f}" font-family="{HEAD_FONT}" font-size="{HEAD_SIZE}" '
    f'font-weight="bold" fill="{INK}">A.</text>')
add(f'<text x="{X0+HEAD_GAP:.2f}" y="{Y0-4.6:.2f}" font-family="{HEAD_FONT}" font-size="{HEAD_SIZE}" '
    f'font-weight="bold" fill="{INK}">The Scale of Global Plastics Production</text>')

for name, v, fill, stroke in LIFE:
    r = rad(v)
    add(f'<circle cx="{CX:.3f}" cy="{BASE-r:.3f}" r="{r:.3f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{SW}"/>')

# plastics beside the nest, same baseline
rp = rad(PLASTIC[1])
pcx = CX + R_MAX + GAP + rp
pcy = BASE - rp
# wedges: area of each slice is proportional to its mass, clockwise from 12 o'clock
tot = sum(v for _, v, _ in PLASTIC_PARTS)
a0 = -math.pi / 2
wedge_mid = []
for name, v, fill in PLASTIC_PARTS:
    a1 = a0 + 2 * math.pi * v / tot
    x0, y0 = pcx + rp * math.cos(a0), pcy + rp * math.sin(a0)
    x1, y1 = pcx + rp * math.cos(a1), pcy + rp * math.sin(a1)
    big = 1 if a1 - a0 > math.pi else 0
    add(f'<path d="M{pcx:.3f},{pcy:.3f} L{x0:.3f},{y0:.3f} '
        f'A{rp:.3f},{rp:.3f} 0 {big} 1 {x1:.3f},{y1:.3f} Z" fill="{fill}" '
        f'stroke="#ffffff" stroke-width="0.35" stroke-linejoin="round"/>')
    wedge_mid.append((a0 + a1) / 2)
    a0 = a1
add(f'<circle cx="{pcx:.3f}" cy="{pcy:.3f}" r="{rp:.3f}" fill="none" '
    f'stroke="{PLASTIC[3]}" stroke-width="{SW*1.6}"/>')

# y axis: one tick per circle at the height of its top, label in the group colour
def darken(hex_, f=0.72):
    """Stroke colour pulled toward black so it reads as text on white."""
    r, g, b = (int(hex_[i:i + 2], 16) for i in (1, 3, 5))
    return "#%02x%02x%02x" % tuple(round(c * f) for c in (r, g, b))


AX = X0 + LABEL_W                    # axis x
TICK = 1.6
add(f'<line x1="{AX:.2f}" y1="{BASE:.2f}" x2="{AX:.2f}" y2="{Y0:.2f}" '
    f'stroke="{INK}" stroke-width="0.35"/>')
add(f'<line x1="{AX-TICK:.2f}" y1="{BASE:.2f}" x2="{AX:.2f}" y2="{BASE:.2f}" '
    f'stroke="{INK}" stroke-width="0.35"/>')
add(f'<text x="{AX-TICK-1.2:.2f}" y="{BASE+1.2:.2f}" text-anchor="end" '
    f'font-family="{FONT}" font-size="3.0" fill="{MUTED}">0</text>')

# (name, value, stroke, circle-top x, height)
ticks = [(n, v, s, CX, 2 * rad(v)) for n, v, _, s in LIFE]
ticks.append((PLASTIC[0], PLASTIC[1], PLASTIC[3], pcx, 2 * rp))
ticks.sort(key=lambda t: t[4])

# spread label baselines so no two sit closer than MIN_GAP mm
MIN_GAP, FS = 4.3, 3.4
ly = []
for t in ticks:
    y = max(t[4], (ly[-1] + MIN_GAP) if ly else 4.0)
    ly.append(y)

for (name, v, stroke, tx, h), lh in zip(ticks, ly):
    col = darken(stroke)
    ty = BASE - h
    # dashed guide from axis to the top of the circle
    add(f'<line x1="{AX:.2f}" y1="{ty:.2f}" x2="{tx:.2f}" y2="{ty:.2f}" '
        f'stroke="{stroke}" stroke-width="0.3" stroke-dasharray="1 0.8"/>')
    add(f'<circle cx="{tx:.2f}" cy="{ty:.2f}" r="0.5" fill="{col}"/>')
    # tick, then a short bend to the (possibly nudged) label
    lx = AX - TICK - 1.8
    add(f'<polyline points="{AX:.2f},{ty:.2f} {AX-TICK:.2f},{ty:.2f} '
        f'{lx+0.6:.2f},{BASE-lh:.2f}" fill="none" stroke="{col}" stroke-width="0.35"/>')
    add(f'<text x="{lx:.2f}" y="{BASE-lh+FS*0.35:.2f}" text-anchor="end" '
        f'font-family="{FONT}" fill="{col}">'
        f'<tspan font-size="{FS}" font-weight="bold">{name}</tspan>'
        f'<tspan font-size="{FS*0.85:.2f}" dx="1.2">{fmt(v)}</tspan></text>')

# breakdown labels beneath the plastics circle, right-aligned to its right edge
# so they run left into the free space under the nest, in each slice's colour
bx = pcx + rp
for i, (name, v, fill) in enumerate(PLASTIC_PARTS):
    by = BASE + 5.5 + i * 4.6
    col = darken(fill, 0.6 if name.startswith("landfill") else 0.85)
    add(f'<circle cx="{bx-1.1:.2f}" cy="{by-1.0:.2f}" r="1.1" fill="{fill}"/>')
    add(f'<text x="{bx-3.4:.2f}" y="{by:.2f}" text-anchor="end" font-family="{FONT}" fill="{col}">'
        f'<tspan font-size="3.2" font-weight="bold">{name}</tspan>'
        f'<tspan font-size="2.9" dx="1.2">{fmt(v)}</tspan></text>')

# axis title
add(f'<text transform="translate({X0+1.5:.2f},{Y0+R_MAX:.2f}) rotate(-90)" '
    f'text-anchor="middle" font-family="{FONT}" font-size="3.4" fill="{INK}">'
    f'circle diameter (area \u221d carbon mass)</text>')

# key
kv, kfill, kstroke = KEY
kr = rad(kv)
ky = BASE + 6.0
add(f'<circle cx="{X0+kr:.3f}" cy="{ky+kr:.3f}" r="{kr:.3f}" '
    f'fill="{kfill}" stroke="{kstroke}" stroke-width="{SW}"/>')
add(f'<text x="{X0+2*kr+2:.2f}" y="{ky+kr+1.2:.2f}" font-family="{FONT}" '
    f'font-size="3.4" fill="{INK}">= 1 billion tonnes of carbon (area)</text>')

W = bx + 8
H = ky + 2 * kr + 10
svg = (f'<?xml version="1.0" encoding="UTF-8"?>\n'
       f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}mm" height="{H:.1f}mm" '
       f'viewBox="0 0 {W:.1f} {H:.1f}">\n'
       f'<rect width="100%" height="100%" fill="#ffffff"/>\n'
       + "\n".join(parts) + "\n</svg>\n")
OUT.write_text(svg)
print(OUT, f"{W:.1f} x {H:.1f} mm")
for name, v, *_ in LIFE + [PLASTIC]:
    print(f"  {name:9s} {v:7.3f} Gt C  radius {rad(v):6.2f} mm")
