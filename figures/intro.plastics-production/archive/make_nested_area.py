"""Nested proportional area chart of global carbon mass (Gt C): life vs plastics.

Numbers from plastics_numbers.py (notebook, 2026, Gt C); colours as 1a.plastics-production.svg. Life groups are nested
squares sharing the bottom-left corner; plastics sits in the bottom-right corner
of the plants square on the same baseline so it can be read against archaea.
Square area is proportional to mass. Units are mm.

Run: python make_nested_area.py  ->  1a.plastics-production.nested.svg
"""
import math
from pathlib import Path

OUT = Path(__file__).with_name("1a.plastics-production.nested.svg")

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # plastics_numbers.py lives one level up
import plastics_numbers as N

# name, Gt C, fill, stroke
LIFE = [
    ("Plants",   N.LIFE_C["Plants"],   "#7de8a1", "#24ca5b"),
    ("Bacteria", N.LIFE_C["Bacteria"], "#c59bed", "#a15ee3"),
    ("Fungi",    N.LIFE_C["Fungi"],    "#f0c5a6", "#e59f6c"),
    ("Archaea",  N.LIFE_C["Archaea"],  "#ecde97", "#ddc23c"),
    ("Protists", N.LIFE_C["Protists"], "#7cd9de", "#33c0c8"),
    ("Animals",  N.LIFE_C["Animals"],  "#f2b3b4", "#e1595b"),
]
# carbon in plastic still in landfill, the environment or in use (2026)
PLASTIC = ("Plastic remaining", N.REMAINING_C, "#7db3e8", "#1b5998")
KEY = (1, "#aeaeae", "#989898")

FONT = "'Ubuntu Sans Mono', monospace"
HEAD_FONT = "'Ubuntu Sans Mono'"   # panel letter + title, same string as the Inkscape panels
HEAD_SIZE = 7.0          # letter and title share one size
HEAD_GAP = 9.2           # title starts this far right of the letter, mm
INK = "#000000"
MUTED = "#555555"

S_MAX = 120.0                      # side of the plants square, mm
K = S_MAX / math.sqrt(LIFE[0][1])  # mm per sqrt(GT)
side = lambda v: K * math.sqrt(v)

X0, Y0 = 12.0, 32.0                # top-left of the plants square
BASE = Y0 + S_MAX                  # shared baseline (bottom edge)
SW = 0.5                           # stroke width

parts = []
add = parts.append


def fmt(v):
    return f"{v:.1f}".rstrip("0").rstrip(".") + " Gt C"


def label(x, y, name, value, anchor="start", size=4.2, bold_name=True):
    w = "bold" if bold_name else "normal"
    add(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-family="{FONT}" fill="{INK}">'
        f'<tspan font-size="{size}" font-weight="{w}">{name}</tspan>'
        f'<tspan font-size="{size*0.82:.2f}" fill="{MUTED}" dx="1.2">{value}</tspan></text>')


# panel letter
add(f'<text x="{X0:.2f}" y="{Y0-4.6:.2f}" font-family="{HEAD_FONT}" font-size="{HEAD_SIZE}" '
    f'font-weight="bold" fill="{INK}">A.</text>')
add(f'<text x="{X0+HEAD_GAP:.2f}" y="{Y0-4.6:.2f}" font-family="{HEAD_FONT}" font-size="{HEAD_SIZE}" '
    f'font-weight="bold" fill="{INK}">The Scale of Global Plastics Production</text>')

# life squares, largest first, anchored bottom-left
for name, v, fill, stroke in LIFE:
    s = side(v)
    add(f'<rect x="{X0:.3f}" y="{BASE-s:.3f}" width="{s:.3f}" height="{s:.3f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{SW}"/>')

# plastics square, bottom-right corner of plants
sp = side(PLASTIC[1])
px = X0 + S_MAX - sp
add(f'<rect x="{px:.3f}" y="{BASE-sp:.3f}" width="{sp:.3f}" height="{sp:.3f}" '
    f'fill="{PLASTIC[2]}" stroke="{PLASTIC[3]}" stroke-width="{SW*1.6}"/>')

# labels: plants and bacteria inside their own top-left corner
label(X0 + 3, Y0 + 7, "Plants", fmt(450), size=5.2)
sb = side(70)
label(X0 + 3, BASE - sb + 6.2, "Bacteria", fmt(70), size=4.6)

# small life squares: leader from each top-right corner to a label column
small = LIFE[2:]
col_x = X0 + side(12) + 8.5          # label column, inside bacteria
label_y = [BASE - side(12) - 1.0, BASE - 14.0, BASE - 8.0, BASE - 2.0]
for (name, v, fill, stroke), ly in zip(small, label_y):
    s = side(v)
    cx, cy = X0 + s, BASE - s
    elbow = X0 + side(12) + 2.0
    add(f'<polyline points="{cx:.2f},{cy:.2f} {elbow:.2f},{cy:.2f} '
        f'{col_x-1.2:.2f},{ly-1.2:.2f}" fill="none" stroke="{MUTED}" '
        f'stroke-width="0.25"/>')
    add(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="0.45" fill="{MUTED}"/>')
    label(col_x, ly, name, fmt(v), size=3.4)

# plastics label above its square, right aligned
label(X0 + S_MAX - 1.0, BASE - sp - 3.0, PLASTIC[0], fmt(PLASTIC[1]),
      anchor="end", size=4.6)

# key: 1 GT reference square below the chart
kv, kfill, kstroke = KEY
ks = side(kv)
ky = BASE + 5.0
add(f'<rect x="{X0:.3f}" y="{ky:.3f}" width="{ks:.3f}" height="{ks:.3f}" '
    f'fill="{kfill}" stroke="{kstroke}" stroke-width="{SW}"/>')
add(f'<text x="{X0+ks+2:.2f}" y="{ky+ks-1.3:.2f}" font-family="{FONT}" '
    f'font-size="3.4" fill="{INK}">= 1 billion tonnes of carbon (area)</text>')

W, H = max(X0 + S_MAX + 12, 160.0), ky + ks + 10   # 160 mm leaves room for the panel title
svg = (f'<?xml version="1.0" encoding="UTF-8"?>\n'
       f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}mm" height="{H:.1f}mm" '
       f'viewBox="0 0 {W:.1f} {H:.1f}">\n'
       f'<rect width="100%" height="100%" fill="#ffffff"/>\n'
       + "\n".join(parts) + "\n</svg>\n")
OUT.write_text(svg)
print(OUT, f"{W:.1f} x {H:.1f} mm")
for name, v, *_ in LIFE + [PLASTIC]:
    print(f"  {name:18s} {v:7.3f} Gt C  side {side(v):6.2f} mm")
