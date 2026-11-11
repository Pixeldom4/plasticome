"""Panel 1A as a treemap of carbon mass (Gt C): life vs plastics, with 2050.

Reproduces the layout mocked up on 2026-11-05. Every box area is
proportional to carbon mass (K mm^2 per Gt C), including the 1 Gt key and
the 2050 projection.

Life: plants on the left; a right column with fungi, then a row of archaea
and protists (the empty slot above protists is where animals would sit),
then bacteria. Animals are pulled out of the treemap and set beside the
plastics, joined to the "Surface" block by a wedge, so the two can be
compared directly.

Plastics (2026, plastics_numbers.py): Surface = landfill, environment or in
use; Recycled; Incinerated. 2050 = cumulative plastic produced by 2050
(Pottinger 2024 projection, same carbon fraction).

Run: python make_treemap.py  ->  1a.plastics-production.svg (the main 1A panel)
"""
import csv
import math
from pathlib import Path

import plastics_numbers as N

HERE = Path(__file__).resolve().parent
OUT = HERE / "1a.plastics-production.svg"

_rows = {int(r["Year"]): r for r in csv.DictReader(open(N.CSV))}
PROJ_YEAR = 2050
PROJ_C = float(_rows[PROJ_YEAR]["Total Produced"]) * N.CARBON_FRACTION

L = N.LIFE_C
FONT = "Arial, Helvetica, sans-serif"
G = 0.8                 # gap between boxes, mm
RX = 0.9                # corner radius, mm

COL = {  # fill, label ink
    "Plants":   ("#7de8a1", "#1f3a29"),
    "Bacteria": ("#c59bed", "#2e2338"),
    "Fungi":    ("#f0c5a6", "#3a2c22"),
    "Archaea":  ("#ecde97", "#c9a227"),
    "Protists": ("#7cd9de", "#2fa3ad"),
    "Animals":  ("#f2b3b4", "#e04848"),
}
BLUE, BLUE_LT, BLUE_DK, BLUE_PALE = "#2a5fa0", "#7db3e8", "#1b5998", "#d2e3f5"
SURFACE_INK = "#8ab4e0"

# ---------------------------------------------------------------- geometry
X0, Y0, H = 5.4, 12.0, 142.6                    # life treemap origin and height
TOTAL = sum(L.values())                          # 545 Gt C
K = 36.9                                         # mm^2 per Gt C
wP = L["Plants"] * K / H
wR = (TOTAL - L["Plants"]) * K / H
xR = X0 + wP
hF = L["Fungi"] * K / wR
hRow = (L["Archaea"] + L["Protists"] + L["Animals"]) * K / wR
wA = L["Archaea"] * K / hRow
wSub = wR - wA
hSlot = L["Animals"] * K / wSub                  # empty slot = animals' place
hPr = L["Protists"] * K / wSub
hB = L["Bacteria"] * K / wR

boxes = {  # outer cell (x, y, w, h) before the gap inset
    "Plants":   (X0, Y0, wP, H),
    "Fungi":    (xR, Y0, wR, hF),
    "Archaea":  (xR, Y0 + hF, wA, hRow),
    "Protists": (xR + wA, Y0 + hF + hSlot, wSub, hPr),
    "Bacteria": (xR, Y0 + hF + hRow, wR, hB),
}

parts = []
add = parts.append


def rect(x, y, w, h, fill, inset=True, extra=""):
    if inset:
        x, y, w, h = x + G / 2, y + G / 2, w - G, h - G
    add(f'<rect x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" rx="{RX}" '
        f'fill="{fill}"{extra}/>')
    return x, y, w, h


def text(x, y, s, size, fill, anchor="middle", weight="normal", rotate=False):
    tr = f' transform="rotate(-90 {x:.2f} {y:.2f})"' if rotate else ""
    add(f'<text x="{x:.2f}" y="{y:.2f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}"{tr}>{s}</text>')


def gt(v):
    return f"{v:.1f}".rstrip("0").rstrip(".") + " Gt"


# no title: the panel letter and title go in the composite / caption

# ---------------------------------------------------------------- life
drawn = {n: rect(*b, COL[n][0]) for n, b in boxes.items()}

x, y, w, h = drawn["Plants"]
text(x + w / 2, y + h / 2 - 1.0, "Plants", 10.4, COL["Plants"][1])
text(x + w / 2, y + h / 2 + 10.6, gt(L["Plants"]), 8.0, COL["Plants"][1])

x, y, w, h = drawn["Bacteria"]
text(x + w / 2 - 1.2, y + h / 2, "Bacteria", 8.8, COL["Bacteria"][1], rotate=True)
text(x + w / 2 + 5.4, y + h / 2, gt(L["Bacteria"]), 6.0, COL["Bacteria"][1], rotate=True)

# empty slot where animals sit in the treemap (they are drawn beside the plastics)
sx_, sy_ = xR + wA + G / 2, Y0 + hF + G / 2
add(f'<rect x="{sx_:.3f}" y="{sy_:.3f}" width="{wSub-G:.3f}" height="{hSlot-G:.3f}" rx="{RX}" '
    f'fill="#ffffff" stroke="#d6d6d6" stroke-width="0.3"/>')

x, y, w, h = drawn["Fungi"]
text(x + w / 2, y + h / 2 + 0.4, "Fungi", 5.2, COL["Fungi"][1])
text(x + w / 2, y + h / 2 + 5.6, gt(L["Fungi"]), 3.6, COL["Fungi"][1])

x, y, w, h = drawn["Archaea"]
text(x + w / 2, y + h / 2 + 1.3, gt(L["Archaea"]), 3.6, COL["Archaea"][1])
lab_end = X0 + wP - 3.2                       # labels sit inside the plants box
EM = {"Archaea": 3.7246, "Protists": 3.334}       # Arial advance widths, em
lab_x = lab_end - EM["Archaea"] * 5.2              # shared left edge, longest label ends at lab_end
ya = y + h / 2 + 1.2
text(lab_x, ya, "Archaea", 5.2, COL["Archaea"][1], anchor="start")

x, y, w, h = drawn["Protists"]
text(x + w / 2, y + h / 2 + 1.3, gt(L["Protists"]), 3.6, COL["Protists"][1])
yp = ya + 7.6
text(lab_x, yp, "Protists", 5.2, COL["Protists"][1], anchor="start")
yl = min(max(yp - 1.6, y + 1.5), y + h - 1.5)  # leader stays within the box
add(f'<line x1="{lab_x + EM['Protists'] * 5.2 + 1.4:.2f}" y1="{yl:.2f}" x2="{x+w*0.35:.2f}" y2="{yl:.2f}" '
    f'stroke="{COL["Protists"][1]}" stroke-width="0.25"/>')

# ---------------------------------------------------------------- plastics, 2026
side = math.sqrt(N.PRODUCED_C * K)
PX = 168.4
PY = sy_ + (hSlot - G) / 2 - side / 2              # plastics centred on the animal box, which sits level with the empty slot
GAPC = 4.4                                       # space between surface and right column
wS = N.REMAINING_C * K / side
wC = (N.RECYCLED_C + N.INCINERATED_C) * K / side
hRe = N.RECYCLED_C * K / wC
hIn = N.INCINERATED_C * K / wC
sx, sy, sw, sh = rect(PX, PY, wS, side, BLUE_LT, inset=False)
cx = PX + wS + GAPC
rect(cx, PY, wC, hRe - G / 2, BLUE_DK, inset=False)
rect(cx, PY + hRe + G / 2, wC, hIn - G / 2, BLUE_PALE, inset=False)

text(sx + sw / 2, sy + sh / 2 + 1.3, gt(N.REMAINING_C), 3.6, BLUE_DK)
text(sx + sw / 2, sy + sh + 4.6, "Surface", 3.6, SURFACE_INK)
lx = cx + wC + 1.4
text(lx, PY + 2.6, gt(N.RECYCLED_C), 3.6, BLUE, anchor="start")
text(lx, PY + 7.4, "Recycled", 3.6, BLUE, anchor="start")
text(lx, PY + hRe + 4.4, gt(N.INCINERATED_C), 3.6, BLUE, anchor="start")
text(lx, PY + hRe + 9.2, "Incinerated", 3.6, BLUE, anchor="start")

# bracket and heading over the plastics block
bx0, bx1, by = PX, cx + wC, PY - 2.0
bmid = (bx0 + bx1) / 2
add(f'<path d="M{bx0:.2f},{by+1.0:.2f} V{by:.2f} H{bx1:.2f} V{by+1.0:.2f} '
    f'M{bmid:.2f},{by:.2f} V{by-0.9:.2f}" fill="none" stroke="{BLUE}" stroke-width="0.25"/>')
text((bx0 + bx1) / 2, PY - 5.0, f"Plastics {N.PRODUCED_C:.1f} Gt", 5.2, BLUE)

# ---------------------------------------------------------------- animals, beside plastics
aw, ah = wSub - G, hSlot - G                     # same shape as the empty slot
ax1 = PX - 7.4
ax0, ay0 = ax1 - aw, sy_                         # same height as the empty slot
add('<defs><linearGradient id="zoom" x1="0" y1="0" x2="1" y2="0">'
    f'<stop offset="0" stop-color="{COL["Animals"][0]}" stop-opacity="0.85"/>'
    f'<stop offset="1" stop-color="{BLUE_LT}" stop-opacity="0.25"/></linearGradient></defs>')
add(f'<polygon points="{ax1:.2f},{ay0:.2f} {sx:.2f},{sy:.2f} {sx:.2f},{sy+sh:.2f} '
    f'{ax1:.2f},{ay0+ah:.2f}" fill="url(#zoom)"/>')
rect(ax0, ay0, aw, ah, COL["Animals"][0], inset=False)
text(ax0 + aw / 2, ay0 + ah / 2 + 1.3, gt(L["Animals"]), 3.6, COL["Animals"][1])
text(ax0 + aw / 2, ay0 - 2.3, "Animal", 5.2, COL["Animals"][1])
# ratio on the wedge: surface plastic carbon over animal carbon
add(f'<text x="{(ax1 + sx) / 2:.2f}" y="{ay0 + ah / 2 + 1.0:.2f}" font-family="{FONT}" '
    f'font-size="2.8" font-style="italic" text-anchor="middle" fill="#6b6b6b">'
    f'{N.REMAINING_C / L["Animals"]:.1f}×</text>')

# ---------------------------------------------------------------- key, bottom right
ks = math.sqrt(K)
kx = 156.0
ky = Y0 + H - 11.0                               # second text line ends on the treemap baseline
rect(kx, ky, ks, ks, "#676767", inset=False)
text(kx + ks + 2.2, ky + 4.6, "= 1 Gigatonne (Gt)", 5.2, "#000000", anchor="start")
text(kx + ks + 2.2 + 4.5, ky + 10.6, "of carbon", 5.2, "#000000", anchor="start")

# ---------------------------------------------------------------- 2050 projection, middle right
ps = math.sqrt(PROJ_C * K)                       # outer edge of the frame = the area
pcx, pbot, sw2 = 185.3, 92.0, 1.6
add(f'<rect x="{pcx-ps/2+sw2/2:.3f}" y="{pbot-ps+sw2/2:.3f}" width="{ps-sw2:.3f}" '
    f'height="{ps-sw2:.3f}" rx="2.0" fill="{BLUE_PALE}" stroke="{BLUE_DK}" '
    f'stroke-width="{sw2}"/>')
text(pcx, pbot - ps / 2 + 1.8, gt(PROJ_C), 5.2, BLUE)
text(pcx, pbot + 6.2, "Projected Plastics", 5.2, BLUE)
text(pcx, pbot + 12.6, str(PROJ_YEAR), 5.2, BLUE)

DY = 9.0                                         # space freed by removing the title
W, Ht = 216.0, 158.0 - DY
OUT.write_text(f'<?xml version="1.0" encoding="UTF-8"?>\n'
               f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{Ht}mm" '
               f'viewBox="0 0 {W} {Ht}">\n<rect width="100%" height="100%" fill="#ffffff"/>\n'
               + f'<g transform="translate(0,{-DY})">\n' + "\n".join(parts) + "\n</g>\n</svg>\n")
print(OUT)
print(f"plastics 2026 {N.PRODUCED_C:.2f} Gt C (side {side:.1f} mm); "
      f"{PROJ_YEAR} {PROJ_C:.2f} Gt C (side {ps:.1f} mm)")
