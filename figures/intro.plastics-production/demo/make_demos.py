"""Two design demos for panel 1A, using the notebook numbers (Gt C, 2026).

demo-1.dot-plot      every group on one log axis; exact, plants no longer dominates
demo-2.focus-circles plastic placed right next to animals; plants and bacteria
                     left out and noted, so the comparison that matters is large

Both: carbon units stated, year and sources in a subtitle, labels in dark ink
(colour only on the marks), title that says what the panel shows.

Run from anywhere: python demo/make_demos.py
"""
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import plastics_numbers as N  # noqa: E402

HEAD_FONT = "'Ubuntu Sans Mono'"            # same string as the other panels
BODY_FONT = "'Ubuntu Sans', Helvetica, Arial, sans-serif"
INK, MUTED, GRID = "#111111", "#5b5b5b", "#d9d9d9"
TITLE = "Plastic Waste Holds More Carbon Than All Animals"
SUB = ("Carbon mass in 2026 (Gt C). Plastics: notebook model from Geyer 2017 and "
       "Pottinger 2024. Life: Bar-On et al. 2018.")

COLORS = {  # fill, stroke (from the original panel)
    "Plants": ("#7de8a1", "#24ca5b"), "Bacteria": ("#c59bed", "#a15ee3"),
    "Fungi": ("#f0c5a6", "#e59f6c"), "Archaea": ("#ecde97", "#ddc23c"),
    "Protists": ("#7cd9de", "#33c0c8"), "Animals": ("#f2b3b4", "#e1595b"),
}
P_FILL, P_STROKE, P_DARK = "#7db3e8", "#1b5998", "#123f6d"
RATIO = N.REMAINING_C / N.LIFE_C["Animals"]


def r1(v):
    return f"{v:.1f}".rstrip("0").rstrip(".")


def svg(w, h, body):
    return (f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}mm" height="{h}mm" '
            f'viewBox="0 0 {w} {h}">\n<rect width="100%" height="100%" fill="#ffffff"/>\n'
            + "\n".join(body) + "\n</svg>\n")


def header(body, x=10.0, y=12.0):
    body.append(f'<text x="{x}" y="{y}" font-family="{HEAD_FONT}" font-size="7" '
                f'font-weight="bold" fill="{INK}">A.</text>')
    body.append(f'<text x="{x+9.2}" y="{y}" font-family="{HEAD_FONT}" font-size="7" '
                f'font-weight="bold" fill="{INK}">{TITLE}</text>')
    body.append(f'<text x="{x+9.2}" y="{y+5.5}" font-family="{BODY_FONT}" font-size="3.1" '
                f'fill="{MUTED}">{SUB}</text>')


def text(body, x, y, s, size=3.4, weight="normal", anchor="start", fill=INK):
    body.append(f'<text x="{x:.2f}" y="{y:.2f}" font-family="{BODY_FONT}" font-size="{size}" '
                f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{s}</text>')


# ---------------------------------------------------------------- demo 1
def dot_plot():
    W, H = 205, 112
    b = []
    header(b)
    x0, x1 = 52.0, 168.0                       # 1 to 1000 Gt C
    X = lambda v: x0 + math.log10(v) / 3 * (x1 - x0)
    rows = [("Plants", 450), ("Bacteria", 70), ("Fungi", 12), ("Archaea", 7),
            ("Plastic", N.REMAINING_C), ("Protists", 4), ("Animals", 2)]
    y_top, step = 30.0, 9.5
    y_axis = y_top + step * (len(rows) - 1) + 7
    # grid and axis
    for v, lab in [(1, "1"), (10, "10"), (100, "100"), (1000, "1,000")]:
        b.append(f'<line x1="{X(v):.2f}" y1="{y_top-5}" x2="{X(v):.2f}" y2="{y_axis}" '
                 f'stroke="{GRID}" stroke-width="0.3"/>')
        text(b, X(v), y_axis + 4.2, lab, size=3.0, anchor="middle", fill=MUTED)
    for v in [2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90,
              200, 300, 400, 500, 600, 700, 800, 900]:
        b.append(f'<line x1="{X(v):.2f}" y1="{y_axis-1.2}" x2="{X(v):.2f}" y2="{y_axis}" '
                 f'stroke="{MUTED}" stroke-width="0.2"/>')
    b.append(f'<line x1="{x0}" y1="{y_axis}" x2="{x1}" y2="{y_axis}" stroke="{MUTED}" stroke-width="0.3"/>')
    text(b, (x0 + x1) / 2, y_axis + 9.0, "Carbon mass, Gt C (log scale)", size=3.1,
         anchor="middle", fill=MUTED)

    for i, (name, v) in enumerate(rows):
        y = y_top + i * step
        if name == "Plastic":
            b.append(f'<rect x="{x0-46}" y="{y-4.4}" width="{x1-x0+52}" height="8.8" rx="1.2" '
                     f'fill="#eaf2fb"/>')
            # produced (hollow) to remaining (filled)
            b.append(f'<line x1="{X(N.REMAINING_C):.2f}" y1="{y}" x2="{X(N.PRODUCED_C):.2f}" y2="{y}" '
                     f'stroke="{P_STROKE}" stroke-width="0.5"/>')
            b.append(f'<circle cx="{X(N.PRODUCED_C):.2f}" cy="{y}" r="2.2" fill="#ffffff" '
                     f'stroke="{P_STROKE}" stroke-width="0.5" stroke-dasharray="0.9 0.6"/>')
            b.append(f'<circle cx="{X(N.REMAINING_C):.2f}" cy="{y}" r="2.2" fill="{P_FILL}" '
                     f'stroke="{P_STROKE}" stroke-width="0.6"/>')
            text(b, x0 - 4, y + 1.2, "Plastic still on Earth", weight="bold", anchor="end")
            text(b, X(N.PRODUCED_C) + 3.4, y - 0.4,
                 f'<tspan font-weight="bold">{r1(N.REMAINING_C)}</tspan> remaining', size=3.1)
            text(b, X(N.PRODUCED_C) + 3.4, y + 3.0,
                 f'{r1(N.PRODUCED_C)} ever produced (dashed)', size=2.7, fill=MUTED)
            y_pl = y
        else:
            fill, stroke = COLORS[name]
            b.append(f'<circle cx="{X(v):.2f}" cy="{y}" r="2.2" fill="{fill}" stroke="{stroke}" '
                     f'stroke-width="0.5"/>')
            text(b, x0 - 4, y + 1.2, name, weight="bold", anchor="end")
            text(b, X(v) + 3.4, y + 1.1, r1(v), size=3.1, fill=MUTED)
        if name == "Animals":
            y_an = y
    # ratio bracket on the left, between plastic and animals
    bx = x0 - 43
    b.append(f'<path d="M{bx+1.6},{y_pl} H{bx} V{y_an} H{bx+1.6}" fill="none" '
             f'stroke="{P_STROKE}" stroke-width="0.4"/>')
    b.append(f'<text transform="translate({bx-1.6:.2f},{(y_pl+y_an)/2:.2f}) rotate(-90)" '
             f'font-family="{BODY_FONT}" font-size="3.1" font-weight="bold" text-anchor="middle" '
             f'fill="{P_STROKE}">{RATIO:.1f}× animals</text>')
    return svg(W, H, b)


# ---------------------------------------------------------------- demo 2
def focus_circles():
    W, H = 205, 118
    b = []
    header(b)
    K = 7.0 / math.sqrt(N.LIFE_C["Animals"])      # animals radius 7 mm
    rad = lambda v: K * math.sqrt(v)
    base = 78.0
    groups = ["Fungi", "Archaea", "Protists", "Animals"]
    x = 12.0
    centres = {}
    for g in groups:
        r = rad(N.LIFE_C[g])
        x += r
        fill, stroke = COLORS[g]
        b.append(f'<circle cx="{x:.2f}" cy="{base-r:.2f}" r="{r:.2f}" fill="{fill}" '
                 f'stroke="{stroke}" stroke-width="0.5"/>')
        text(b, x, base + 5.0, g, weight="bold", anchor="middle")
        text(b, x, base + 9.0, f"{r1(N.LIFE_C[g])} Gt C", size=3.0, anchor="middle", fill=MUTED)
        centres[g] = (x, r)
        x += r + 6.0
    # plastics: dashed ring = ever produced, filled = still on Earth
    x += 8.0
    rp, rr = rad(N.PRODUCED_C), rad(N.REMAINING_C)
    px = x + rp
    b.append(f'<circle cx="{px:.2f}" cy="{base-rp:.2f}" r="{rp:.2f}" fill="none" '
             f'stroke="{P_STROKE}" stroke-width="0.5" stroke-dasharray="1.2 0.8"/>')
    b.append(f'<circle cx="{px:.2f}" cy="{base-rr:.2f}" r="{rr:.2f}" fill="{P_FILL}" '
             f'stroke="{P_STROKE}" stroke-width="0.8"/>')
    text(b, px, base + 5.0, "Plastic still on Earth", weight="bold", anchor="middle")
    text(b, px, base + 9.0, f"{r1(N.REMAINING_C)} Gt C", size=3.0, anchor="middle", fill=MUTED)
    text(b, px, base - 2 * rp - 2.2, f"dashed: all ever produced, {r1(N.PRODUCED_C)} Gt C",
         size=2.7, anchor="middle", fill=MUTED)
    # ratio arrow from animals to plastic
    ax, ar = centres["Animals"]
    ya = base - 2 * rp - 9.5
    b.append(f'<path d="M{ax:.2f},{base-2*ar-1.5:.2f} V{ya:.2f} H{px:.2f} V{base-2*rp-5.5:.2f}" '
             f'fill="none" stroke="{P_STROKE}" stroke-width="0.4"/>')
    b.append(f'<path d="M{px-1.1:.2f},{base-2*rp-7.0:.2f} L{px:.2f},{base-2*rp-5.2:.2f} '
             f'L{px+1.1:.2f},{base-2*rp-7.0:.2f}" fill="none" stroke="{P_STROKE}" stroke-width="0.4"/>')
    text(b, (ax + px) / 2, ya - 1.6, f"{RATIO:.1f}× the carbon in all animals",
         size=3.3, weight="bold", anchor="middle", fill=P_STROKE)
    # key and context note
    kr = rad(1)
    b.append(f'<circle cx="{12+kr:.2f}" cy="{base+20:.2f}" r="{kr:.2f}" fill="#aeaeae" '
             f'stroke="#989898" stroke-width="0.4"/>')
    text(b, 12 + 2 * kr + 2, base + 21.0, "= 1 Gt of carbon (area)", size=3.0, fill=MUTED)
    text(b, 12, base + 33.0,
         f"Not shown: plants ({r1(N.LIFE_C['Plants'])} Gt C) and bacteria "
         f"({r1(N.LIFE_C['Bacteria'])} Gt C), whose circles would be "
         f"{math.sqrt(N.LIFE_C['Plants']/N.LIFE_C['Animals']):.0f} and "
         f"{math.sqrt(N.LIFE_C['Bacteria']/N.LIFE_C['Animals']):.0f} times as wide as the animals circle.",
         size=2.8, fill=MUTED)
    return svg(W, H, b)


if __name__ == "__main__":
    for name, fn in [("demo-1.dot-plot", dot_plot), ("demo-2.focus-circles", focus_circles)]:
        out = HERE / f"{name}.svg"
        out.write_text(fn())
        print(out)
