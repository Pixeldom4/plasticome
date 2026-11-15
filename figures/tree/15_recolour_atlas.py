#!/usr/bin/env python3
"""
Step 15: recolour the fold palette to the hues the PETadex atlas uses.

The atlas (petadex.net/atlas) does not store a hex per fold. It stores a base
hue per CATH superfamily, in `frontend/src/utils/cathColors.js` as `CATH_HUE`,
and builds every colour from it: HSL(hue, 60, 38) for the superfamily and
HSL(hue, 70, 35 to 70) for the components inside it. So "track the atlas" means
take its hue assignment, not its hex values, which are tuned for light points
on a dark map and would swamp the substrate coloured tips if used as slice
backgrounds.

Each fold therefore keeps the pale background role it has in the figure and
only its hue moves:

    new tint = HSL(atlas hue for the fold's CATH, S, L)

with S and L the means of the ten tints the figure already uses, 52 and 89, so
the figure is no paler or stronger than before and hue is the only thing that
changed.

The atlas defines six of the ten. The other four CATH superfamilies are not in
the PETadex component set, and the atlas would draw them in its unknown grey,
which is no use here where four grey slices would be indistinguishable from
each other and from the "Other folds" slice. They are given new hues in the
atlas's own idiom, placed in the widest gaps left by the six fixed ones, so the
ten stay as far apart as the fixed assignments allow.

Recolouring is done on the SVG rather than by redrawing, so hand editing in
Inkscape survives. Every colour the figure derives from a fold tint is
computed, old and new, and substituted in one pass: the tint itself, the
`shade(tint, 0.955)` of the alternating component band, and the 26 depth bands
of the structure cartoons with their outlines. Substituting in several passes
would chain, since one fold's new hue can equal another's old one.

Usage: 15_recolour_atlas.py INPUT.svg [-o OUTPUT.svg]
"""

import argparse
import colorsys
import importlib
import os
import re
import subprocess

pie = importlib.import_module("08_draw_pie")        # shade()

HERE = os.path.dirname(os.path.abspath(__file__))

# The ten folds of the figure and the catalytic CATH superfamily of each.
FOLD_CATH = {
    "g00": "3.40.50.1820", "g01": "3.60.70.12", "g02": "3.90.1300.10",
    "g03": "none", "g04": "3.40.710.10", "g05": "2.40.10.10",
    "g06": "3.40.50.200", "g07": "1.10.760.10", "g08": "2.60.40.420",
    "g09": "3.40.50.1110",
}

# petadex.io/frontend/src/utils/cathColors.js, CATH_HUE
ATLAS_HUE = {
    "3.40.50.1820": 210, "3.40.710.10": 140, "3.60.70.12": 25,
    "3.90.1300.10": 0, "2.40.10.10": 270, "2.60.40.420": 50,
}

# Not in the atlas. Placed in the widest gaps left by the six above
# (50 to 140, 140 to 210, 210 to 270, 270 to 360).
EXTRA_HUE = {
    "none": 90,            # g03 oxygenase MpaB family
    "3.40.50.200": 175,    # g06 peptidase S8/S53
    "3.40.50.1110": 240,   # g09 SGNH hydrolase
    "1.10.760.10": 315,    # g07 cytochrome c like
}

OLD_TINT = {"g00": "#f3ead6", "g01": "#f5dde2", "g02": "#e4dff3",
            "g03": "#dcedd6", "g04": "#d3ece8", "g05": "#fbe2cf",
            "g06": "#d8e6f5", "g07": "#f4efc2", "g08": "#eddcee",
            "g09": "#e2e8cc"}

SAT, LIGHT = 51.9, 89.0      # means of the ten tints above
SECOND_TINT = 0.955          # 08_draw_pie.py, the alternating component band
BANDS, DEPTH, OUTLINE = 26, (0.74, 1.10), 0.52      # 12_render_structures.py


def hsl_hex(h, s, l):
    r, g, b = colorsys.hls_to_rgb(h / 360.0, l / 100.0, s / 100.0)
    return "#%02x%02x%02x" % (round(r * 255), round(g * 255), round(b * 255))


def hex_hls(h):
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return colorsys.rgb_to_hls(r, g, b)


def text_recolour(svg, new_hue, old_hue, tol=3.0, min_sat=0.15):
    """Hand coloured text, matched to a fold by hue and moved with it.

    The figure's own labels are grey. A coloured `<text>` is something added
    by hand in Inkscape, in practice a fold label tinted to match its slice.
    Those are matched to a fold on hue alone and recoloured to the fold's new
    hue, keeping the lightness and saturation they were given, so the styling
    chosen by hand survives and only the hue follows the slice.
    """
    out = {}
    for c in set(re.findall(r'<text[^>]*fill:\s*(#[0-9a-fA-F]{6})', svg)):
        h, l, sat = hex_hls(c.lower())
        if sat < min_sat:
            continue                         # a grey, one of the figure's own
        h *= 360
        g = min(old_hue, key=lambda k: min(abs(old_hue[k] - h), 360 - abs(old_hue[k] - h)))
        if min(abs(old_hue[g] - h), 360 - abs(old_hue[g] - h)) > tol:
            continue                         # not one of the fold hues
        out[c.lower()] = hsl_hex(new_hue[g], sat * 100, l * 100)
    return out


def derived(tint):
    """Every colour the figure makes from one fold tint."""
    out = [tint, pie.shade(tint, SECOND_TINT)]
    for b in range(BANDS):
        light = DEPTH[0] + ((b + 0.5) / BANDS) * (DEPTH[1] - DEPTH[0])
        out.append(pie.shade(tint, light))
        out.append(pie.shade(tint, light * OUTLINE))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("svg")
    ap.add_argument("-o", "--out")
    ap.add_argument("--no-png", action="store_true")
    ap.add_argument("--background", default="#fcfcfb",
                    help="PNG background. The SVG keeps its own; this only "
                         "stops transparent margins rendering black.")
    args = ap.parse_args()
    out = args.out or args.svg.replace(".svg", "-atlas.svg")

    new_tint = {g: hsl_hex(ATLAS_HUE.get(c, EXTRA_HUE.get(c)), SAT, LIGHT)
                for g, c in FOLD_CATH.items()}

    # A colour two folds both produce cannot be recoloured, since there is no
    # way to tell which fold an instance of it belongs to. The top depth band
    # saturates to white for every pale tint, so white is always one of these.
    # They are dropped from the map and left as they are, which is right: they
    # are the highlights, and a highlight carries no fold identity.
    mapping, ambiguous = {}, set()
    for g, old in OLD_TINT.items():
        for a, b in zip(derived(old), derived(new_tint[g])):
            if a in mapping and mapping[a] != b:
                ambiguous.add(a)
            mapping[a] = b
    for a in ambiguous:
        del mapping[a]

    svg = open(args.svg).read()
    hue_of = {g: ATLAS_HUE.get(c, EXTRA_HUE.get(c)) for g, c in FOLD_CATH.items()}
    old_hue = {g: hex_hls(t)[0] * 360 for g, t in OLD_TINT.items()}
    hand = text_recolour(svg, hue_of, old_hue)
    mapping.update(hand)
    hit = set()

    def sub(m):
        c = m.group(0).lower()
        if c in mapping:
            hit.add(c)
            return mapping[c]
        return m.group(0)

    new = re.sub(r'#[0-9a-fA-F]{6}\b', sub, svg)
    with open(out, "w") as fh:
        fh.write(new)

    print(f"{len(mapping)} colours in the map, {len(hit)} of them present and "
          f"replaced; {len(ambiguous)} left alone as shared highlights "
          f"({', '.join(sorted(ambiguous))})")
    if hand:
        print(f"{len(hand)} hand coloured text fills moved with their fold: "
              + ", ".join(f"{a}->{b}" for a, b in sorted(hand.items())))
    print(f"{'slice':6} {'CATH':15} {'hue':>4}  {'old':9} -> {'new':9}  source")
    for g in sorted(FOLD_CATH):
        c = FOLD_CATH[g]
        h = ATLAS_HUE.get(c)
        print(f"{g:6} {c:15} {h if h is not None else EXTRA_HUE[c]:>4}  "
              f"{OLD_TINT[g]:9} -> {new_tint[g]:9}  "
              f"{'atlas' if h is not None else 'new, atlas idiom'}")
    print(f"\nwrote {out}")

    if not args.no_png:
        png = out.replace(".svg", ".png")
        w = int(float(re.search(r'viewBox="0 0 ([\d.]+)', new).group(1)) * 4)
        subprocess.run(["rsvg-convert", "-w", str(w), "-b", args.background,
                        "-o", png, out], check=True)
        print(f"wrote {png}")


if __name__ == "__main__":
    main()
