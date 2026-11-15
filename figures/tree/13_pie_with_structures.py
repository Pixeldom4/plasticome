#!/usr/bin/env python3
"""
Step 13: the pie with a strip of representative structures under it.

Takes `figure-pie-hand-edited/figure_pie.svg`, the hand cleaned version of the figure,
and writes a new file with the step 12 structures added as a labelled strip
along the bottom. The original is copied to a backup first and is never
edited in place, so the cleaning done in Inkscape survives.

The strip rather than the rim. Nine of the ten slices are narrow and sit
inside a 100 degree arc of the pie, and their fold labels already run down
the left with leader lines, so there is no room outside the rim for ten
thumbnails without moving work that was positioned by hand. A strip along the
bottom leaves every existing element where it is. The tie back to the slice is
the colour, which is the slice's own tint, and the slice id under each
structure.

The structures go in their own Inkscape layer, `structures`, so the strip can
be moved, restyled or deleted as one object. Each structure is a nested svg
with its ids prefixed, since a nested svg shares the id namespace of the file
it is dropped into and matplotlib names its clip paths the same way every
time.

Writes figure-pie-hand-edited/figure_pie_structures.svg and .png.
"""

import os
import re
import shutil
import subprocess

from common import OUT, read_tsv

HERE = os.path.dirname(os.path.abspath(__file__))
PIE = os.path.join(HERE, "figure-pie-hand-edited", "figure_pie.svg")
BACKUP = os.path.join(HERE, "figure-pie-hand-edited", "backups",
                      "figure_pie_backup_before_structures.svg")
NEW = os.path.join(HERE, "figure-pie-hand-edited", "figure_pie_structures.svg")
SDIR = os.path.join(HERE, "..", "structures", "fold-representatives", "structures")

STRIP = 186.0          # points of canvas added at the bottom
BOX = 108.0            # the square a structure is drawn in
TOP_PAD = 14.0         # from the old canvas edge to the top of the boxes
CAPTION_CHARS = 25     # the fold name wraps to two lines at about this width
FONT = ("font-family: 'DejaVu Sans', 'Bitstream Vera Sans', 'Computer Modern "
        "Sans Serif', 'Lucida Grande', 'Verdana', 'Geneva', 'Lucid', 'Arial', "
        "'Helvetica', 'Avant Garde', sans-serif")
INK = "http://www.inkscape.org/namespaces/inkscape"


def shrink(svg, prefix):
    """Namespace the ids and drop pointless coordinate precision."""
    ids = set(re.findall(r'id="([^"]+)"', svg))
    for i in sorted(ids, key=len, reverse=True):
        svg = svg.replace(f'id="{i}"', f'id="{prefix}{i}"')
        svg = svg.replace(f'url(#{i})', f'url(#{prefix}{i})')
        svg = svg.replace(f'href="#{i}"', f'href="#{prefix}{i}"')
    return re.sub(r'(-?\d+\.\d{3,})',
                  lambda m: f"{float(m.group(1)):.2f}".rstrip("0").rstrip("."), svg)


def inner(svg):
    """The body and the viewBox of a standalone svg file."""
    vb = re.search(r'viewBox="([^"]+)"', svg).group(1)
    body = svg[svg.index(">", svg.index("<svg")) + 1:svg.rindex("</svg>")]
    body = re.sub(r"<metadata>.*?</metadata>", "", body, flags=re.S)
    return body, vb


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def wrap(text, width=CAPTION_CHARS, lines=2):
    """Greedy wrap of a fold name to at most `lines` lines."""
    out, cur = [], ""
    for word in text.split():
        trial = f"{cur} {word}".strip()
        if len(trial) <= width or not cur:
            cur = trial
        else:
            out.append(cur)
            cur = word
    out.append(cur)
    if len(out) > lines:
        out = out[:lines]
        out[-1] = out[-1][:width - 1].rstrip() + "\u2026"
    return out


def main():
    reps = read_tsv(os.path.join(OUT, "fold_representatives.tsv"))
    fill = {r["slice"]: r["fill"] for r in read_tsv(os.path.join(OUT, "cath_folds.tsv"))}
    src = {r["slice"]: r["source"] for r in read_tsv(os.path.join(SDIR, "sources.tsv"))}

    pie = open(PIE).read()
    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', pie)
    W, H = float(m.group(1)), float(m.group(2))
    if not os.path.exists(BACKUP):
        shutil.copyfile(PIE, BACKUP)
        print(f"backup: {os.path.relpath(BACKUP, HERE)}")
    else:
        print(f"backup already there, left alone: {os.path.relpath(BACKUP, HERE)}")

    files = {f.split("_")[0]: os.path.join(SDIR, "svg", f)
             for f in os.listdir(os.path.join(SDIR, "svg")) if f.endswith(".svg")}
    n = len(reps)
    step = W / n
    x0 = (step - BOX) / 2
    top = H + TOP_PAD

    parts = [f'<g id="layer-structures" {{ink}}:groupmode="layer" '
             f'{{ink}}:label="structures">'.format(ink="inkscape")]
    for i, rep in enumerate(reps):
        g = rep["slice"]
        body, vb = inner(open(files[g]).read())
        x = x0 + i * step
        parts.append(f'<g id="struct-{g}" inkscape:label="{esc(rep["representative"])}">')
        parts.append(f'<svg x="{x:.2f}" y="{top:.2f}" width="{BOX}" height="{BOX}" '
                     f'viewBox="{vb}" overflow="visible">'
                     + shrink(body, f"{g}-") + '</svg>')
        cx = x + BOX / 2
        ty = top + BOX + 13
        note = {"rcsb": rep["rep_pdb"], "afdb": "AlphaFold", "esmfold": "ESMFold"}
        caption = [(0, 8.6, "bold", "#1a1a1a", f'{g}  {rep["representative"]}')]
        for k, line in enumerate(wrap(rep["fold"])):
            caption.append((11.0 + 8.8 * k, 7.0, "normal", "#4a4a4a", line))
        caption.append((11.0 + 8.8 * len(wrap(rep["fold"])) + 1.5, 6.8, "normal",
                        "#6b6b6b", f'{rep["catalytic_cath"]}  ({note[src[g]]})'))
        for dy, size, weight, colour, text in caption:
            parts.append(
                f'<text x="{cx:.2f}" y="{ty + dy:.2f}" style="font-size: {size}px; '
                f'{FONT}; font-weight: {weight}; text-anchor: middle; fill: {colour}">'
                f'{esc(text)}</text>')
        parts.append('</g>')
    parts.append('</g>')

    # the figure's own background path stops at the old height, so the strip
    # would come out transparent. Extend it rather than lay a rect over it.
    body_src = pie[pie.index(">", pie.index("<svg")) + 1:pie.rindex("</svg>")]
    bg = re.search(r'(<path d="M 0 )([\d.]+)( +L +[\d.]+ )([\d.]+)( .*?z " style="fill: #[a-f0-9]{6}" />)',
                   body_src)
    if not bg:
        raise SystemExit("could not find the figure background path to extend")
    body_src = body_src.replace(
        bg.group(0),
        f"{bg.group(1)}{H + STRIP}{bg.group(3)}{H + STRIP}{bg.group(5)}", 1)

    head_end = pie.index(">", pie.index("<svg")) + 1
    head = pie[:head_end]
    head = re.sub(r'width="[\d.]+pt"', f'width="{W:.6f}pt"', head)
    head = re.sub(r'height="[\d.]+pt"', f'height="{H + STRIP:.6f}pt"', head)
    head = re.sub(r'viewBox="0 0 [\d.]+ [\d.]+"', f'viewBox="0 0 {W} {H + STRIP}"', head)
    if 'xmlns:inkscape' not in head:
        head = head.replace("<svg ", f'<svg xmlns:inkscape="{INK}" ', 1)

    out = head + body_src + "\n" + "\n".join(parts) + "\n</svg>\n"
    with open(NEW, "w") as fh:
        fh.write(out)

    png = NEW[:-4] + ".png"
    subprocess.run(["rsvg-convert", "-w", str(int(W * 4)), "-o", png, NEW], check=True)
    print(f"wrote {os.path.relpath(NEW, HERE)}  "
          f"{os.path.getsize(NEW) / 1e6:.1f} MB, canvas {W:.0f} x {H + STRIP:.0f} pt")
    print(f"wrote {os.path.relpath(png, HERE)}  {os.path.getsize(png) / 1e6:.1f} MB")
    print(f"{n} structures in layer 'structures', {BOX:.0f} pt boxes")


if __name__ == "__main__":
    main()
