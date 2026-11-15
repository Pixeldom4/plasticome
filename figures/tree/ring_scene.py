"""The ring tree as plain data, for the HTML editor (step 18).

18_ring_tree.py draws with matplotlib and, at the same call sites, records
every element here in figure inches (y up), the same coordinates matplotlib
uses. write_html() embeds the record, the structure renders and the editor
template (templates/ring_editor.html) into one file that opens from disk.

What the editor can change is written back as a layout file (Export layout),
which 18_ring_tree.py --layout reads, so a structure moved or a label
reworded in the browser stays that way in the next PDF:

  {"structures": {fold_key: {"x": in, "y": in, "size": in, "view": "full"|"zoom"}},
   "text": {text id: "new text"},
   "legends": {legend id: {"dx": in, "dy": in}},
   "hidden": [layer, ...]}
"""

import base64
import io
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "templates", "ring_editor.html")


def css(c):
    """Any matplotlib colour ("#hex", "C0", a name, an RGB(A) tuple) as CSS."""
    from matplotlib.colors import to_hex, to_rgba
    r, g, b, a = to_rgba(c)
    return to_hex((r, g, b)) if a >= 1 else f"rgba({round(r * 255)},{round(g * 255)},{round(b * 255)},{a:.3f})"


class Scene:
    def __init__(self):
        self.items = []          # paint order
        self.structures = []
        self.legends = []
        self.box = None

    def _add(self, layer, group, kind, **kw):
        self.items.append({"layer": layer, "group": group, "kind": kind, **kw})

    def wedge(self, layer, group, r0, w, a0, a1, fill, hatch=False, title=None):
        """Annular sector between radii r0 and r0 + w, angles in radians."""
        self._add(layer, group, "wedge", r0=r0, w=w, a0=min(a0, a1), a1=max(a0, a1),
                  fill=css(fill), hatch=hatch, title=title)

    def rect(self, layer, group, x0, y0, x1, y1, fill, title=None):
        """Axis aligned box in figure inches (y up), for the vertical layout."""
        self._add(layer, group, "rect", x0=x0, y0=y0, x1=x1, y1=y1, fill=css(fill), title=title)

    def poly(self, layer, group, pts, fill, title=None):
        """Filled closed path in figure inches (y up). Each point is (x, y) for a
        straight segment, or (c1x, c1y, c2x, c2y, x, y) for a cubic curve."""
        self._add(layer, group, "poly", pts=[[round(v, 4) for v in p] for p in pts],
                  fill=css(fill), title=title)

    def line(self, layer, group, pts, color, lw, dash=None):
        self._add(layer, group, "line", pts=[[round(x, 4), round(y, 4)] for x, y in pts],
                  color=css(color), lw=lw, dash=dash)

    def text(self, layer, group, x, y, s, size, color, ha="center", va="center",
             rotation=0.0, weight="normal", box=False, tid=None):
        self._add(layer, group, "text", x=x, y=y, s=s, size=size, color=css(color),
                  ha=ha, va=va, rotation=rotation, weight=weight, box=box, id=tid)

    def marker(self, layer, group, x, y, shape, color, s, title=None, hollow=False, lw=0.3):
        """s is matplotlib's marker area, in points squared."""
        self._add(layer, group, "marker", x=x, y=y, shape=shape, color=css(color), s=s,
                  title=title, hollow=hollow, lw=lw)

    def structure(self, key, x, y, size, route, name, lines, color, has_image):
        """route: {"angle": radians of the slice, "r0": where the leader leaves
        the rings, "lane": the radius it runs along}."""
        self.structures.append({"key": key, "x": x, "y": y, "size": size, "route": route,
                                "name": name, "lines": lines,
                                "color": css(color), "image": has_image})

    def legend(self, lid, x, y, title, entries, ncol=1, scale=1.0):
        """entries: {"label", "kind": square|line|marker, "color", "dash",
        "shape", "size", "hollow"}."""
        self.legends.append({"id": lid, "x": x, "y": y, "title": title,
                             "entries": entries, "ncol": ncol, "scale": scale})

    def write_html(self, path, title, box, images, scale_bar):
        data = {"title": title, "box": box, "items": self.items,
                "structures": self.structures, "legends": self.legends,
                "images": images, "scaleBar": scale_bar}
        with open(TEMPLATE) as fh:
            html = fh.read()
        html = html.replace("__TITLE__", title).replace(
            "/*__SCENE__*/null", json.dumps(data, separators=(",", ":")))
        with open(path, "w") as fh:
            fh.write(html)


def png_data_uri(path, px=1000):
    """A render, downsized for embedding."""
    from PIL import Image
    im = Image.open(path)
    im.thumbnail((px, px), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "WEBP", quality=90, method=6)
    return "data:image/webp;base64," + base64.b64encode(buf.getvalue()).decode()
