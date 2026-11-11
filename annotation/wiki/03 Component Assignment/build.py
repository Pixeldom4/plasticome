"""Inline the page data into template.html -> component-assignment.html.

  /*DATA*/  data.json from extract_data.py
"""
import pathlib

here = pathlib.Path(__file__).parent
html = (here / "template.html").read_text().replace("/*DATA*/null", (here / "data.json").read_text())
(here / "component-assignment.html").write_text(html)
print(here / "component-assignment.html")
