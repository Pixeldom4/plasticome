"""Single source for the numbers in the 1A figures.

Plastics: input/product-to-year.csv (copied from the notebook's input folder;
Geyer 2017 to 2015, Pottinger 2024 after). Carbon content uses the method in
the notebook's section 2.3: polymer shares from Geyer 2017 (PE 36.3, PP 21.0,
PVC 11.8, PET 10.2, PS 7.6; PUR and other polymers left out) times the carbon
mass fraction of each repeat unit. That gives 0.771 and reproduces the
notebook's 2024 values (9.331 Gt C produced, 6.347 Gt C remaining).

Life: Bar-On, Phillips & Milo 2018 (PNAS), Gt of carbon, as in the notebook's
"2026 Update" table (animals 2).

Run this file to print the table.
"""
import csv
from pathlib import Path

YEAR = 2026
CSV = Path(__file__).with_name("input") / "product-to-year.csv"

_C, _H, _O, _CL = 12.011, 1.008, 15.999, 35.45
_UNITS = {"PE": (2, 4, 0, 0), "PP": (3, 6, 0, 0), "PVC": (2, 3, 0, 1),
          "PET": (10, 8, 4, 0), "PS": (8, 8, 0, 0)}
_SHARE = {"PE": 36.3, "PP": 21.0, "PVC": 11.8, "PET": 10.2, "PS": 7.6}


def _cfrac(c, h, o, cl):
    return c * _C / (c * _C + h * _H + o * _O + cl * _CL)


CARBON_FRACTION = (sum(_SHARE[k] * _cfrac(*_UNITS[k]) for k in _SHARE)
                   / sum(_SHARE.values()))

_rows = {int(r["Year"]): r for r in csv.DictReader(open(CSV))}
_r = _rows[YEAR]
# plastic mass, Gt
PRODUCED = float(_r["Total Produced"])
RECYCLED = float(_r["Total Recycled"])
INCINERATED = float(_r["Total Incinerated"])
REMAINING = PRODUCED - RECYCLED - INCINERATED   # landfill, environment and in use

# carbon in that plastic, Gt C
PRODUCED_C = PRODUCED * CARBON_FRACTION
RECYCLED_C = RECYCLED * CARBON_FRACTION
INCINERATED_C = INCINERATED * CARBON_FRACTION
REMAINING_C = REMAINING * CARBON_FRACTION

# Bar-On et al. 2018, Gt C
LIFE_C = {"Plants": 450, "Bacteria": 70, "Fungi": 12, "Archaea": 7,
          "Protists": 4, "Animals": 2}

if __name__ == "__main__":
    print(f"year {YEAR}, carbon fraction {CARBON_FRACTION:.4f}")
    print(f"{'':28s}{'plastic Gt':>12s}{'carbon Gt C':>13s}")
    for name, m in [("produced", PRODUCED), ("recycled", RECYCLED),
                    ("incinerated", INCINERATED),
                    ("recycled + incinerated", RECYCLED + INCINERATED),
                    ("remaining", REMAINING)]:
        print(f"{name:28s}{m:12.3f}{m * CARBON_FRACTION:13.3f}")
    print(f"remaining carbon / animals: {REMAINING_C / LIFE_C['Animals']:.2f}x")
