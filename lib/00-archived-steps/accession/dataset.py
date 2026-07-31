"""Loading and writing the accession sheet."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

BASE = Path(__file__).parent
DATA = BASE / "input-data" / "accession-validation.csv"
OUT = BASE / "out"

# Columns the pipeline requires. The sheet may carry more (organism,
# ncbi_taxonomy_id, ...) which we read when present and ignore otherwise.
REQUIRED_COLUMNS = ["enzyme_name", "accession", "direct_or_indirect_activity", "pazy_id", "doi", "aa_sequence"]


@dataclass
class Row:
    index: int
    enzyme_name: str
    accession: str
    direct_or_indirect_activity: str
    pazy_id: str
    doi: str
    aa_sequence: str = ""
    organism: str = ""
    ncbi_taxonomy_id: str = ""
    flags: dict = field(default_factory=dict)


def load_rows(path: Path = DATA) -> list[Row]:
    """Read the sheet by column name.

    Reads columns by header rather than position, so extra curator columns
    (organism, ncbi_taxonomy_id) are picked up when present and absent ones
    default to empty. Fully blank rows are skipped.
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        header = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise ValueError(f"sheet is missing required columns {missing}; header={header}")

        rows = []
        for i, rec in enumerate(reader):
            cell = lambda k: (rec.get(k) or "").strip()
            if not any(cell(c) for c in REQUIRED_COLUMNS):
                continue
            rows.append(Row(
                index=i,
                enzyme_name=cell("enzyme_name"),
                accession=cell("accession"),
                direct_or_indirect_activity=cell("direct_or_indirect_activity"),
                pazy_id=cell("pazy_id"),
                doi=cell("doi"),
                aa_sequence=cell("aa_sequence"),
                organism=cell("organism"),
                ncbi_taxonomy_id=cell("ncbi_taxonomy_id"),
            ))
    return rows


def write_json(name: str, obj) -> Path:
    OUT.mkdir(exist_ok=True)
    p = OUT / name
    p.write_text(json.dumps(obj, indent=2, sort_keys=True))
    return p


def read_json(name: str):
    return json.loads((OUT / name).read_text())


def write_csv(name: str, fieldnames: list[str], records: list[dict]) -> Path:
    OUT.mkdir(exist_ok=True)
    p = OUT / name
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)
    return p
