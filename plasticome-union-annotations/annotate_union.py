"""Annotate the plasticome union with fields available from both source tables.

Reads plasticome-v1-v260701-union.csv, joins each row to plasticome.v1.1 and/or
pazy-260701_pull (per its `source`), and writes
plasticome-v1-v260701-union.annotated.csv.

Merge rules for rows present in both sources:
  substrate          union of abbreviations from both
  organism           PAZy, falling back to v1.1
  superkingdom       PAZy, falling back to v1.1
  metagenomic        yes if either source marks it metagenomic
  uniprot/pdb        union of distinct valid IDs from both

PDB IDs are corrected through PDB_FIXES (typos/obsolete entries in the sources,
each verified by aligning the replacement entry to the row's sequence). Rows
with no usable pazy_id fall back to matching PAZy on accession.
"""
import csv
import hashlib
import re

import pandas as pd

WORKBOOK = "plasticome.v2.xlsx"
UNION = "plasticome-v1-v260701-union.csv"
OUT = "plasticome-v1-v260701-union.annotated.corrected.csv"
DESCRIPTION = (
    "This table contains the union set annotated with substrate, organism, superkingdom, "
    "metagenomic origin, UniProt and PDB accessions, merged from both \"pazy.v1.1\" and "
    "\"pazy-260701_pull\" (PAZy preferred where they disagree). PDB typos/obsolete IDs corrected "
    "and verified against sequence."
)
# source PDB ID -> correct ID; each replacement aligns 100% to the row's sequence
PDB_FIXES = {
    "7CWG": "7CWQ",  # ApLip-Ec: 7CWG is an unrelated phosphodiesterase
    "4EBO": "4EB0",  # LCC; letter O for zero, 4EBO does not exist
    "501L": "5O1L",  # LCP_STRK3; zero for letter O, 501L does not exist
    "6IQ9": "6M43",  # E4; 6IQ9 obsoleted, superseded by 6M43
}
NEW_COLS = ["substrate", "organism", "superkingdom", "metagenomic",
            "uniprot_accession", "pdb_accession", "has_structure"]

SUPERKINGDOM = {
    "Bacteria": "Bacteria", "Bacillati": "Bacteria", "Pseudomonadati": "Bacteria",
    "Thermotogati": "Bacteria",
    "Archaea": "Archaea", "Thermoproteati": "Archaea",
    "Eukarya": "Eukaryota", "Eukaryota": "Eukaryota", "Eukaryotic": "Eukaryota",
    "Fungi": "Eukaryota", "Metazoa": "Eukaryota", "Viridiplantae": "Eukaryota",
}
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
PDB_RE = re.compile(r"^[0-9][A-Z0-9]{3}$")
METAGENOME_RE = re.compile(r"metagenom|uncultured|environmental", re.IGNORECASE)


def blank(x):
    return pd.isna(x) or str(x).strip() in ("", "-", "NA", "nan", "None")


def md5(seq):
    return hashlib.md5(str(seq).strip().upper().encode()).hexdigest()


def strip_version(acc):
    return str(acc).strip().split(".")[0]


def split(x):
    return [] if blank(x) else [t.strip() for t in re.split(r"[;/]", str(x)) if t.strip()]


def uniprot_ids(x):
    # v1.1 sometimes holds entry names (A0A0G3BI90_9BURK); keep the accession part
    return [t.split("_")[0] for t in split(x) if UNIPROT_RE.match(t.split("_")[0])]


def pdb_ids(x):
    ids = [t.split("_")[0].upper() for t in split(x) if PDB_RE.match(t.split("_")[0].upper())]
    return [PDB_FIXES.get(i, i) for i in ids]


def join(values):
    out = list(dict.fromkeys(v for v in values if v))
    return "; ".join(out) if out else "NA"


def main():
    union = pd.read_csv(UNION, header=1, dtype=str, keep_default_na=False)
    v1 = pd.read_excel(WORKBOOK, "plasticome.v1.1", dtype=str)
    pazy = pd.read_excel(WORKBOOK, "pazy-260701_pull", dtype=str)

    # v1.1: match accession against any of its ID columns, fall back to sequence
    v1_by_acc = {}
    for i, r in v1.iterrows():
        for col in ("pazy_accession", "uniprot_accession", "pdb_accession", "BLAST_link"):
            for acc in split(r[col]):
                v1_by_acc.setdefault(strip_version(acc), i)
    v1_by_seq = {md5(s): i for i, s in v1.aa_sequence.items()}
    pazy_by_id = {strip_version(x): i for i, x in pazy.pazy_id.items()}
    pazy_by_acc = {strip_version(x): i for i, x in pazy.accession.items() if not blank(x)}

    rows, unmatched = [], []
    for _, u in union.iterrows():
        v = p = None
        if u.source in ("v1.1", "both"):
            i = v1_by_acc.get(strip_version(u.accession), v1_by_seq.get(md5(u.aa_sequence)))
            v = v1.loc[i] if i is not None else None
        if u.source in ("v260701", "both"):
            i = pazy_by_id.get(strip_version(u.pazy_id), pazy_by_acc.get(strip_version(u.accession)))
            p = pazy.loc[i] if i is not None else None
        if (u.source != "v260701" and v is None) or (u.source != "v1.1" and p is None):
            unmatched.append((u.plasticome_id, u.source, u.accession))

        organism = next((x for x in (p is not None and p.organism_name, v is not None and v.host)
                         if x and not blank(x)), "NA")
        kingdoms = [k for k in (p is not None and p.kingdom, v is not None and v.kingdom)
                    if k and not blank(k)]
        superkingdom = next((SUPERKINGDOM[k] for k in kingdoms if k in SUPERKINGDOM), "NA")
        metagenomic = ("Metagenome" in kingdoms) or bool(METAGENOME_RE.search(organism))
        substrates = sorted({s for src in (p is not None and p.substrate_abbrev,
                                           v is not None and v.substrate)
                             if src for s in split(src)})
        uniprot = [a for src in (p, v) if src is not None for a in uniprot_ids(src.uniprot_accession)]
        pdb = [a for src in (p, v) if src is not None for a in pdb_ids(src.pdb_accession)]

        row = u.to_dict()
        seq = row.pop("aa_sequence")
        row.update(substrate=join(substrates), organism=organism, superkingdom=superkingdom,
                   metagenomic="yes" if metagenomic else "no",
                   uniprot_accession=join(uniprot), pdb_accession=join(pdb),
                   has_structure="yes" if pdb else "no", aa_sequence=seq)
        rows.append(row)

    columns = [c for c in union.columns if c != "aa_sequence"] + NEW_COLS + ["aa_sequence"]
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([DESCRIPTION] + [""] * (len(columns) - 1))
        w.writerow(columns)
        for r in rows:
            w.writerow([r[c] for c in columns])

    out = pd.DataFrame(rows)
    print(f"wrote {len(out)} rows to {OUT}")
    for c in NEW_COLS:
        filled = (~out[c].isin(["NA", "no"])).sum()
        print(f"  {c:18s} {filled:4d}/{len(out)}  top: {out[c].value_counts().head(4).to_dict()}")
    print(f"unmatched to a source row: {unmatched}")


if __name__ == "__main__":
    main()
