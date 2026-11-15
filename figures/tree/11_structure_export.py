#!/usr/bin/env python3
"""
Step 11: the fold representatives as a folder to hand to a folding server.

Step 10 picks one enzyme per non-singleton fold. This lays those ten out as a
self-contained folder: the sequences as FASTA, the tables behind them, and a
note of which ones already have a solved structure and which have to be
folded.

Three FASTA files, so the server can be pointed at whichever is wanted:

  fold_representatives.faa   all ten, to fold every representative with one
                             method and get models that are comparable
  folds_3plus.faa            the eight folds that hold three or more distinct
                             sequences, dropping the two that are a pair
  to_fold.faa                only the two with no structure at all, PLAase I
                             and LMC021

Eight of the ten have an experimental structure already. Their PDB entries,
with download URLs, are in tables/existing_structures.tsv, so nothing needs
refolding unless uniform models are wanted.

The centroid set is the primary one here, since the targets are centroids. The
all-sequence set differs in one slice: Savinase is not a 90 percent centroid,
so the centroid run takes Subtilisin BL for the subtilisin fold. Both are
included.

Makes no network calls and folds nothing. Writes
../structures/fold-representatives/.
"""

import csv
import os
import shutil

from common import CSV, OUT, OUT_ALL, read_fasta, read_tsv

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.path.join(HERE, "..", "structures", "fold-representatives")
PDB_URL = "https://files.rcsb.org/download/{acc}.pdb"
MIN_SEQS_SUBSET = 3


def wrap(fh, header, seq, width=60):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def main():
    if os.environ.get("PLASTICOME_TREE_SET") != "centroids":
        raise SystemExit("run with PLASTICOME_TREE_SET=centroids")

    reps = read_tsv(os.path.join(OUT, "fold_representatives.tsv"))
    seqs = read_fasta(os.path.join(OUT, "fold_representatives.faa"))
    by_slice = {h.split("|")[0]: (h, s) for h, s in seqs.items()}

    for sub in ("fasta", "tables"):
        os.makedirs(os.path.join(BUNDLE, sub), exist_ok=True)

    def fasta(name, rows):
        with open(os.path.join(BUNDLE, "fasta", name), "w") as fh:
            for r in rows:
                h, s = by_slice[r["slice"]]
                wrap(fh, h, s)
        return rows

    all_ten = fasta("fold_representatives.faa", reps)
    subset = fasta("folds_3plus.faa",
                   [r for r in reps if int(r["sequences"]) >= MIN_SEQS_SUBSET])
    gaps = fasta("to_fold.faa",
                 [r for r in reps if r["rep_has_structure"].strip().lower() != "yes"])

    # --- tables -------------------------------------------------------------
    for name, src in (
            ("fold_representatives.tsv", os.path.join(OUT, "fold_representatives.tsv")),
            ("fold_representatives_all_sequences.tsv",
             os.path.join(OUT_ALL, "fold_representatives.tsv")),
            ("fold_table.tsv", os.path.join(OUT, "fold_table.tsv")),
            ("cath_folds.tsv", os.path.join(OUT, "cath_folds.tsv")),
            ("cath_domains.tsv", os.path.join(OUT, "cath_domains.tsv")),
            ("tree_tips.tsv", os.path.join(OUT, "tree_tips.tsv"))):
        shutil.copyfile(src, os.path.join(BUNDLE, "tables", name))
    shutil.copyfile(os.path.join(OUT_ALL, "fold_representatives.faa"),
                    os.path.join(BUNDLE, "fasta", "fold_representatives_all_sequences.faa"))
    shutil.copyfile(CSV, os.path.join(BUNDLE, "tables",
                                      "plasticome-v2-union-annotated.csv"))

    with open(os.path.join(BUNDLE, "tables", "existing_structures.tsv"),
              "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["slice", "representative", "pdb_accession", "uniprot", "url"])
        for r in reps:
            if r["rep_has_structure"].strip().lower() == "yes":
                w.writerow([r["slice"], r["representative"], r["rep_pdb"],
                            r["rep_uniprot"], PDB_URL.format(acc=r["rep_pdb"])])

    # --- the folder's own note ---------------------------------------------
    rows = "\n".join(
        f"| {r['slice']} | {r['fold']} | {r['catalytic_cath']} | "
        f"{r['representative']} | {r['rep_accession']} | {r['rep_length_aa']} | "
        f"{r['rep_pdb'] if r['rep_has_structure'].strip().lower() == 'yes' else 'to fold'} | "
        f"{r['sequences']} |" for r in reps)
    with open(os.path.join(BUNDLE, "README.md"), "w") as fh:
        fh.write(NOTE.format(
            n_all=len(all_ten), n_subset=len(subset), n_gaps=len(gaps),
            n_struct=len(all_ten) - len(gaps), rows=rows,
            subset_slices=", ".join(r["slice"] for r in subset),
            gap_names=" and ".join(r["representative"] for r in gaps),
            min_len=min(int(r["rep_length_aa"]) for r in reps),
            max_len=max(int(r["rep_length_aa"]) for r in reps)))

    print(f"bundle: {os.path.normpath(BUNDLE)}")
    print(f"  fasta/fold_representatives.faa  {len(all_ten)} seqs, all non-singleton folds")
    print(f"  fasta/folds_3plus.faa           {len(subset)} seqs ({', '.join(r['slice'] for r in subset)})")
    print(f"  fasta/to_fold.faa               {len(gaps)} seqs, no structure "
          f"({', '.join(r['representative'] for r in gaps)})")
    print(f"  {len(all_ten) - len(gaps)} of {len(all_ten)} already have a PDB entry")


NOTE = """# Plasticome fold representatives

One representative enzyme for each fold of the plasticome tree figure, as
sequences to fold and the tables they came from. Built by
`figures/tree/11_structure_export.py` in the plasticome repo.

A fold is a slice of the pie in that figure: the SSN components whose
catalytic CATH superfamily is the same and whose homology was confirmed on
structure. There are 16, of which {n_all} hold more than one distinct
sequence. The other six are a single sequence that already stands for itself,
so they are not here.

The representative of each fold was picked on, in order: whether it is a
curated seed enzyme, whether it has a solved structure, the size of its SSN
component, then medoid centrality within that component. The `chosen_for`
column of `tables/fold_representatives.tsv` records which of those the winner
actually stood on.

| Slice | Fold | Catalytic CATH | Representative | Accession | Length | PDB | Sequences in fold |
|---|---|---|---|---|---|---|---|
{rows}

{n_struct} of the {n_all} already have an experimental structure. Only
{gap_names} do not.

## What to fold

| File | Sequences | Use |
|---|---|---|
| `fasta/fold_representatives.faa` | {n_all} | every representative, for one uniform set of models |
| `fasta/folds_3plus.faa` | {n_subset} | only the folds with three or more distinct sequences ({subset_slices}) |
| `fasta/to_fold.faa` | {n_gaps} | only the ones with no structure at all |

Lengths run {min_len} to {max_len} residues, all single chains, so any of
AlphaFold2, ColabFold, ESMFold or Boltz handles them. FASTA headers are
`slice|representative|accession|catalytic CATH|fold`, for example

```
>g00|TfCut_2|CBY05530.1|3.40.50.1820|Alpha/beta hydrolase
```

The slice is the stable key back into every table here, so name the output
models after it.

`tables/existing_structures.tsv` lists the {n_struct} solved entries with
RCSB download URLs, in case the existing structures are wanted beside the
predictions rather than refolded.

## Structures

`structures/` holds the representatives drawn as cartoons, one SVG and one PNG
each, every structure tinted with its fold's own colour from the pie. The
source of each is in `structures/sources.tsv`, with the chain used, its
identity to the representative, and how the secondary structure was assigned.
Eight come from RCSB, LMC021 from AlphaFold DB and PLAase I from ESMFold, so a
prediction is never mistaken for a deposition. `structures/pdb/` keeps the
coordinates each image was drawn from, and `folds_contact_sheet.{{svg,png}}` is
all ten on one sheet.

Predicted models have their dangling termini trimmed, a propeptide or signal
peptide that hangs off the folded core, so the picture is of the domain. The
residue counts in `sources.tsv` are after trimming.

## Tables

| File | Contents |
|---|---|
| `fold_representatives.tsv` | the representatives with their fold's categories: CATH class and architecture, mechanism, counts by substrate and origin, accessory domains, homology verdict, and the representative's own metadata |
| `fold_representatives_all_sequences.tsv` | the same over every union record rather than the 90 percent centroids |
| `fold_table.tsv` | fold, CATH accession, components, sequences, records, homology verdict |
| `cath_folds.tsv` | one row per fold, the fuller version from step 9 |
| `cath_domains.tsv` | one row per CATH domain per fold it occurs in, catalytic and accessory |
| `tree_tips.tsv` | one row per distinct sequence: fold, component, substrate, origin |
| `plasticome-v2-union-annotated.csv` | the source table everything is derived from |

The centroid set is primary here, since the targets are centroids. The
all-sequence set differs in one slice: Savinase is not a 90 percent centroid,
so the centroid run takes Subtilisin BL for the subtilisin fold. The two
`fold_representatives*` files let you compare them.

Note also that fold group numbering is by fold size, so the slice ids are not
interchangeable between the two sets: the subtilisin and beta-lactamase folds
swap places.
"""


if __name__ == "__main__":
    main()
