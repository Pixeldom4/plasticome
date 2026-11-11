# Plasticome fold representatives

One representative enzyme for each fold of the plasticome tree figure, as
sequences to fold and the tables they came from. Built by
`figures/tree/11_structure_export.py` in the plasticome repo.

A fold is a slice of the pie in that figure: the SSN components whose
catalytic CATH superfamily is the same and whose homology was confirmed on
structure. There are 16, of which 10 hold more than one distinct
sequence. The other six are a single sequence that already stands for itself,
so they are not here.

The representative of each fold was picked on, in order: whether it is a
curated seed enzyme, whether it has a solved structure, the size of its SSN
component, then medoid centrality within that component. The `chosen_for`
column of `tables/fold_representatives.tsv` records which of those the winner
actually stood on.

| Slice | Fold | Catalytic CATH | Representative | Accession | Length | PDB | Sequences in fold |
|---|---|---|---|---|---|---|---|
| g00 | Alpha/beta hydrolase | 3.40.50.1820 | TfCut_2 | CBY05530 | 261 | 4CG1 | 309 |
| g01 | L-amino peptidase D-ALA esterase/amidase | 3.60.70.12 | NylC | WP_012476895.1 | 355 | 5XYG | 45 |
| g02 | Amidase signature (AS) domain | 3.90.1300.10 | NylA | WP_012476897.1 | 493 | 3A2P | 19 |
| g03 | Oxygenase MpaB family protein | none | LCP_STRK3 | Q3L8N0.2 | 407 | 5O1L | 10 |
| g04 | DD-peptidase/beta-lactamase superfamily | 3.40.710.10 | NylB | CAA24927.1 | 392 | 2E8I | 5 |
| g05 | Trypsin-like serine proteases | 2.40.10.10 | PLAase I | WP_037320891.1 | 245 | to fold | 5 |
| g06 | Peptidase S8/S53 domain | 3.40.50.200 | Subtilisin BL | P29599 | 269 | 1ST3 | 5 |
| g07 | Cytochrome c-like domain | 1.10.760.10 | RoxA | AAP41848.1 | 678 | 4B2N | 3 |
| g08 | Cupredoxins - blue copper proteins | 2.60.40.420 | LMC021 | AII11185.1 | 516 | to fold | 2 |
| g09 | SGNH hydrolase | 3.40.50.1110 | E4 | 7C2A_A | 191 | 6M43 | 2 |

8 of the 10 already have an experimental structure. Only
PLAase I and LMC021 do not.

## What to fold

| File | Sequences | Use |
|---|---|---|
| `fasta/fold_representatives.faa` | 10 | every representative, for one uniform set of models |
| `fasta/folds_3plus.faa` | 8 | only the folds with three or more distinct sequences (g00, g01, g02, g03, g04, g05, g06, g07) |
| `fasta/to_fold.faa` | 2 | only the ones with no structure at all |

Lengths run 191 to 678 residues, all single chains, so any of
AlphaFold2, ColabFold, ESMFold or Boltz handles them. FASTA headers are
`slice|representative|accession|catalytic CATH|fold`, for example

```
>g00|TfCut_2|CBY05530.1|3.40.50.1820|Alpha/beta hydrolase
```

The slice is the stable key back into every table here, so name the output
models after it.

`tables/existing_structures.tsv` lists the 8 solved entries with
RCSB download URLs, in case the existing structures are wanted beside the
predictions rather than refolded.

## Structures

`structures/` holds the representatives drawn as cartoons, one SVG and one PNG
each, every structure tinted with its fold's own colour from the pie. The
source of each is in `structures/sources.tsv`, with the chain used, its
identity to the representative, and how the secondary structure was assigned.
Eight come from RCSB, LMC021 from AlphaFold DB and PLAase I from ESMFold, so a
prediction is never mistaken for a deposition. `structures/pdb/` keeps the
coordinates each image was drawn from, and `folds_contact_sheet.{svg,png}` is
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
