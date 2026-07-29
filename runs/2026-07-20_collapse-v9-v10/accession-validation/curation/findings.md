# Curation findings (human review)

Findings from human review that the automated pipeline cannot act on by itself.
Per-accession verdicts that *do* change a row's label live in `verdicts.csv` and are
applied by Tier 3; this file records the rest.

## Data-completeness gaps (source sheet)

- **Missing row — amidase paper strain C.** The amidase-screening paper
  `10.1016/j.jbiosc.2021.05.003` screened seven amidases (strains A–G). Six are in the
  sheet (their organism assignments were checked against the paper's Materials & Methods
  and all match). The seventh — ***Klebsiella oxytoca* 1686, `CP003218.1`** (strain C) —
  has **no row**. It should be added to the source sheet if the plasticome is meant to
  cover the full screen.

## Accession-provenance caveats (now auto-verified by Tier P)

The screening paper reports several sources by **nucleotide** accession; the sheet
carries a **protein** accession. These were not verifiable from the paper text, so
Tier P (`tierP_provenance.py`) checked each against NCBI by fetching the cited
nucleotide record's coding sequences and confirming the protein is among them
(`verdict=verify_cds` in verdicts.csv):

| enzyme | organism | paper gives | sheet has | status |
|---|---|---|---|---|
| Amidase | *R. erythropolis* MP50 | `AY026386` (nt) | `AAK11724.1` (prot) | **auto-verified** — CDS of AY026386 |
| Amidase | *M. hydrocarbonoxydans* | `GU116480` (nt) | `ACY56506.1` (prot) | **auto-verified** — CDS of GU116480 |
| Amidase | *B. phytofirmans* DSM17436 | `ACD16728.1` (prot) | `WP_012433325.1` (RefSeq) | verified correct (same organism) |

All three amidase provenance items are now resolved (`verified`); none remain in the
review queue.
