# Legend — `cleaned_pazy-v260701_accession-validation_tierS.csv`

Two columns were added by the Tier S run (sequence identity to the cited paper).
Counts are over all 446 rows.

## `tier_s_label` — what Tier S could say about the row

Ordered strongest evidence → silence. The two `*_in_supp` labels are positive
confirmations that the citation graph (PubMed links) structurally cannot produce
for mining/screening papers; the rest are degrees of "no evidence", never evidence
of error.

| `tier_s_label` | Rows | Meaning | Read as |
|---|---|---|---|
| `seq_in_supp` | 88 | Row's amino-acid sequence found in a cited paper's supplementary material (exact, or length-guarded substring to absorb signal-peptide/tag trimming). | **Confirmed** — strongest |
| `acc_in_supp` | 10 | Row's accession appears in a cited paper's SI (sequence not matched — e.g. MGnify rows with no sequence available). | **Confirmed** by accession |
| `acc_in_maintext` | 7 | Accession appears only in the cited paper's main-text body, not its SI. | Weak support |
| `supp_no_match` | 2 | A cited paper's SI was fetched and parsed, but neither accession nor sequence matched. | Silence, not contradiction |
| `supp_unavailable` | 335 | No reachable/parseable SI for any cited DOI (not open-access, PDF-only SI, or publisher blocked the fetch). The row may still have a sequence — there is just nothing to check it against. | Silence, not error |
| `no_sequence` | 4 | No sequence obtainable from NCBI, UniProt/UniParc, or MGnify: 1 MGnify id (service in maintenance) + 3 `OP*` nucleotide records mis-filed in the protein column. | Cannot speak |

## `sequence_source` — where the filled `aa_sequence` came from

| `sequence_source` | Rows | Meaning |
|---|---|---|
| `ncbi` | 412 | NCBI protein, `efetch(rettype=fasta)` by accession. |
| `uniparc` | 22 | Inactive/deleted UniProt (TrEMBL) accession; sequence recovered from the UniParc archive. Mostly the Nyl* nylonases. |
| `uniprot` | 1 | Active UniProt entry whose accession is absent from NCBI (HiC, `A0A075B5G4`). |
| `ncbi_uid` | 1 | PDB entry/chain efetch won't take by accession; fetched by the uid Tier 0 resolved (PHL-7, `7NEI`). |
| *(empty)* | 10 | No sequence obtained — 7 MGnify (`MGYP…`, service down for maintenance) + 3 `OP*` nucleotide records. |

Rerunning once MGnify's maintenance ends fills the 7 `MGYP…` sequences (the cache
only re-hits the ones that failed).
