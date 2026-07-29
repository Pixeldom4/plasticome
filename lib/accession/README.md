# Plasticome accession validation

Validates the 446 plastic-degrading-enzyme rows in `accession-validation.csv`: the
accessions are real and cited but not trusted, and the dominant failure mode is an
accession that resolves fine yet points to the *wrong* protein. See the plan
(`Accession Validation Plan.md`) for the full rationale.

## Run

```
NCBI_EMAIL=you@example.com python3 run_all.py
```

Credentials live in a gitignored `.env` (`NCBI_API_KEY`, `NCBI_EMAIL`); source it
first so the run picks them up (a non-interactive shell does not inherit your
profile exports):

```
set -a; . ./.env; set +a; python3 run_all.py
```

Every NCBI/publisher response is cached under `.cache/`, `.supp_cache/` and
`.blastcache/`, so a second run is near-instant. Delete a cache dir to re-query
from scratch. An `NCBI_ENV_VAR` (or `NCBI_API_KEY`) API key lifts the E-utilities
rate limit from 3 to 10 req/s.

Outputs land in `out/`:
- `validation_full.csv` — every row with its derived `confidence` label + `action`.
- `review_queue.csv` — only the rows a curator must look at, most-urgent first.
- `accession-validation-enriched.csv` — the input sheet with `aa_sequence` filled
  from NCBI (403/446 rows).
- `summary.json` — counts per tier.

## Tiers

| Tier | File | What it decides |
|------|------|-----------------|
| 0 | `tier0_resolve.py` | Existence/format: does the accession resolve to a live NCBI *protein*? |
| 1 | `tier1_names.py` | Name concordance: does the record's family match `enzyme_name`? (abstains on opaque codes) |
| 2 | `tier2_doi.py` | Citation graph: is the accession linked to the cited DOI in PubMed? |
| O | `tierO_organism.py` | Organism concordance: does the record's `taxid` agree with PAZy's, rank-aware? |
| S | `tierS_sequence.py` | **Sequence identity to the cited paper** — see below. |
| A | `tierA_remap.py` | **Sequence-driven re-mapping**: for stuck rows, BLAST a candidate set and rank it by organism + cited-author overlap — see below. |
| P | `tierP_provenance.py` | Protein↔nucleotide provenance: is the protein a CDS of a paper-cited nucleotide accession? |
| 3 | `tier3_report.py` | Merges the tier flags + curator verdicts into one `confidence` label. |

## Curator verdicts (`curation/verdicts.csv`)

Human review beats the heuristics without loosening any tier's logic. Each adjudicated
row is keyed by accession with a `verdict`:
- `correct` — row is right despite its flag → relabelled `verified`, dropped from the
  queue (its automated label kept in `auto_confidence`).
- `verify_cds` + an `evidence_accession` (a nucleotide accession the paper cites) →
  **Tier P** confirms the sheet's protein is a CDS of that record (by `protein_id` or by
  sequence). Confirmed → `verified`; contradicted → `provenance_mismatch` (top of queue).
- `note` — leaves the row in the queue with a caveat attached for the next reviewer.

Re-run `python3 tier3_report.py` (no network) after editing verdicts to fold them in;
`verify_cds` rows also need `python3 tierP_provenance.py` first.

## Tier O — organism concordance

For each row, compares what PAZy claims the source organism is (`organism` /
`ncbi_taxonomy_id`) against what the NCBI record itself reports (`taxid`, captured by
Tier 0). Comparison is on **taxid, walked to a common rank** (species, then genus) via
`ncbi.efetch_taxonomy` — not free-text names, which diverge on reclassification. Emits
`tierO_flag`: `exact` / `same_species` / `same_genus` / `mismatch` / `uncomparable`.

A genus-level `mismatch` on a row Tier 2/S could not confirm becomes the
`organism_mismatch` review label — a wrong-protein / mis-deposit candidate (9 rows).
Everything genus-or-better is organism-consistent (necessary, not sufficient): 166 of
the 186 remaining `unconfirmed` rows pass this check, and none of the 115 confirmations
fail it. The PAZy taxid is itself a claim, so `mismatch` means "sources disagree", never
proof.

## Tier S — the citation-graph ceiling breaker

Tier 2 can only confirm a row when NCBI links the accession to the cited paper,
which never happens for a mining/screening paper that merely *analysed* a sequence
someone else deposited. Those papers instead list the accessions and sequences they
studied in their **supplementary files**. Tier S:

1. pulls every resolvable protein's `aa_sequence` from NCBI (`ncbi.efetch_fasta`);
2. downloads + parses the cited papers' supplementary material (`supp.py`) and
   matches each row's accession/sequence against it — `seq_in_supp` / `acc_in_supp`
   are high-confidence confirmations Tier 2 structurally cannot produce;
3. recovers a live replacement accession for suppressed records via the Identical
   Protein Group (`ncbi.efetch_ipg`), since their `replacedby` pointer is empty;
4. flags rows sharing an identical sequence (duplicate curation).

This moved 82 rows out of the `unconfirmed` bucket (270 → 192) and gave 27 of the
54 suppressed records a recovered accession. The bulk (79 rows) comes from one 2025
ACS Catalysis mining paper whose SI spreadsheet lists 472 accessions / 954
sequences.

### Reproducibility caveat — supplementary fetching

Supplementary files are fetched **direct from the publisher** because NCBI's FTP OA
package and the PMC `/bin/` path are both unreachable from the build sandbox (FTP
blocked, HTTPS mirror 404s). Only ACS (`pubs.acs.org`) serves SI unauthenticated,
and it **rate-limits/blocks bursts**: after a few rapid requests it answers 200 with
an HTML challenge page instead of the file. `supp._curl` rejects those pages (magic-
byte check) rather than caching garbage, and paces requests, but a cold run can
still be throttled. The one load-bearing file — `cs5c03460_si_002.xlsx` — is
therefore committed under `.supp_cache/files/` so the headline result reproduces
offline. `supp_unavailable` in the output means "no parseable SI reachable for any
cited DOI", i.e. silence, never evidence of error.

## Tier A — sequence-driven re-mapping with an author tiebreaker

Every tier above *flags* a bad accession; none says what to point at instead. Tier A
closes that gap for the rows that are stuck with no candidate set:

- `obsolete_record` rows whose Identical Protein Group is itself entirely suppressed,
  so Tier S recovered no live replacement;
- `organism_mismatch` / `title_unmapped` rows — live accessions where the record's
  organism or annotation disagrees with PAZy (candidate wrong-protein).

For each it turns the row's sequence into a ranked candidate set and picks a best
re-map target:

1. candidates come cheap-first — Identical Protein Group members (exact, `ncbi.efetch_ipg`),
   then BLAST near-identical neighbours (`blast.blastp`, `refseq_protein` then `nr`,
   gated at ≥98% identity / ≥90% coverage);
2. each candidate is resolved (`ncbi.efetch_docsum`), non-live/suppressed ones dropped,
   then scored on **organism concordance** (Tier O's comparator) and **cited-author
   overlap** — does the candidate's own linked paper share an author surname with the
   row's cited DOI, the tiebreaker that separates the right protein from a look-alike;
3. it also measures the *existing* accession's author overlap, which settles the
   organism_mismatch / title_unmapped rows the other way: if the flagged record's own
   literature shares an author with the cited paper, the flag is benign
   (`self_corroborated`).

`tierA_flag`: `remap_strong` (live ≥99% id, organism genus-or-better, author-corroborated)
/ `remap_plausible` (a live ≥98% candidate, corroboration incomplete) / `self_corroborated`
/ `remap_none`. Tier A **never re-labels** a row — an obsolete accession must still be
re-pointed by a human — it attaches `remap_accession` + evidence and leaves the row
reviewable; a `remap_strong` row simply floats to the top of its label in the queue.

### Reproducibility caveat — BLAST

BLAST is a live remote service (the NCBI URL API, queued jobs that cost minutes each).
Every search is cached under `.blastcache/` keyed by `sha256(sequence)+database`, so a
rerun is free and offline, but a **cold run needs network + time** (~40 searches). An
empty search is cached too: `remap_none` means "no near-identical live neighbour found",
i.e. silence, never proof the row is unfixable. `python3 blast.py` runs a one-sequence
self-test; `TIERA_LIMIT=N python3 tierA_remap.py` caps Tier A to the first N stuck rows
for a smoke test.
