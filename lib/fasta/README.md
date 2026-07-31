# lib/fasta

TSV → protein FASTA generation, plus the accession-column order checker and its
fixer. Two generators, split by *function* rather than by input schema — both
read either source sheet through flags.

| Script | Role |
|---|---|
| `tsv_to_fasta.py` | Generic engine. Any TSV, any columns. Adds `--dedup`, `--desc-cols`/`--sep`, `--width`. |
| `tsv_to_fasta_no_nucleotide.py` | **Production.** Same output conventions, but `{accession}` is filtered to protein/stable identifiers — GenBank/RefSeq *nucleotide* accessions are stripped. This is the `plasticome.v1.1` convention. |
| `validate_accession_order.py` | QC: checks the `;`-separated accession field is ordered protein → cross-ref → PDB. |
| `fix_accession_order.py` | Applies the hand-verified fixes `validate_accession_order.py` flags. Run it from inside this directory. |

Shared conventions: `normalize()` strips every non-letter and upper-cases (gaps,
`*` stops and whitespace never reach a record), sequences wrap at 60 columns,
blank-sequence rows are dropped rather than emitted as empty records, and
`--header-template` keeps empty placeholders so literal separators — the
always-empty genbank pipe in `PL18|WP_015787089.1||18` — survive.

`tsv_to_fasta_pl_headers.py` was retired on 2026-07-30; it duplicated the
header format with three behavioral differences that were never wanted
downstream (no nucleotide filtering, no sequence normalization, no line
wrapping).

## Canonical invocations

### `source-data/plasticome.v260701/` — cleaned PAZy sheet

Columns: `enzyme_name, accession, direct_or_indirect_activity, pazy_id, doi,
aa_sequence, organism, ncbi_taxonomy_id, tier_s_label, sequence_source,
dna_acc_to_pro_acc, comment`.

```bash
python3 lib/fasta/tsv_to_fasta_no_nucleotide.py \
  source-data/plasticome.v260701/cleaned_pazy-260701_retrieving_from_accession.tsv \
  --header-template "{enzyme_name}|{accession}||{pazy_id}" \
  -o out.fasta
```

473 rows. 27 have an empty `accession`, 5 an empty `pazy_id` — both render as
empty header fields, not as dropped records.

> **This sheet currently ships a duplicate `aa_sequence` column** (column 6 holds
> the sequences; column 12 is empty on 469 of 473 rows). `csv.DictReader` keeps
> only the *last* column of a repeated name, so reading it yields 4 records out
> of 473. Both generators now detect duplicate headers and exit 1 rather than
> write a truncated FASTA. Drop column 12 from the source to run.

### `source-data/plasticome.v1.1/` — re-annotated v1.1 sheet

Columns: `rowid, retrieved, component, gene, domain, cath, kingdom, host,
enzyme_name, substrate, EC, reference, pazy_accession, BLAST_link,
pdb_accession, uniprot_accession, retrieved_name, retrieved_host, synthesized,
enzymatic_activity, aa_sequence`.

There is no column named `accession`; the NCBI accession lives in `retrieved`.
Point `--acc-col` at it — the filtered list is still exposed to the template as
`{accession}`.

```bash
python3 lib/fasta/tsv_to_fasta_no_nucleotide.py \
  source-data/plasticome.v1.1/plasticome.v1.1.tsv \
  --acc-col retrieved \
  --header-template "{rowid}|{accession}||" \
  -o out.fasta
```

212 rows.

> **Neither schema carries a `plasticome_id` or `pazy_id` usable as the PL
> number**, so the `>PL18|WP_015787089.1||18` header used by the clustering runs
> cannot be produced from these sheets directly — it needs a back-join to the
> 611-row union table (`runs/2026-07-23_plasticome.v1.v260701-union/fasta/`).
> That join is ambiguous: the union lists the same protein twice, once per
> source. Cut190 is `PL18 | WP_015787089.1 | pazy 18` on the v2 side and
> `PL485 | BAO42836.1 | (no pazy_id)` on the v1 side, and the v1.1 sheet's
> `retrieved` value matches the v1 side. Decide which side is authoritative
> before generating PL-numbered headers.

### Reproducing the shipped `plasticome.v1.1.fasta`

The 611-record FASTA came from the *committed* v1.1 table (`plasticome_id,
enzyme_name, accession, pazy_id, aa_sequence, source`), not from the
re-annotated sheet above:

```bash
python3 lib/fasta/tsv_to_fasta_no_nucleotide.py plasticome.v1.1.tsv \
  --header-template "{plasticome_id}|{accession}||{pazy_id}" \
  -o plasticome.v1.1.fasta
# wrote 611 records from 611 rows (0 blank seq; dropped 3 nucleotide
# accessions; 3 rows left with an empty accession field)
```

Verified byte-identical to the committed FASTA. See
`runs/2026-07-23_plasticome.v1.v260701-union/fasta/README.md` for the union
build and the seed-removal step.

## Accession-order QC

Both default their positional argument to `cleaned_pazy_final.tsv`, which no
longer exists at that bare path after the reorg — pass the TSV explicitly.

```bash
python3 lib/fasta/validate_accession_order.py <tsv> --acc-col accession --strict
cd lib/fasta && python3 fix_accession_order.py <tsv> --dry-run
```

`fix_accession_order.py` imports `validate_accession_order` and re-invokes it by
bare filename, so it only works with `lib/fasta/` as the working directory. Its
edits are hardcoded and hand-verified; it aborts if a target row no longer
matches the recorded `before` value.
