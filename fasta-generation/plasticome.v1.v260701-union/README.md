# plasticome.v1.v260701-union

The v1 ∪ v2 union set of the plasticome (611 enzymes), plus the protein FASTA
derived from it and a seed-free subset for downstream homology search.

## Files

| File | Records | Description |
|---|---|---|
| `plasticome_v1.v260701-union.csv` | 611 | **Source of truth.** Columns: `plasticome_id, enzyme_name, accession, pazy_id, aa_sequence, source`. |
| `plasticome.v1.v260701-union.tsv` | 611 | Tab-separated copy of the CSV (same columns), used as input to the FASTA generator. |
| `plasticome.v1.v260701-union.fasta` | 611 | Protein FASTA. Header: `>{plasticome_id}\|{accession}\|\|{pazy_id}`, sequences wrapped at 60 columns. |
| `pazy-seed-subset/putative-curated-seed-set.csv` | 13 (+4 empty rows) | Hand-curated seed enzymes to exclude. Columns: `enzyme_name, accession, pazy_id, justification, aa_sequence`. |
| `pazy-seed-subset/plasticome.v1.v260701-union.no-seeds.fasta` | 598 | The master FASTA with the 13 seeds removed (611 − 13). |

## Header schema

```
>{plasticome_id}|{accession}||{pazy_id}
```

- `accession` is a `;`-joined list of **protein/stable** identifiers only —
  GenBank/RefSeq **nucleotide** accessions are stripped (see the generator's
  `is_nucleotide`). The empty third field is a reserved, always-empty `genbank`
  slot kept for column alignment across the repo.
- This matches the `plasticome.v1.1` convention (nucleotides stripped), **not**
  the `singleton-cleaned-union.v1.v2` convention (nucleotides kept).
- 42 records carry an empty accession field: 39 were already empty in the
  source CSV, and 3 (PL137/138/139) had only nucleotide accessions, which were
  stripped.

## Workflow

### 1. CSV → TSV

The FASTA generator consumes a TSV. Convert with quote-safe handling:

```python
import csv
with open('plasticome_v1.v260701-union.csv', newline='') as f, \
     open('plasticome.v1.v260701-union.tsv', 'w', newline='') as o:
    w = csv.writer(o, delimiter='\t', lineterminator='\n')
    for row in csv.reader(f):
        w.writerow(row)
```

### 2. TSV → FASTA

Uses the repo's nucleotide-stripping generator, one level up in
`fasta-generation/`. `plasticome_id` is the header identifier:

```bash
python3 ../tsv_to_fasta_no_nucleotide.py \
  plasticome.v1.v260701-union.tsv \
  --header-template "{plasticome_id}|{accession}||{pazy_id}" \
  -o plasticome.v1.v260701-union.fasta
```

Reports: `wrote 611 records from 611 rows (0 blank seq; dropped 3 nucleotide
accessions; 3 rows left with an empty accession field)`.

### 3. FASTA → seed-free subset

Remove the curated seeds in `pazy-seed-subset/putative-curated-seed-set.csv`
from the master FASTA. Matching is by **exact accession**: each seed's
`accession` is matched against the header's `;`-separated accession list, and
matching records are dropped. This removes exactly 13 records (one per seed) →
598 kept.

```bash
python3 make_no_seeds.py    # see snippet below
```

```python
import csv
# 1. collect seed accessions (skip the empty rows)
seedset = set()
for r in csv.DictReader(open('pazy-seed-subset/putative-curated-seed-set.csv', newline='')):
    if (r.get('aa_sequence') or '').strip() and (r.get('accession') or '').strip():
        seedset.add(r['accession'].strip())

# 2. read master FASTA
def read_fasta(p):
    hdr, buf = None, []
    for line in open(p):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if hdr is not None: yield hdr, ''.join(buf)
            hdr, buf = line[1:], []
        else: buf.append(line)
    if hdr is not None: yield hdr, ''.join(buf)

def header_accs(hdr):                       # 2nd pipe-field, ';'-split
    parts = hdr.split('|')
    return {t.strip() for t in (parts[1] if len(parts) > 1 else '').split(';') if t.strip()}

def wrap(s, w=60):
    return '\n'.join(s[i:i+w] for i in range(0, len(s), w))

# 3. write complement (records whose accessions don't intersect the seed set)
with open('pazy-seed-subset/plasticome.v1.v260701-union.no-seeds.fasta', 'w') as o:
    for hdr, seq in read_fasta('plasticome.v1.v260701-union.fasta'):
        if not (header_accs(hdr) & seedset):
            o.write(f'>{hdr}\n{wrap(seq)}\n')
```

## Notes & caveats

- **Matching is by exact accession, not sequence.** 9 of the 13 seeds also have
  a *sequence-identical duplicate* elsewhere in the master (under a different
  accession, e.g. LCC's `AEV21261.1` at PL16 duplicated by `G9BY57.1` at PL522).
  Those duplicates **remain** in the no-seeds set. Switch to normalized-sequence
  matching if a fully seed-sequence-free reference is required (that removes 22).
- Cut190's seed accession `BAO42836.1` matched **PL485**, not PL18 — PL18 carries
  the RefSeq `WP_015787089.1` for the same enzyme.
- Version suffixes matter under exact matching: seed `CBY05530` (TfCut_2) matched
  PL24 but **not** PL602 `CBY05530.1`. Version-insensitive matching would catch
  both.

## Removed seeds (13)

ISPETase (PL1), LCC (PL16), TfCut_2 (PL24), Cut190 (PL485), HiC (PL43),
PueA (PL48), NylA (PL51), NylB (PL52), NylC (PL53), PLAase I (PL59),
PhaZ1 (PL85), Mipa-P (PL239), Kubu-P (PL245).
