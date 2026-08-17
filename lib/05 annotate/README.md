# Step 5 — the union annotated with its cluster and component

Steps 2 to 4 are cluster-level: 411 rows standing in for 609 sequences, with the
198 non-centroid rows visible only as labels inside `member_labels`. Step 5 puts
that back at union-row level, so every sequence states where it landed.

```bash
./pipeline_run.bash --run-dir runs/<run> --only 5
python "lib/05 annotate/annotate_union.py" runs/<run>/03-alignment.tsv \
    --union runs/<run>/01-union.tsv -o runs/<run>/05-union-with-components.tsv
```

In: `01-union.tsv` (609 rows) and `03-alignment.tsv` (411 rows). Out:
`05-union-with-components.tsv`, 609 rows. Opt-in, like step 6, not in the driver's
default 1-4 range.

Row order is the union's, so 05 and 01 are row-aligned. The first five columns are
`01-union.tsv` verbatim.

## Columns

| column | source |
|---|---|
| `plasticome_id` `enzyme_name` `accession` `pazy_id` `source` | the union, verbatim |
| `seq_len` | length of the normalized sequence |
| `cluster_id` `component_id` | the step-3 table |
| `is_centroid` | `yes` / `no` |
| `pct_id_to_centroid` | `member_pct_ids`, as the aligner reported it |
| `centroid_identifier` `centroid_label` `centroid_accession` | the centroid this row sits under |
| `aa_sequence` | the union, verbatim |

Two values that are conventions rather than measurements:

**`pct_id_to_centroid` is `100.0` on a centroid.** The clusters table stores one
`member_pct_ids` entry per non-representative member and nothing for the
representative, whose identity to itself the aligner does not report. Exact rather
than approximate: a centroid is its own sequence.

**The same column holds the literal string `dup` on 9 rows**, marking seed
duplicates folded in by step 2. All 9 are byte-identical to their centroid, so
`dup` means 100%, but **`float()` on this column raises**. Map it first.

**`aa_sequence` passes through unnormalized**, unlike step 6, which emits the
normalized form because a hash is only meaningful over the bytes hashed. On the
current union the two are identical anyway.

## Identifier resolution

`centroid_identifier` is `PL` plus the `plasticome_id` of the union row the
centroid's `U####` label points at, so it joins to `04.identifier` exactly.

Resolved **positionally**, not from `rep_plasticome_id` (blank on all 411 rows of
runs built before 2026-08-12) and not from the accession (the route that renamed
union row 18 to `PL17`). Membership comes from
[`lib/common/membership.py`](../common/membership.py), shared with step 6, which
verifies the positional map is a partition of the union before returning it and
reports every fault together rather than one run at a time.

## Validation

`05-union-with-components.tsv` on `runs/2026-08-06.final-usearch.2` was built
outside this repo before `lib/05` existed. This script reproduces it **byte for
byte** from `01-union.tsv` and `03-alignment.tsv` alone, all 14 columns and 609
rows. That is what established the reverse-engineered schema is correct.

609 rows over 411 clusters and 46 components: 411 centroids, 198 members.

## Open

The step-6 design note calls for an `alt_accessions` column here too, at union-row
level, so all 609 rows carry the alternates for their own `seq_md5` group. Not
added: it would break the byte-for-byte reproduction above, currently the only
proof the schema is right. It wants doing alongside the step-4 fill, since the two
share the `06-nr.tsv` join and the same naming decision.
