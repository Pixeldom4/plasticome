# Step 5 — the union annotated with its cluster and component

Steps 2 to 4 are cluster-level: 411 rows standing in for 609 sequences, with the
198 non-centroid rows visible only as labels inside `member_labels`. Step 5 puts
that back at union-row level, so every sequence states where it landed.

```
python "lib/05 annotate/annotate_union.py" runs/<run>/03-alignment.tsv \
    --union runs/<run>/01-union.tsv -o runs/<run>/05-union-with-components.tsv
```

or through the driver:

```
./pipeline_run.bash --run-dir runs/<run> --only 5
```

Like step 6, step 5 is in the driver's `1-6` range but not in the default `1-4`,
so a plain run still produces exactly steps 1 to 4.

| in | out |
| --- | --- |
| `01-union.tsv`, 609 rows, one per sequence, no cluster facts | `05-union-with-components.tsv`, 609 rows, one per sequence, with them |
| `03-alignment.tsv`, 411 rows, one per cluster, membership as labels | |

Row order is the union's, unchanged, so `05` and `01` are row-aligned and can be
pasted side by side. The first five columns are `01-union.tsv` verbatim.

## Columns

| column | source |
| --- | --- |
| `plasticome_id` `enzyme_name` `accession` `pazy_id` `source` | the union, verbatim |
| `seq_len` | length of the normalized sequence |
| `cluster_id` `component_id` | the step-3 table |
| `is_centroid` | `yes` / `no` |
| `pct_id_to_centroid` | `member_pct_ids`, as the aligner reported it |
| `centroid_identifier` `centroid_label` `centroid_accession` | the centroid this row sits under |
| `aa_sequence` | the union, verbatim |

Two conventions worth knowing, because neither is a measured value:

**`pct_id_to_centroid` is `100.0` on a centroid.** The clusters table stores one
`member_pct_ids` entry per *non-representative* member and nothing for the
representative, since its identity to itself is not something the aligner
reports. `100.0` is this script's convention for that empty slot. It is exact
rather than approximate: a centroid is its own sequence.

**`aa_sequence` is passed through unnormalized**, unlike step 6, which emits the
normalized form because a hash is only meaningful over the bytes that were
hashed. Step 5 is the union with columns added, so its sequence column stays the
union's. On the current union the two are identical anyway, since step 1 already
writes clean sequences; the distinction only matters if that stops being true.

## Identifier resolution

`centroid_identifier` is `PL` plus the `plasticome_id` of the union row that the
centroid's `U####` label points at, so `05.centroid_identifier` joins to
`04.identifier` exactly.

It is resolved **positionally**, not from `rep_plasticome_id`, which is blank on
all 411 rows of the current clusters table, and not from the centroid's
accession, which is the route that renamed union row 18 to `PL17`.

Membership comes from [`lib/common/membership.py`](../common/membership.py),
shared with step 6, which verifies the positional map is a partition of the union
before returning it: every union row in exactly one cluster, no row in two, no
label without a `U####` prefix, no label pointing outside the union, each
cluster's label count equal to its stated `size`, and one `member_pct_ids` entry
per non-representative member. Every fault it finds is reported together rather
than one run at a time.

## Validation

`05-union-with-components.tsv` on `runs/2026-08-06.final-usearch.2` was built
outside this repo before `lib/05` existed. This script reproduces that file
**byte for byte** from `01-union.tsv` and `03-alignment.tsv` alone, all 14
columns and all 609 rows, which is what established that the reverse-engineered
schema is the right one.

609 rows over 411 clusters and 46 components: 411 centroids, 198 non-centroid
members.

## Open

- The design note for step 6 calls for an `alt_accessions` column here too, at
  union-row level, so all 609 rows carry the alternates for their own `seq_md5`
  group rather than only the centroids. Not added: it would change the column set
  and break the byte-for-byte reproduction above, which is currently the only
  proof the schema is right. It wants doing at the same time as the step-4 fill,
  since the two share the `06-nr.tsv` join and the same naming decision.
