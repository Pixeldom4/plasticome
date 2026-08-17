# Step 6 — the 100% non-redundant set

The union collapsed to one row per distinct amino-acid sequence, so a downstream
`nr` search can ask, once per sequence, which other accessions retrieve it.

```bash
./pipeline_run.bash --run-dir runs/<run> --only 6     # runs both scripts
python "lib/06 nr/build_nr.py"  runs/<run>/01-union.tsv
python "lib/06 nr/crosswalk.py" runs/<run>/06-nr.tsv
```

## Why it is a parallel branch

| set | file | rows | purpose |
|---|---|---|---|
| A | `01-union.tsv` | 609 | starting set of record |
| B | `04-<run>.fasta` | 411 | 90% centroids, node set for the component partition |
| C | `06-nr.tsv` | 493 | 100% non-redundant, for resolving alternative accessions |

C's only input is `01-union.tsv`, and it never feeds steps 2 or 3. Building B on
top of C would shrink the database those steps search, and a smaller database
gives smaller e-values for the same score, so the fixed `evalue < 1e-5` post-filter
would become more permissive: more edges, fewer and larger components. That is the
confound that made `runs/2026-08-05.align-before-cluster` uninterpretable, and with
`EVALUE_SCALE = 1.0` nothing holds the cutoff fixed across the two sizes. B is the
delivered artifact and its lineage stays frozen; C is built beside it.

Step 6 is therefore opt-in, not in the driver's default 1-4 range.

## Why hashing, not clustering at 100%

`cluster_fast -id 1.0` and `cd-hit -c 1.0` are not exact-identity operators: both
measure identity over an alignment with a coverage model, so a 300 aa protein
wholly inside a 320 aa one collapses at 1.0. That is a domain match, not the same
sequence, and it is the wrong relation for an accession lookup. Hashing is O(n),
deterministic, and needs no binary or Docker.

`normalize()` is imported from `lib/03 alignment/scripts/step1_nodes.py`, not
redefined. Two identity spaces would diverge on exactly the rows differing by a
trailing `*`, lowercase masking or whitespace, and the join would fail silently.
`provenance.json` records the md5 of that function's own source text.

## Files

**`06-nr.tsv`** —
`rep_plasticome_id, rep_accession, n_members, member_plasticome_ids, member_accessions, seq_len, aa_sequence, seq_md5`

- Representative is the **smallest `plasticome_id`** in the md5 group, matching
  step 3's convention.
- The two `member_*` lists are **positionally parallel**, `;`-delimited, ascending
  `plasticome_id`, so the representative is first and `zip()` recovers which
  accession belongs to which row. A row with no accession leaves an **empty slot**;
  41 of 609 have none, so dropping blanks would misalign the lists. A naive split
  therefore yields empty tokens.
- **Accession versions are kept.** Step 2 folds seeds version-insensitively, which
  is right for seeds and wrong here: `.1` and `.2` can be different sequences.
- `aa_sequence` is the *normalized* one, so what was hashed and what gets queried
  are the same bytes.

**`06-nr.fasta`** — headed on `seq_md5` alone, the only handle true of a group
independent of any labelling choice. The handoff artifact for the nr search.

**`06-nr-to-clusters.tsv`** (optional) —
`seq_md5, cluster_id, component_id, is_centroid, centroid_identifier, engine`

Branch-B facts are deliberately kept out of `06-nr.tsv`, since they depend on an
unpinned `--engine` while C's row count does not. `crosswalk.py` exits 0 with a
note when steps 2 and 3 are absent, which keeps `--only 6` runnable in a directory
holding only step 1.

Two joins, for two reasons. The **centroid is matched on `seq_md5`**, bypassing the
`rep_plasticome_id` / `rep_label` / accession cascade that once renamed union row
18 to `PL17`. **Membership is matched on label position**, because 02/03 carry
`rep_aa_sequence` for the centroid only, so there is no sequence to hash for a
non-centroid row.

**Sidecars** — `provenance.json` (input md5, `normalize` source md5, counts) and
`crosswalk.json` (branch B built against, engine, identity, counts, checks passed).
Identity lives here rather than as a column: it is a per-run constant and was `0.9`
on every run to date.

## Checks

All hard failures.

| script | check |
|---|---|
| `build_nr.py` | no blank `aa_sequence`; `plasticome_id` unique; no `;` in an accession; `sum(n_members)` equals union rows; optional `--expect N` |
| `crosswalk.py` | labels are a partition of the union; **md5 containment**; centroid md5s distinct; every C row gets a cluster |

**md5 containment** means every group falls inside exactly one cluster and one
component. A violation means the clustering split identical sequences, which
invalidates B rather than relabelling it. It properly belongs as a step-2
postcondition so it fails where the fault originates; it is re-asserted here
because it is free. Adding it to step 2 is still open.

The mechanism for a real miss is not the threshold, since identical sequences are
trivially above 0.90, but usearch abandoning its seed search after `-maxrejects`
consecutive non-matches. The diamond arm is structurally safer, its phase-0
substrate being all-vs-all under `--query-or-subject-cover 90`.

## Fixtures

609 union rows → **493** distinct sequences. 109 duplicate groups covering 225
rows, largest 5 (225 − 109 = 116 eliminated; 609 − 116 = 493). Crosswalk: 493 over
411 clusters and 46 components, 411 centroid and 82 not. Zero groups span more than
one cluster or component. All 411 `centroid_identifier` values, resolved purely
through md5, agree with the delivered `04-*.tsv`.

**493 is a fixture number.** A run against this union not producing it means the
union or `normalize()` changed.

## Scope

No network lookup. Alternates come only from other union rows sharing a `seq_md5`,
so every accession reported is already in `01-union.tsv`. Step 6 is a pure function
of that file: deterministic, offline, byte-reproducible, sub-second.

The nr search is downstream and out of scope. Its results merge in as a **second
tranche**, which belongs in a separate step so the offline half stays reproducible.
Note for whoever runs it: `blastp` has no `-perc_identity`, and `pident == 100`
alone is insufficient since a 100% HSP over partial coverage is a domain match. It
needs `pident == 100 AND qcovhsp == 100 AND scovhsp == 100`, plus
`-max_target_seqs 0`. NCBI Identical Protein Groups is exact by construction and
far cheaper for anything with an NCBI protein accession.

## Open

- Fill `genbank_accessions` in step 4 via `--nr-alternates 06-nr.tsv`. Not done: it
  is a cross-branch dependency and gives `04-<run>.fasta` two byte-states under one
  name.
- Field name settled as `alt_accessions`, not yet applied; `04-<run>.tsv` still says
  `genbank_accessions`.
- md5 containment as a step-2 postcondition.
- Step 2 comparison run with `-maxrejects 0`, diffed against the delivered 411.
- `normalize()` exists twice, byte-identical, in `step1_nodes.py` and
  `lib/04 fasta/tsv_to_fasta.py`. Step 4 imports the second, step 6 the first.
