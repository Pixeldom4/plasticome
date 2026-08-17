# Step 6 — the 100% non-redundant set

The union collapsed to one row per distinct amino-acid sequence, so a downstream
`nr` search can ask, once per sequence, which other accessions retrieve it.

```
python "lib/06 nr/build_nr.py"  runs/<run>/01-union.tsv     # 06-nr.tsv, 06-nr.fasta
python "lib/06 nr/crosswalk.py" runs/<run>/06-nr.tsv        # 06-nr-to-clusters.tsv
```

or through the driver, which runs both:

```
./pipeline_run.bash --run-dir runs/<run> --only 6
```

## The three sets

| set | file | rows | built by | purpose |
| --- | --- | --- | --- | --- |
| A | `01-union.tsv` | 609 | step 1 | the starting set of record |
| B | `04-<run>.fasta` | 411 | steps 2–4 | 90% centroids, the node set for the component partition |
| C | `06-nr.tsv` | 493 | step 6 | 100% non-redundant, for resolving alternative accessions |

**B and C are parallel branches off step 1, not a chain.** C's only input is
`01-union.tsv`; it never feeds steps 2 or 3. Building B on top of C would shrink
the database those steps search, and a smaller database means smaller e-values
for the same alignment score, so the fixed `evalue < 1e-5` post-filter would
become effectively more permissive: more edges, fewer and larger components. That
is the same database-size confound that made
`runs/2026-08-05.align-before-cluster` uninterpretable, and with
`EVALUE_SCALE = 1.0` nothing holds the cutoff fixed across the two sizes. The
exposure differs by engine — usearch's step 2 is identity-only so only step 3 is
affected, while the diamond arm's phase-0 substrate carries DIAMOND's implicit
e-value default and is exposed at step 2 as well.

B is the delivered artifact and its lineage stays frozen. C is built beside it.

Step 6 is therefore not in the driver's default range. A plain run still produces
exactly steps 1–4; ask for C with `--only 6` or `--to 6`.

## Why hashing rather than clustering at 100%

`usearch cluster_fast -id 1.0` and `cd-hit -c 1.0` are not exact-identity
operators. Both measure identity over an alignment with a coverage model, so a
300 aa protein wholly contained in a 320 aa protein collapses at 1.0 under
defaults. That is a domain match, not "the same amino acid sequence", and it is
the wrong relation for an accession lookup. A hash groupby is O(n),
deterministic, needs no binary, and removes the Docker dependency for this step
entirely.

`normalize()` is imported from `lib/03 alignment/scripts/step1_nodes.py` rather
than redefined. Two identity spaces would diverge on exactly the rows that differ
by a trailing `*`, lowercase masking or stray whitespace, and the join between C
and everything else would fail silently. `provenance.json` records the md5 of
that function's own source text, so an edit to it surfaces there rather than as
an unexplained change in the group count.

## Files

### `06-nr.tsv`

`rep_plasticome_id, rep_accession, n_members, member_plasticome_ids, member_accessions, seq_len, aa_sequence, seq_md5`

- Representative is the **smallest `plasticome_id`** in the md5 group, matching
  step 3's smallest-PL-represents convention.
- `member_plasticome_ids` and `member_accessions` are positionally parallel
  `;`-delimited lists in ascending `plasticome_id`, so the representative is
  first and `zip()` recovers which accession belongs to which union row. A row
  with no accession leaves an **empty slot** rather than shortening the accession
  list; 41 of the 609 union rows have none, so dropping blanks would silently
  misalign the two lists. A naive split of `member_accessions` therefore yields
  empty tokens.
- **Accession versions are kept.** Step 2 folds seeds version-insensitively,
  which is right for seeds and wrong here: `WP_012345.1` and `.2` can be
  different sequences, and step 6 is precisely about accessions that retrieve an
  identical one.
- `aa_sequence` is the *normalized* sequence, so what was hashed and what gets
  queried are the same bytes.

### `06-nr.fasta`

Headed on `seq_md5` alone. The hash is the only handle that is true of a group
independent of any labelling choice, which is the property that makes C worth
building on a parallel branch. Everything else is one join away in the TSV, and
the TSV stays authoritative. This is the handoff artifact for the nr search.

### `06-nr-to-clusters.tsv` (optional)

`seq_md5, cluster_id, component_id, is_centroid, centroid_identifier, engine, identity`

`cluster_id`, `component_id` and `is_centroid` are deliberately **not** in
`06-nr.tsv`. They are branch-B facts and branch B depends on an unpinned
`--engine` choice, whereas C's single most valuable property is that it is the
one set in the pipeline whose row count does not move when the engine does. They
live here instead, in a file that carries the provenance of the branch it was
built against. `crosswalk.py` exits 0 with a note when steps 2 and 3 are absent,
which is what keeps `--only 6` runnable in a directory holding nothing but step 1.

Two joins, for two different reasons:

- **The centroid is matched on `seq_md5`.** Both sides carry the sequence, so
  both get hashed. This bypasses the `rep_plasticome_id` / `rep_label` /
  accession resolution cascade in `clusters_to_fasta.py` entirely, including the
  accession route that renamed union row 18 to `PL17`.
- **Membership is matched on label position.** `02`/`03` carry `rep_aa_sequence`
  for the centroid only, so there is no sequence on the B side to hash for a
  non-centroid row. Membership comes from `rep_label` plus `member_labels`, whose
  `U####` prefix is the union row's 1-based position. That positional map is
  checked to be a partition of the union before it is trusted.

### `06-nr.intermediates/provenance.json`

Input name and file md5, the source path and source-text md5 of `normalize`, row
and group counts, largest group, and the two settled policy flags.

## Checks

Every one of these is a hard failure. None is cosmetic.

**`build_nr.py`**

| check | why |
| --- | --- |
| no blank `aa_sequence` | step 1 leaves zero blanks on current inputs, so a blank means `fetch_sequences.py` did not complete |
| `plasticome_id` unique | it is the join key for every downstream set |
| no `;` in any accession | the delimiter is a hard error, not something to escape |
| `sum(n_members)` == union rows | the universe must be conserved |
| `--expect N` (optional) | turns the fixture number into an enforced contract |

**`crosswalk.py`**

| check | why |
| --- | --- |
| labels are a partition of the union | if the positional map is not a partition, containment means nothing |
| **md5 containment** | every group falls inside exactly one cluster and one component; a violation means the clustering split identical sequences, which invalidates B rather than relabelling it |
| centroid md5s distinct | what makes the B-to-C join one-to-one on the centroid side |
| every C row gets a cluster | coverage |

Containment properly belongs as a **step 2 postcondition**, so that it fails
where the fault is introduced rather than three steps downstream in an optional
file. Making it an `06` precondition would give C a dependency on
`02-clusters.tsv`, and therefore on the engine choice, destroying the property
that made C worth building on a parallel branch. It is re-asserted here because
it is free. Adding it to step 2 is still open.

The mechanism for a real containment miss is not the threshold — identical
sequences are trivially above 0.90 — but usearch's heuristic seed search
terminating early: `cluster_fast` abandons the search after `-maxrejects`
consecutive non-matches, so a query can be made a new centroid without ever
reaching its identical twin. `crosswalk.py`'s error message says so. The diamond
arm is structurally safer, since its phase-0 substrate is all-vs-all under
`--query-or-subject-cover 90` and identical sequences have 100% coverage in both
directions.

## Fixtures, measured against `runs/2026-08-06.final-usearch.2`

- 609 union rows → **493** distinct sequences.
- **109** duplicate groups covering **225** union rows; largest group 5.
  225 − 109 = 116 eliminated, and 609 − 116 = 493.
- Crosswalk: 493 sequences over 411 clusters and 46 components; 411 are branch-B
  centroids, 82 are not.
- **0** groups spanning more than one cluster or component.
- All 411 `centroid_identifier` values, resolved purely through the md5 route,
  agree with the delivered `04-final-usearch.2.tsv`.

**493 is a fixture number.** Any future run against this union that does not
produce 493 rows means the union or `normalize()` changed. Treat it the way the
42-component figure is treated in the pipeline notes.

## Scope: internal alternates only

Step 6 does no network lookup. Alternates come **only** from other union rows
sharing the same `seq_md5`, so every accession it reports is already in
`01-union.tsv`. No IPG, no UniParc, no blastp, no HTTP client, no cache, no rate
limiting, no query date. Step 6 is a pure function of `01-union.tsv`:
deterministic, offline, byte-reproducible, sub-second. Nothing about it goes
stale, so there is no snapshot question.

The nr search is downstream and out of scope; `06-nr.fasta` is the handoff. When
its results come back they merge in as a **second tranche**, and that merge
belongs in a separate step so the offline half stays reproducible on its own.

### Note for whoever runs the nr search

`blastp` has no `-perc_identity`, so exact retrieval is a post-filter, and
`pident == 100` alone is not sufficient: a 100% identity HSP over partial
coverage is a domain match or a fragment. It needs
`pident == 100 AND qcovhsp == 100 AND scovhsp == 100`, plus `-max_target_seqs 0`.
This is the same coverage trap already documented in steps 2 and 3. NCBI
Identical Protein Groups is exact by construction and far cheaper for anything
carrying an NCBI protein accession.

## Open

- Fill `genbank_accessions` in step 4 via a `--nr-alternates 06-nr.tsv` option on
  `clusters_to_fasta.py`, joining on the md5 of the centroid's
  `rep_aa_sequence`. Not done: it is a cross-branch dependency, C into B, and it
  gives `04-<run>.fasta` two possible byte-states under one name.
- **Field name.** Settled as `alt_accessions` rather than `genbank_accessions`,
  since the field will hold UniProt and PDB accessions too. Not yet applied —
  `04-<run>.tsv` still says `genbank_accessions`. The FASTA header is positional,
  so only the TSV column name moves.
- md5 containment as a step 2 postcondition.
- Comparison run of step 2 with `-maxrejects 0`, diffed against the delivered 411
  centroids, to decide whether containment is structural or stays contingent.
- `rep_plasticome_id` is blank on all 411 rows of `02-clusters.tsv`, so route 1
  in `clusters_to_fasta.py` is dead and resolution silently falls through to the
  `U####` prefix. Populate it in step 2 or delete route 1.
- `normalize()` exists twice, byte-identical, in `step1_nodes.py` and
  `lib/04 fasta/tsv_to_fasta.py`. Step 4 imports the second, step 6 the first, and
  step 4 is where the two sets are eventually meant to join on md5.
- Step 5 has no code. Step 6 stays opinion-free about delivered-table schema so a
  future `lib/05 annotate/` can consume `06-nr.tsv` rather than duplicate it.
