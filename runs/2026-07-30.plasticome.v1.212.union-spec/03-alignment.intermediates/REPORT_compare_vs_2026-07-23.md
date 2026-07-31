# 2026-07-30 union-spec alignment vs the 2026-07-23 singleton-cleaned passes

**This run:** `runs/2026-07-30.plasticome.v1.212.union-spec/03 alignment.tsv.work/`
(deliverable: `03 alignment.tsv`)
411 cluster centroids → 411 md5-unique nodes (0 collapsed), db 144,131 residues.

**Compared against:** `runs/2026-07-23_singleton-cleaned-union.v1.v260701/pass1-alignment`
and `pass3-alignment` (412 centroids each, usearch **and** diamond in both).

Method is identical across all three — usearch11 v11.0.667 `-allpairs_local -acceptall`
via Docker linux/amd64, both orientations, best HSP per unordered pair by max bits,
post-filter `pct_id ≥ 30 AND e-value < 1e-5`, single-linkage components. Only the input
node set changed. Everything below is keyed on **`sequence_md5`**, since `m####` node
ids are assigned in md5-sorted order within a run and are meaningless across runs.

Regenerate with `python compare_vs_2026-07-23.py` →
`comparison_vs_2026-07-23.json`, `node_disagreements_vs_2026-07-23.tsv`.

## Headline — the partition did not move

| metric | **this run** | pass1/usearch | pass3/usearch | pass1/diamond | pass3/diamond |
|---|---:|---:|---:|---:|---:|
| nodes | **411** | 412 | 412 | 412 | 412 |
| passing edges | **22,211** | 22,207 | 22,410 | 22,321 | 22,541 |
| components | **46** | 47 | 46 | 47 | 46 |
| singletons | **25** | 26 | 25 | 25 | 24 |
| largest | **232** | 232 | 233 | 228 | 229 |
| top-10 sizes | 232, 45, 19, 16, 14, 10, 5, 5, 5, 5 | same | 233 then same | 228 then same | 229 then same |

Ranks 2–10 are **identical in all five partitions**: 45, 19, 16, 14, 10, 5, 5, 5, 5.

### Agreement on the nodes each pair actually shares

| vs | shared nodes | edge Jaccard | Rand | **ARI** | nodes in a non-corresponding component |
|---|---:|---:|---:|---:|---:|
| pass1 / usearch | 396 | 0.99995 | 1.000000 | **1.000000** | **0** |
| pass3 / usearch | 397 | 0.99995 | 1.000000 | **1.000000** | **0** |
| pass1 / diamond | 396 | 0.97996 | 0.985884 | 0.967988 | 6 |
| pass3 / diamond | 397 | 0.97863 | 0.985892 | 0.968081 | 6 |

**Against both usearch passes the partition is exactly reproduced — ARI 1.000, zero
node disagreements.** Every component-count difference (47 vs 46, 26 vs 25 singletons,
233 vs 232) is accounted for by nodes that are not in both inputs, not by any node
being regrouped.

The one edge that differs from *both* passes is the same edge each time:

```
CL0094|U0464|WKW63607.1 -- CL0272|U0186|SHM40309    pct_id 33.8   e = 9.990e-06
```

That is 0.1% under the 1e-5 cutoff, and usearch's e-value scales with database size —
144,131 residues here vs 144,774 (pass1) / 145,005 (pass3). A slightly smaller db
pushes it just inside the threshold. It is redundant anyway: both endpoints already
sit in C001, which is why ARI stays at 1.000.

## What changed in the node set

Vs pass3 (14 nodes here are new, 15 of its nodes are absent):

| category | pass1 | pass3 | what it is |
|---|---:|---:|---|
| same accession, different sequence | 3 | 4 | the union rebuild changed the residues |
| centroid representative swap | 11 | 9 | same cluster, a different member is now the rep |
| genuinely added | 1 | 1 | `GuaPA`, 299 aa, no accession |
| genuinely dropped | 2 | 2 | `PL458` (324 aa, no acc); `PL14`/`PL444` |

**Every one of these landed in the component corresponding to its predecessor's** — the
rep swaps and re-sequenced nodes all stayed put (11 of 14 in the C001 hydrolase hub, the
rest in their same-size satellites). `C###` labels themselves differ between runs because
numbering derives from identifier order, and this run's identifiers are `CL####|…` rather
than `PL###`; the labels are cosmetic, the grouping is what the ARI measures.

### The four re-sequenced accessions

| accession | 2026-07-23 | this run | Δ | note |
|---|---:|---:|---:|---|
| EGD44994 (`503`) | 294 | 267 | −27 | Erickson manual-assignment exception, per the run readme |
| WP_093412886 (`611`) | 293 | 269 | −24 | same |
| WP_085690612 (`EstB`) | 231 | 218 | −13 | pass3 only |
| **NP_999368 (porcine pancreatic lipase)** | 448 | **112** | **−336** | see below |

Three of the four are the documented Erickson/v1.1 manual assignments — expected, and
the run readme predicts them. The fourth is not.

> **`NP_999368.1` is a 112-aa fragment in this run's input.** It comes straight from
> `sources/plasticome.v260701/cleaned_pazy-260701_retrieving_from_accession.tsv`, which
> carries only 112 residues for that accession; the union passed it through unchanged
> (`fetched_sequences.json` is empty — Stage 2 resolved nothing). The 2026-07-23 runs had
> the full 448-aa record. Porcine pancreatic lipase is ~450 aa, so 112 is a truncation in
> the **source table**, not something the alignment or the union introduced. It is a
> singleton either way (C046 here, C044 in pass1), so it does not affect this partition —
> but the source row should be repaired before the sequence is used for anything else.

## usearch vs diamond — unchanged from 2026-07-23

The 6 disagreeing nodes against diamond are the *same six* the
`pass3-alignment/REPORT_compare_2026-07-23.md` identified, under this run's labels:

| this run | 2026-07-23 | usearch | diamond |
|---|---|---|---|
| `CL0359\|U0118` | PL122 | hub C001 | satellite (5) |
| `CL0370\|U0040` (RLI42440.1) | PL40 | hub C001 | satellite (5) |
| `CL0379\|U0117` | PL117 | hub C001 | satellite (5) |
| `CL0409\|U0122` | PL118 | hub C001 | satellite (5) |
| `CL0360\|U0195` (RGD93181.1) | PL195 | hub C001 | singleton |
| `CL0201\|U0454` (WDT94443) | PL454 | small comp (4) | singleton |

Same nodes, same direction, same magnitude (ARI 0.968 here vs 0.963 reported on the
412-node set). The aligner-choice question is exactly where it was left.

## Bottom line

The union rebuild and re-clustering changed **~3.5% of the node set** (14–16 sequences,
mostly representative swaps within unchanged clusters) and changed **nothing about the
partition**: ARI 1.000 against both prior usearch runs, identical size spectrum at ranks
2–10, zero regrouped nodes, one redundant threshold-boundary edge. The single finding
that warrants action is the truncated `NP_999368.1` source row.
