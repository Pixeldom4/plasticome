# PETadex clustering — 90% identity

Sequence clustering of the plasticome union set at **90% amino-acid identity** using
**USEARCH v11.0.667 `cluster_fast`**.

## Tool & parameters

```
USEARCH v11.0.667  -cluster_fast  -id 0.90  -sort length
```

- `-cluster_fast` — greedy centroid-based clustering; sequences are ordered before clustering.
- `-id 0.90` — 90% identity threshold to join a centroid.
- `-sort length` — order by decreasing length (longest sequences become centroids first).

### Running on this machine (arm64 macOS)

`usearch11` is a **Linux x86-64 ELF** binary and will not run natively on this arm64 Mac
(`exec format error`). It is run inside a Docker `linux/amd64` container. Docker Desktop must
be running first (`open -a Docker`).

Binary location: `../petadex-alignment-usearch/bin/usearch11` (v11.0.667).

```bash
docker run --rm --platform linux/amd64 \
  -v "/Users/Pixel/Documents/projects/plasticome/petadex-alignment-usearch/bin:/b" \
  -v "$PWD:/d" debian:stable-slim \
  /b/usearch11 -cluster_fast /d/<input>.fasta -id 0.90 \
    -centroids /d/centroids.fasta \
    -uc /d/clusters.uc \
    -sort length \
    -log /d/usearch.log
```

## Inputs

Two runs were performed on equivalent content (the plasticome v1–v2 union, 611 records):

| Run dir | Input | Source |
|---|---|---|
| `./` | `input.fasta` (built from TSV) | `../v1-v2-union/plasticome.v1-v2-union.tsv` |
| `./plasticome.v1.1/` | `plasticome.v1.1.fasta` (pre-built FASTA) | provided directly |

### TSV → FASTA conversion (root run only)

The union TSV has columns `enzyme_name, accession, pazy_id, aa_sequence, source`.
Neither `accession` (39 empty) nor `pazy_id` (143 empty, not unique) is a safe standalone
label, so each record was given a unique row-indexed label and a mapping table was kept:

```
label  = seqNNNN|<accession or NA>|pazy<pazy_id or NA>|<source>
```

- `*` stop characters stripped from sequences.
- `label_map.tsv` maps each FASTA label back to its TSV fields.

The `plasticome.v1.1/` run used the supplied FASTA as-is (already unique, clean labels,
no whitespace or stop characters), so no conversion or label map was needed.

## Pipeline

1. **(root run only)** Convert union TSV → `input.fasta`, writing `label_map.tsv`.
2. Run `usearch11 -cluster_fast` at `-id 0.90 -sort length` in the Docker container,
   emitting `centroids.fasta`, `clusters.uc`, and `usearch.log`.
3. Parse `clusters.uc` (`S` = centroid/seed, `H` = hit/member, `C` = cluster summary) into
   `cluster_membership.tsv` (cluster_id, size, role, label), sorted by cluster size desc.
4. Verify `S + H` record counts equal the 611 input sequences.

## Results (identical for both runs)

- **611 sequences → 411 clusters** at 90% identity.
- Singletons: **270** (65.7% of clusters).
- Largest cluster: **21** members.
- USEARCH collapses 495 unique sequences internally (109 exact-duplicate sequences in the
  union), but all 611 originals are recorded in `clusters.uc` (411 `S` + 200 `H` = 611).

Cluster size distribution (size : count):

| size | 1 | 2 | 3 | 4 | 5 | 6 | 21 |
|---|---|---|---|---|---|---|---|
| clusters | 270 | 116 | 14 | 5 | 4 | 1 | 1 |

## Output files

Per run directory:

| File | Description |
|---|---|
| `centroids.fasta` | 411 cluster representative sequences |
| `clusters.uc` | Raw USEARCH cluster records (`S`/`H`/`C`) |
| `cluster_membership.tsv` | Tidy membership: cluster_id, size, role (centroid/member), label |
| `usearch.log` | Full USEARCH run log |
| `input.fasta` + `label_map.tsv` | (root run only) FASTA built from TSV and its label map |

## Notes

- The union contains 109 exact-duplicate sequences. Clustering handles them (identical
  sequences co-cluster). To dereplicate first, or to prefer a centroid by source priority
  (e.g. v1 over v2), the run would need adjusting — not done here.
- Related memory: `usearch11/usearch12` are Docker-only on this arm64 workstation.
