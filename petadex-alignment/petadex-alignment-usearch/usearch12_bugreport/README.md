# usearch12.0-beta — evidence package (two issues vs usearch v11.0.667)

Two reproducible discrepancies between **usearch v12.0-beta1** and **usearch v11.0.667**,
found while reproducing an all-vs-all protein similarity network (213 sequences, edges at
≥30% amino-acid identity AND e-value < 1e-5).

Binaries, SHA256, provenance, and version banners: [`binaries.txt`](binaries.txt).
- v11.0.667 `i86linux64` — sha256 `a2b981fb…590ea8` — `rcedgar/usearch_old_binaries`
- v12.0-beta1 `linux_x86` — sha256 `4193abea…27f3dc` — `rcedgar/usearch12` release `v12.0-beta1`

---

## Issue 1 — missing exhaustive / e-value search commands
`claim1_missing_commands/` · rerun: `bash run_claim1.sh` · raw per-command logs in `logs/`

v12.0-beta returns **`Unknown command-line option`** for commands that v11.0.667 implements:

| command | v11.0.667 | v12.0-beta |
|---|---|---|
| `-allpairs_local` | recognized | **Unknown command-line option** |
| `-allpairs_global` | recognized | **Unknown command-line option** |
| `-ublast` | recognized | **Unknown command-line option** |
| `-calc_distmx` | recognized | **Unknown command-line option** |
| `-usearch_local` | recognized | recognized |
| `-usearch_global` | recognized | recognized |

Impact: the exhaustive all-vs-all (`allpairs_local`) and the sensitive e-value search
(`ublast`) are the natural way to build an all-vs-all similarity graph. Their absence means a
v11 pipeline cannot be reproduced on v12; only the heuristic `usearch_local`/`_global` remain
— and those disagree with the exhaustive result (and carry Issue 2).

---

## Issue 2 — `usearch_local` e-value inflated by a constant factor ~13,500 (10^4.13)
`claim2_evalue_prefactor/` · rerun: `bash run_claim2.sh` (runs both binaries + `analyze_claim2.py`)

**Controlled design (version is the only variable):** same command `-usearch_local`, same input
FASTA (`inputs/v1_nodes_213.fasta`, query = db), same parameters
`-id 0.05 -evalue 1000 -maxaccepts 0 -maxrejects 0 -fulldp -wordlength 2`, same post-filter.

1. **Alignments are identical**, not merely similar — on 7,903 shared directed pairs the median
   difference (v12−v11) is **0.000 in percent identity, 0.000 in bit score, 0.000 in alignment
   length** (`01_alignments_identical.txt`). So the e-value difference is not a scoring/alignment
   difference.
2. **The e-value offset is flat across score** — median Δlog10(E) = **+4.129** at bits<100,
   100–200, 200–300, and ≥300 alike (`02_flat_offset.txt`). A λ / substitution-matrix difference
   would scale with score; flatness points to a multiplicative constant.
3. **Model selection localizes it to one constant** — both versions fit `E = K·m·N·2^(−bits)`
   (m = query length, N = 74,254 db letters), selecting the **same** `m·N` search-space form over
   8 competitors with **identical residual scatter σ=0.087** (`03_model_selection.txt`). Only the
   prefactor differs: **K_v11 = 1.01** (log10K +0.004) vs **K_v12 = 13,472** (log10K +4.132),
   IQR ±0.08 log across all ~7,900 edges.
4. **Minimal reproducible example** (`04_mre_n84_n35.txt`): pair **n84 vs n35**, identical
   alignment **30.4% id, alnlen 125, 54 bits** in both binaries, but **E_v11 = 1.13e-9 vs
   E_v12 = 1.53e-5** — ratio 13,540. Reproduce with a single command on 2 query sequences vs the
   provided db.

Impact: any criterion phrased in e-values (here, `e < 1e-5`) silently changes meaning between
versions. In our case the 10^4.13 inflation flips borderline edges and re-partitions the graph;
v11's e-value (K≈1, the textbook Karlin–Altschul form) is the correct reference.

Per-edge data for both issues: `claim2_evalue_prefactor/evalue_per_edge_v11_vs_v12.csv`
(`evalue_v11`, `evalue_v12`, `log10_Eratio…`, `impliedK_v11`, `impliedK_v12`, `pass_*`, `flips`).

## Layout
```
binaries.txt
claim1_missing_commands/  run_claim1.sh  claim1_output.txt  logs/
claim2_evalue_prefactor/  run_claim2.sh  analyze_claim2.py
                          pairs_u11.tsv pairs_u12.tsv  mre_u11.tsv mre_u12.tsv
                          01_alignments_identical.txt 02_flat_offset.txt
                          03_model_selection.txt 04_mre_n84_n35.txt
                          evalue_per_edge_v11_vs_v12.csv
inputs/  v1_nodes_213.fasta  mre_n84_n35.fasta
```
