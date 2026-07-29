# Curated Seed Set — All-vs-All Alignment

**Date:** 2026-07-23
**Input:** `petadex-clustering/curated-seed-validation/putative-curated-seed-set.csv`
**Method:** frozen PETadex component pipeline (usearch11 `-allpairs_local`, both orientations,
best-HSP by max-bits, post-filter **≥30% amino-acid identity AND e-value < 1e-5**).

This run validates the curated benchmark seed set — the 14 hand-picked reference enzymes
proposed as anchors for the plasticome clustering — by asking which of them are mutually
homologous under the paper-faithful similarity criterion.

---

## 1. Input

19 CSV rows → **14 sequences** (5 blank spacer rows dropped; every named enzyme retained),
14 unique md5s, DB = 4,738 residues.

| node | enzyme | accession | pazy_id | len | note (from `justification`) |
|------|--------|-----------|---------|-----|------|
| m0009 | ISPETase  | WP_054022242.1 | 1   | 290 | PET hydrolase reference |
| m0002 | LCC       | AEV21261.1     | 16  | 293 | PET hydrolase reference |
| m0008 | TfCut_2   | CBY05530       | 25  | 261 | long-standing cutinase reference |
| m0007 | Cut190    | BAO42836.1     | 18  | 304 | Ca²⁺-activated cutinase benchmark |
| m0005 | BhrPETase | GBD22443.1     | 17  | 293 | thermophilic, high-activity PETase |
| m0010 | HiC       | A0A075B5G4     | 44  | 194 | fungal cutinase outgroup |
| m0014 | Mipa-P    | SDZ16714.1     | 270 | 294 | **new** — putative unexplored PETase lineage |
| m0013 | Kubu-P    | WP_116180173.1 | 276 | 288 | **new** — putative unexplored PETase lineage |
| m0006 | NylA      | WP_012476897.1 | 59  | 493 | nylonase (distinct mechanism) |
| m0012 | NylB      | CAA24927.1     | 60  | 392 | nylonase |
| m0011 | NylC      | WP_012476895.1 | 61  | 355 | nylonase |
| m0003 | PhaZ1     | WP_011614907.1 | 100 | 419 | PHA depolymerase representative |
| m0004 | PLAase I  | WP_037320891.1 | 68  | 245 | PLA-ase representative |
| m0001 | PueA      | WP_011061486.1 | 50  | 617 | polyurethanase representative |

> Δ vs. the previous 12-seq run: adds **Mipa-P** and **Kubu-P**.

---

## 2. Summary

| metric | value |
|--------|-------|
| nodes | 14 |
| raw HSPs (both orientations) | 102 |
| unique aligned pairs | 51 |
| **passing edges** (≥30% id ∧ e<1e-5) | **21** |
| connected components | **8** (1 cluster + 7 singletons) |
| largest component | **7** |

**Headline:** all six PET-active references (ISPETase, LCC, BhrPETase, TfCut_2, Cut190) **plus
both new candidate lineages (Mipa-P, Kubu-P)** collapse into a single connected component
(**C002, n=7**). The two newcomers are confirmed homologous to the established PETase/cutinase
clade — Mipa-P and Kubu-P are each other's closest partner (56.1% id) and both link to every
established member above threshold. The nylonases, PHA depolymerase, PLA-ase and polyurethanase
stay as isolated singletons, as expected for their distinct folds/mechanisms.

---

## 3. Components

| component | n | members |
|-----------|---|---------|
| **C002** | **7** | LCC, BhrPETase, TfCut_2, Cut190, ISPETase, **Mipa-P**, **Kubu-P** |
| C001 | 1 | HiC |
| C003 | 1 | NylB |
| C004 | 1 | PueA |
| C005 | 1 | PhaZ1 |
| C006 | 1 | NylC |
| C007 | 1 | NylA |
| C008 | 1 | PLAase I |

---

## 4. All passing edges (21), sorted by % identity

| A | B | % id | e-value | bits |
|---|---|------|---------|------|
| LCC | BhrPETase | 92.2 | 2.4e-159 | 547 |
| Cut190 | TfCut_2 | 65.6 | 8.1e-103 | 359 |
| LCC | TfCut_2 | 59.9 | 8.1e-87 | 306 |
| BhrPETase | TfCut_2 | 59.9 | 5.3e-86 | 304 |
| **Kubu-P** | **Mipa-P** | 56.1 | 4.6e-91 | 320 |
| LCC | Cut190 | 54.4 | 4.7e-83 | 294 |
| BhrPETase | Cut190 | 53.3 | 1.6e-83 | 295 |
| TfCut_2 | ISPETase | 51.3 | 9.7e-72 | 256 |
| LCC | **Mipa-P** | 47.9 | 1.7e-69 | 249 |
| LCC | ISPETase | 47.3 | 5.0e-69 | 247 |
| BhrPETase | **Mipa-P** | 46.9 | 6.2e-67 | 240 |
| TfCut_2 | **Mipa-P** | 46.2 | 3.9e-65 | 234 |
| BhrPETase | ISPETase | 44.7 | 2.2e-64 | 232 |
| Cut190 | ISPETase | 44.5 | 1.8e-66 | 239 |
| TfCut_2 | **Kubu-P** | 44.1 | 1.6e-58 | 212 |
| BhrPETase | **Kubu-P** | 43.3 | 2.5e-60 | 218 |
| LCC | **Kubu-P** | 42.7 | 7.7e-62 | 223 |
| Cut190 | **Mipa-P** | 41.7 | 2.1e-62 | 225 |
| Cut190 | **Kubu-P** | 41.2 | 6.1e-59 | 214 |
| ISPETase | **Mipa-P** | 40.6 | 1.8e-58 | 212 |
| ISPETase | **Kubu-P** | 39.4 | 3.1e-55 | 201 |

- **LCC–BhrPETase (92.2%)** is the only pair ≥90% id — BhrPETase is essentially a
  thermostabilized LCC-family PETase. The next real hit drops to 65.6%.
- Mipa-P and Kubu-P join at **39–48% id** to the established members — clearly homologous,
  clearly a distinct sub-lineage (they are more similar to each other, 56%, than to any
  established reference). Weakest cluster edge: ISPETase–Kubu-P at 39.4%.

---

## 5. Caveat — HiC remains a singleton

HiC (fungal cutinase, 194 aa) is biologically a cutinase and would be *expected* inside C002,
but it falls out as a singleton. Its best hits are all fragmentary short local windows that fail
both thresholds:

| pair | % id | e-value | bits |
|------|------|---------|------|
| HiC–PLAase I | 58.8 | 5.6e+00 | 17 |
| HiC–NylC | 43.5 | 4.3e+00 | 18 |
| HiC–Kubu-P | 33.3 | 4.3e+00 | 18 |
| HiC–PueA | 25.9 | 9.6e+00 | 17 |

These are 17–18 bit, e≈4–10 fragments — high nominal %id over a handful of residues, **no**
significant global alignment. This is a known sensitivity limit of usearch's local seed-based
`-allpairs_local` for short, divergent sequences, **not** a data problem. If HiC needs to cluster
with the PET clade, use a profile/HMM search or a global aligner rather than usearch allpairs.
(The high-%id `HiC–PLAase I` row is meaningless without its bit score — a reminder that %id must
always be read alongside e-value/bits.)

---

## 6. Files (this folder)

| file | contents |
|------|----------|
| `results/combined_nodes.tsv` | 14-node roster (id, md5, enzyme, accession, length) |
| `results/combined_pairs.tsv` | raw usearch HSPs, forward orientation (headered) |
| `results/combined_rev_pairs.tsv` | raw usearch HSPs, reverse orientation (headered) |
| `results/component_edges_*.csv` | the 21 passing edges |
| `results/component_members_*.csv` | per-component membership |
| `results/component_assignment_*.csv` | per-node component assignment |
| `results/stats_*.json` | run provenance + summary counts |

Pairs-file columns: `query  target  pct_id  evalue  bits` (usearch `-userfields
query+target+id+evalue+bits`). Node IDs (`m0001`…) map to enzymes via `combined_nodes.tsv`.

## 7. Reproduce

```bash
cd petadex-alignment-usearch/analysis
# CSV -> TSV, then the frozen 3-step pipeline into outputs/curated_seed_putative/
python scripts/step1_nodes.py --tsv outputs/curated_seed_putative/input.tsv \
  --v1 outputs/curated_seed_putative/_empty_overlay.csv \
  --outprefix outputs/curated_seed_putative/combined
# usearch (Docker linux/amd64) on combined.fasta and combined_rev.fasta -> *_pairs.tsv
python scripts/step23_graph.py --prefix outputs/curated_seed_putative/combined \
  --outdir outputs/curated_seed_putative/results --tag curatedseed --date 2026-07-23
```

Method is frozen (v4–v8 paper-faithful); constants (`ID_MIN=30`, `EVALUE_MAX=1e-5`,
both-orientation search) live in `analysis/config.py` and must not change without a version bump.
