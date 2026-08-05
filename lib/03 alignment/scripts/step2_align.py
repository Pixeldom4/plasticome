#!/usr/bin/env python3
"""Step 2 - all-vs-all alignment of the node FASTA, permissive at search time.

The one place the pipeline's aligner is chosen. Everything about *what counts as an
edge* lives in step23_graph.py (>=30% aaid AND e<1e-5); this step only has to produce
`query, target, pct_id, evalue, bits` rows and drop nothing that might pass that gate.

  diamond (default)  one `blastp` of the FASTA against itself. Query == database, and
                     DIAMOND reports each pair in both directions inside that single
                     run, so no reversed-orientation second pass is needed.
  usearch            `-allpairs_local -acceptall` per orientation. -allpairs_local is
                     asymmetric in which sequence it treats as query, so the reversed
                     FASTA is searched too and the best HSP per unordered pair is taken
                     by max bits downstream (the B3 convention; keep both together).

Writes <prefix>_pairs.tsv, plus <prefix>_rev_pairs.tsv under usearch.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True, type=Path,
                    help="node-set stem from step1_nodes.py (expects <prefix>.fasta)")
    ap.add_argument("--engine", default=config.ALIGNER, choices=["diamond", "usearch"])
    ap.add_argument("--log", type=Path, default=None, help="aligner log (default: beside the pairs)")
    args = ap.parse_args()

    fasta = Path(f"{args.prefix}.fasta")
    if not fasta.exists():
        sys.exit(f"error: no node FASTA at {fasta} (run step1_nodes.py first)")

    if args.engine == "diamond":
        log = args.log or fasta.parent / "diamond.log"
        out = config.diamond(fasta, Path(f"{args.prefix}_pairs.tsv"), log)
        print(f"[step2:diamond] {fasta.name} -> {out.name}")
    else:
        rev = Path(f"{args.prefix}_rev.fasta")
        if not rev.exists():
            sys.exit(f"error: usearch needs both orientations; missing {rev}")
        config.usearch(fasta, Path(f"{args.prefix}_pairs.tsv"))
        config.usearch(rev, Path(f"{args.prefix}_rev_pairs.tsv"))
        print(f"[step2:usearch] {fasta.name} + {rev.name} -> {args.prefix.name}_pairs.tsv (x2)")


if __name__ == "__main__":
    main()
