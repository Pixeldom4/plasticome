"""Run the full validation pipeline: Tier 0 -> 1 -> 2 -> 3.

All NCBI responses (successes and failures alike) are cached under .cache/, so a
second run is near-instant. Delete .cache/ to re-query NCBI from scratch.

    python3 run_all.py

Set NCBI_ENV_VAR to an NCBI API key to lift the request ceiling from 3/s to 10/s.
"""

from __future__ import annotations

import sys
import time

import tier0_resolve
import tier1_names
import tier2_doi
import tier3_report
import tierA_remap
import tierO_organism
import tierP_provenance
import tierS_sequence

STAGES = [
    ("Tier 0  existence and format", tier0_resolve.main),
    ("Tier 1  name concordance", tier1_names.main),
    ("Tier 2  DOI cross-reference", tier2_doi.main),
    ("Tier O  organism concordance", tierO_organism.main),
    ("Tier S  sequence identity to cited paper", tierS_sequence.main),
    ("Tier A  sequence-driven re-mapping (BLAST + author tiebreaker)", tierA_remap.main),
    ("Tier P  protein<->nucleotide provenance", tierP_provenance.main),
    ("Tier 3  adjudication queue", tier3_report.main),
]


def main() -> int:
    for title, fn in STAGES:
        print(f"\n{'=' * 64}\n{title}\n{'=' * 64}")
        start = time.monotonic()
        rc = fn()
        if rc != 0:
            print(f"\n{title} failed (rc={rc})", file=sys.stderr)
            return rc
        print(f"-- {time.monotonic() - start:.1f}s")
    print("\nDone. See out/review_queue.csv and out/validation_full.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
