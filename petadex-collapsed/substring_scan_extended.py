"""
Extended: near-substring scan.

Same pipeline as the base exact scan (`substring_scan.py`) — same loading,
union-find grouping, and sheet output — but the pairwise test is relaxed so it
ALSO catches minor variants that exact containment misses: point substitutions
and internal indels (e.g. an engineered mutant vs its wildtype). Writes to a
separate sheet ('near-substring-scan') so the exact scan's results are untouched.

How the relaxed test works (`near_substring_match`):
  difflib matching blocks between the shorter and longer key are summed; if they
  cover at least MIN_COVERAGE of the shorter sequence, the shorter is treated as
  a near-substring of the longer. A handful of substitutions/indels break the
  match into a few blocks that still cover most of the shorter sequence, so they
  pass; unrelated sequences share little and fall well below the threshold.

Notes / caveats:
  - This is a heuristic, not a real alignment. Tune MIN_COVERAGE (or swap in a
    Biopython pairwise aligner) if you want stricter/looser calls.
  - It is O(n^2) difflib comparisons, heavier than the exact `in` test, but fine
    at the scale of unique sequences in this dataset.
  - Exact containments still get precise prefix/suffix/internal labels; genuinely
    near (non-exact) matches are labelled 'near'.

Run:  python substring_scan_extended.py
"""

import difflib

from substring_scan import run, exact_position, load_config


def make_near_match(min_coverage):
    """Build a near-substring predicate for the given coverage threshold.

    Exact containment short-circuits to True; otherwise the difflib matching
    blocks between the two keys must cover at least `min_coverage` of the
    shorter sequence.
    """
    def near_substring_match(short_key, long_key):
        if short_key in long_key:
            return True
        sm = difflib.SequenceMatcher(None, short_key, long_key, autojunk=False)
        matched = sum(block.size for block in sm.get_matching_blocks())
        return matched / len(short_key) >= min_coverage
    return near_substring_match


def near_position(core_key, s):
    """Precise label when the core is exactly contained; else 'near'."""
    return exact_position(core_key, s) if core_key in s else 'near'


if __name__ == '__main__':
    cfg = load_config()
    run(match=make_near_match(cfg['near']['min_coverage']),
        position_fn=near_position,
        sheet_name=cfg['near']['sheet_name'],
        config=cfg)
