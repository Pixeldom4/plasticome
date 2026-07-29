"""
Run the exact-substring scan on the erickson-collapsed TSV.

Reuses the importable pipeline from `substring_scan.py` (normalization,
union-find grouping over exact containment, group -> row assembly); only the
record loader changes, because this TSV has its own schema:

    identifier, enzyme_name, pazy_id, accession, organism, aa_sequence

(one row per enzyme after the 'Enzyme N [like]' collapse). Accession is used as
the output `identifier`; sequences are normalized (whitespace-stripped,
uppercased) and de-duplicated on the normalized key, coalescing enzyme_name /
accession when identical sequences collapse.

Output: a TSV in the same 11-column schema as substring-scan_*.tsv, with the
manual-curation columns (change / comments / most-wildtype) left blank.
"""

import argparse
import os

import pandas as pd

from substring_scan import (
    norm, coalesce_join, build_groups, assemble_rows,
    exact_substring_match, exact_position,
)
from run_substring_scan_csv import TSV_COLUMNS, to_tsv_frame


def load_records_erickson(path):
    """Load unique-sequence records from the erickson-collapsed TSV."""
    df = pd.read_csv(path, sep='\t', dtype=str)
    recs = []
    for _, row in df.iterrows():
        key = norm(row['aa_sequence'])
        if not key:
            continue
        recs.append({
            'key': key,
            'seq': str(row['aa_sequence']).strip(),
            'len': len(key),
            'enzyme_name': coalesce_join(row.get('enzyme_name')),
            'accessions': coalesce_join(row.get('accession')),
            'source': coalesce_join(row.get('organism')),
        })
    out = pd.DataFrame(recs)
    # Collapse identical sequences, coalescing the identifying fields.
    agg = (out.groupby('key', as_index=False)
              .agg({'seq': 'first', 'len': 'first',
                    'enzyme_name': lambda s: coalesce_join(*s),
                    'accessions': lambda s: coalesce_join(*s),
                    'source': lambda s: coalesce_join(*s)}))
    return agg.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('tsv', help='input TSV (erickson-collapsed)')
    ap.add_argument('-o', '--out', help='output TSV (default: substring-scan_<stem>.tsv)')
    args = ap.parse_args()

    df = load_records_erickson(args.tsv)
    groups = build_groups(df, exact_substring_match)
    scan = assemble_rows(df, groups, exact_position)
    tsv = to_tsv_frame(scan)

    out = args.out or f"substring-scan_{os.path.splitext(os.path.basename(args.tsv))[0]}.tsv"
    tsv.to_csv(out, sep='\t', index=False)

    n_groups = scan['group'].nunique()
    print(f"records={len(df)} groups={n_groups} involved={len(scan)} "
          f"cores={(scan.role=='core').sum()} supersets={(scan.role=='superset').sum()}")
    print(f"wrote {out}")


if __name__ == '__main__':
    main()
