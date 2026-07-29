"""
Run the exact-substring scan on a CSV export (e.g. the cleaned pazy
retrieving-from-accession table) instead of the joined Excel workbook.

Reuses the importable pipeline from `substring_scan.py` (normalization,
union-find grouping over exact containment, and group -> row assembly); only
the record loader changes, because this CSV has a different schema than the
`joined` sheet:

    enzyme_name, accession, ..., aa_sequence, ..., sequence_source, ...

Sequences are normalized (whitespace-stripped, uppercased) and de-duplicated on
the normalized key; when identical sequences collapse, their enzyme_name and
accession values are coalesced into ';'-joined lists (same convention as the
workbook loader).

Output: a TSV in the same 11-column schema as substring-scan_result.tsv, with
the manual-curation columns (change / comments / most-wildtype) left blank.
"""

import argparse
import os

import pandas as pd

from substring_scan import (
    norm, coalesce_join, build_groups, assemble_rows,
    exact_substring_match, exact_position,
)


def load_records_csv(path):
    """Load unique-sequence records from the cleaned pazy CSV."""
    df = pd.read_csv(path)
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
            'source': row.get('sequence_source'),
        })
    out = pd.DataFrame(recs)
    # Collapse identical sequences, coalescing the identifying fields.
    agg = (out.groupby('key', as_index=False)
              .agg({'seq': 'first', 'len': 'first',
                    'enzyme_name': lambda s: coalesce_join(*s),
                    'accessions': lambda s: coalesce_join(*s),
                    'source': lambda s: coalesce_join(*s)}))
    return agg.reset_index(drop=True)


# Output schema mirrors substring-scan_result.tsv.
TSV_COLUMNS = ['group', 'identifier', 'enzyme_name', 'length_aa', 'role',
               'delta_aa', 'position', 'change', 'aa_sequence',
               'comments', 'most-wildtype']


def to_tsv_frame(scan):
    """Re-shape assemble_rows() output into the deliverable TSV schema."""
    out = pd.DataFrame({
        'group': scan['group'],
        'identifier': scan['accessions'],
        'enzyme_name': scan['enzyme_name'],
        'length_aa': scan['length_aa'],
        'role': scan['role'],
        'delta_aa': scan['delta_aa'].apply(lambda d: f'+{d}' if d > 0 else str(d)),
        'position': scan['core_position'],
        'change': '',
        'aa_sequence': scan['aa_sequence'],
        'comments': '',
        'most-wildtype': '',
    })
    return out[TSV_COLUMNS]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('csv', help='input CSV (cleaned pazy retrieving-from-accession)')
    ap.add_argument('-o', '--out', help='output TSV (default: substring-scan_<csv-stem>.tsv)')
    args = ap.parse_args()

    df = load_records_csv(args.csv)
    groups = build_groups(df, exact_substring_match)
    scan = assemble_rows(df, groups, exact_position)
    tsv = to_tsv_frame(scan)

    out = args.out or f"substring-scan_{os.path.splitext(os.path.basename(args.csv))[0]}.tsv"
    tsv.to_csv(out, sep='\t', index=False)

    n_groups = scan['group'].nunique()
    print(f"records={len(df)} groups={n_groups} involved={len(scan)} "
          f"cores={(scan.role=='core').sum()} supersets={(scan.role=='superset').sum()}")
    print(f"wrote {out}")


if __name__ == '__main__':
    main()
