#!/usr/bin/env python3
"""Step 1 (v5) - assemble the md5-unique node set from the updated cleaned PAZy
final CSV (enzyme_name, accession, organism, aa_sequence).

Flat de-novo run (same policy as v4): there is no v1 anchor and no delta. Node
identity is the md5 of the frozen normalize() (letters only, uppercased), matching
the paper's "component" definition. v1 component/cath are attached only as an
*overlay* for Step 3 interpretation, never used to seed or gate the partition.

Differences from v4's input, handled here:
  * The updated CSV is NOT pre-deduplicated. Blank-sequence rows are DROPPED
    (cannot be aligned) and md5-duplicate rows are COLLAPSED to one node, with all
    contributing accessions / enzyme names preserved in *_all provenance columns.
    md5 = node identity is the paper-faithful reading (see workflow doc Step 1).
  * v1 overlay is read directly from the updated v1 reference CSV (component, cath),
    keyed by md5, instead of v2/canonical_nodes.tsv.

Representative-accession rule (G2): the lexicographically smallest accession among
the collapsed rows - deterministic. Full list preserved in `accession_all`.

Outputs (petadex-alignment/v5/):
  combined.fasta / combined_rev.fasta   node FASTA in both orientations (B3)
  combined_nodes.tsv                    node roster
  combined_stats.json                   input provenance
"""
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent                 # v5/
CSV_IN = ROOT / "cleaned_pazy_final.csv"
V1_IN = ROOT / "plasticome_v1.csv"
OUT = ROOT / "combined"


def normalize(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def md5_of(seq: str) -> str:
    return hashlib.md5(normalize(seq).encode()).hexdigest()


def load_v1_overlay():
    """md5 -> (v1_component, cath) from the updated v1 reference (overlay only)."""
    overlay = {}
    for r in csv.DictReader(V1_IN.open()):
        seq = normalize(r.get("aa_sequence", ""))
        if not seq:
            continue
        # first writer wins; v1 md5 map is de-facto unique on (component, cath)
        overlay.setdefault(md5_of(seq), (r.get("component", "").strip(),
                                         r.get("cath", "").strip()))
    return overlay


def main():
    overlay = load_v1_overlay()
    rows = list(csv.DictReader(CSV_IN.open()))

    # collapse to md5-unique nodes, preserving provenance
    groups = defaultdict(list)
    n_blank = 0
    for r in rows:
        seq = normalize(r["aa_sequence"])
        if not seq:
            n_blank += 1
            continue
        groups[md5_of(seq)].append(r)

    nodes = {}
    for h, mem in groups.items():
        seq = normalize(mem[0]["aa_sequence"])
        accs = sorted({(m["accession"] or "").strip() for m in mem if (m["accession"] or "").strip()})
        names = sorted({(m["enzyme_name"] or "").strip() for m in mem if (m["enzyme_name"] or "").strip()})
        orgs = sorted({(m["organism"] or "").strip() for m in mem if (m["organism"] or "").strip()})
        v1c, cath = overlay.get(h, ("", ""))
        nodes[h] = {
            "sequence_md5": h,
            "sequence": seq,
            "accession": accs[0] if accs else "",
            "accession_all": ";".join(accs),
            "enzyme_name": names[0] if names else "",
            "enzyme_name_all": ";".join(names),
            "organism": orgs[0] if orgs else "",
            "sequence_length": len(seq),
            "n_source_rows": len(mem),
            "v1_component": v1c,
            "cath": cath,
        }

    # deterministic node ids by md5 order
    for i, h in enumerate(sorted(nodes), start=1):
        nodes[h]["node_id"] = f"m{i:04d}"
    ordered = [nodes[h] for h in sorted(nodes)]

    def write_fasta(path, seq_nodes):
        with Path(path).open("w") as fh:
            for n in seq_nodes:
                fh.write(f">{n['node_id']}\n")
                for i in range(0, len(n["sequence"]), 60):
                    fh.write(n["sequence"][i:i + 60] + "\n")

    write_fasta(f"{OUT}.fasta", ordered)
    write_fasta(f"{OUT}_rev.fasta", list(reversed(ordered)))  # B3 orientation symmetrization

    cols = ["node_id", "sequence_md5", "accession", "accession_all", "enzyme_name",
            "enzyme_name_all", "organism", "sequence_length", "n_source_rows",
            "v1_component", "cath"]
    with Path(f"{OUT}_nodes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for n in ordered:
            w.writerow(n)

    lens = sorted(n["sequence_length"] for n in ordered)
    stats = {
        "input_csv": str(CSV_IN.relative_to(ROOT)),
        "input_md5": hashlib.md5(CSV_IN.read_bytes()).hexdigest(),
        "v1_overlay_csv": str(V1_IN.relative_to(ROOT)),
        "n_rows": len(rows),
        "n_blank_seq_dropped": n_blank,
        "n_rows_with_seq": len(rows) - n_blank,
        "n_nodes": len(nodes),
        "n_md5_collapsed_rows": (len(rows) - n_blank) - len(nodes),
        "n_multi_row_nodes": sum(1 for n in ordered if n["n_source_rows"] > 1),
        "n_unique_md5": len(set(n["sequence_md5"] for n in ordered)),
        "db_letters": sum(lens),
        "len_min": lens[0], "len_median": lens[len(lens) // 2], "len_max": lens[-1],
        "n_lt_200aa": sum(l < 200 for l in lens),
        "n_lt_100aa": sum(l < 100 for l in lens),
        "n_with_v1_overlay": sum(1 for n in ordered if n["v1_component"]),
    }
    Path(f"{OUT}_stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
