#!/usr/bin/env python3
"""Step 1 (v4) - assemble the flat 476-node set from the merged, non-redundant
plasticome final CSV (enzyme_name, accession, organism, aa_sequence).

This is a FLAT de-novo run (decision 1): there is no v1 anchor and no delta.
Every row is one node; the file is already md5-unique under the frozen normalize()
(letters only, uppercased). v1 component/cath are attached only as an *overlay*
(interpretation for Step 3), never used to seed or gate the partition.

Representative-accession rule (G2): the lexicographically smallest accession in the
`;`-joined cell - deterministic. The full list is preserved in `accession_all`.

Outputs (petadex-alignment/v4/):
  combined.fasta / combined_rev.fasta   node FASTA in both orientations (B3)
  combined_nodes.tsv                    node roster
  combined_stats.json                   input provenance
"""
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # petadex-alignment/
V2 = ROOT / "v2"
CSV_IN = ROOT / "v4" / "plasticome.v2.cleaned_plasticome_final.csv"
OUT = ROOT / "v4" / "combined"


def normalize(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def md5_of(seq: str) -> str:
    return hashlib.md5(normalize(seq).encode()).hexdigest()


def rep_accession(cell: str) -> str:
    """Lexicographically smallest accession in a `;`-joined cell (G2)."""
    parts = [p.strip() for p in cell.split(";") if p.strip()]
    return min(parts) if parts else ""


def load_v1_overlay():
    """md5 -> (v1_component, cath) from the validated v1 reference (overlay only)."""
    overlay = {}
    for r in csv.DictReader((V2 / "canonical_nodes.tsv").open(), delimiter="\t"):
        if r["source"] != "v1":
            continue
        overlay[r["sequence_md5"]] = (r["v1_component"], r["cath"])
    return overlay


def main():
    overlay = load_v1_overlay()
    rows = list(csv.DictReader(CSV_IN.open()))

    nodes = {}
    for r in rows:
        seq = normalize(r["aa_sequence"])
        if not seq:
            continue
        h = md5_of(seq)
        assert h not in nodes, f"input not md5-unique: {h}"
        v1c, cath = overlay.get(h, ("", ""))
        nodes[h] = {
            "sequence_md5": h,
            "sequence": seq,
            "accession": rep_accession(r["accession"]),
            "accession_all": r["accession"].strip(),
            "enzyme_name": (r["enzyme_name"] or "").strip(),
            "organism": (r["organism"] or "").strip(),
            "sequence_length": len(seq),
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
            "organism", "sequence_length", "v1_component", "cath"]
    with Path(f"{OUT}_nodes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for n in ordered:
            w.writerow(n)

    lens = sorted(n["sequence_length"] for n in ordered)
    stats = {
        "input_csv": str(CSV_IN.relative_to(ROOT)),
        "input_md5": hashlib.md5(CSV_IN.read_bytes()).hexdigest(),
        "n_rows": len(rows),
        "n_nodes": len(nodes),
        "n_unique_md5": len(set(n["sequence_md5"] for n in ordered)),
        "db_letters": sum(lens),
        "len_min": lens[0], "len_median": lens[len(lens) // 2], "len_max": lens[-1],
        "n_lt_200aa": sum(l < 200 for l in lens),
        "n_lt_100aa": sum(l < 100 for l in lens),
        "n_multi_accession": sum(1 for n in ordered if ";" in n["accession_all"]),
        "n_with_v1_overlay": sum(1 for n in ordered if n["v1_component"]),
    }
    Path(f"{OUT}_stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
