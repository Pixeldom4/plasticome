#!/usr/bin/env python3
"""Step 0 (v5) - engine sanity: rebuild the 213 v1 nodes from the updated v1
reference CSV and emit both-orientation FASTA + node roster for usearch. All 213
rows are kept as DISTINCT nodes (md5-duplicates included), matching the paper's
Step 0 (213 rows, 205 unique md5, 8 redundant kept distinct). Target: 42 components.
"""
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
V1_IN = ROOT / "plasticome_v1.csv"
OUTDIR = ROOT / "step0_check"
OUT = OUTDIR / "v1_213"


def normalize(seq):
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def rep_acc(r):
    for k in ("pazy_accession", "uniprot_accession", "pdb_accession"):
        v = (r.get(k) or "").strip()
        if v and v != "-":
            return v
    return (r.get("gene") or r.get("rowid") or "").strip()


def main():
    OUTDIR.mkdir(exist_ok=True)
    rows = [r for r in csv.DictReader(V1_IN.open()) if normalize(r["aa_sequence"])]
    nodes = []
    for i, r in enumerate(sorted(rows, key=lambda r: int(r["rowid"])), start=1):
        seq = normalize(r["aa_sequence"])
        nodes.append({
            "node_id": f"v{i:03d}",
            "sequence": seq,
            "sequence_md5": hashlib.md5(seq.encode()).hexdigest(),
            "accession": rep_acc(r),
            "accession_all": rep_acc(r),
            "enzyme_name": (r.get("enzyme_name") or "").strip(),
            "organism": (r.get("host") or r.get("retrieved_host") or "").strip(),
            "sequence_length": len(seq),
            "v1_component": (r.get("component") or "").strip(),
            "cath": (r.get("cath") or "").strip(),
        })

    def write_fasta(path, ns):
        with Path(path).open("w") as fh:
            for n in ns:
                fh.write(f">{n['node_id']}\n")
                for i in range(0, len(n["sequence"]), 60):
                    fh.write(n["sequence"][i:i + 60] + "\n")

    write_fasta(f"{OUT}.fasta", nodes)
    write_fasta(f"{OUT}_rev.fasta", list(reversed(nodes)))
    cols = ["node_id", "sequence_md5", "accession", "accession_all", "enzyme_name",
            "organism", "sequence_length", "v1_component", "cath"]
    with Path(f"{OUT}_nodes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for n in nodes:
            w.writerow(n)
    print(f"wrote {len(nodes)} v1 nodes ({len({n['sequence_md5'] for n in nodes})} unique md5)")


if __name__ == "__main__":
    main()
