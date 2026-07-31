#!/usr/bin/env python3
"""Step 0 - engine sanity. Rebuild the 213 v1 reference nodes and emit
both-orientation FASTA + a node roster for usearch.

All 213 v1 rows are kept as DISTINCT nodes (md5-duplicates included), matching the
paper's Step 0 (213 rows, 205 unique md5, 8 redundant kept distinct). Running the
identical Step 2/3 engine over these must land on the paper's **42 components /
3,178 edges** -- that equality is the whole point of this step. It validates the
aligner + partition, nothing biological about the 484-node set.

The roster carries the SAME columns as step1_nodes.py (identifier/identifier_all/
pazy_id present but blank for v1 nodes, accession populated) so the one unified
step23_graph.py runs on it unchanged. Component COUNT and edge COUNT are label-
independent, so the accession fallback labeling here still reproduces the 42.
"""
import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402


def normalize(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def rep_acc(r: dict) -> str:
    for k in ("pazy_accession", "uniprot_accession", "pdb_accession"):
        v = (r.get(k) or "").strip()
        if v and v != "-":
            return v
    return (r.get("gene") or r.get("rowid") or "").strip()


def write_fasta(path: Path, ns: list) -> None:
    with Path(path).open("w") as fh:
        for n in ns:
            fh.write(f">{n['node_id']}\n")
            for i in range(0, len(n["sequence"]), 60):
                fh.write(n["sequence"][i:i + 60] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--v1", type=Path, default=config.V1_CSV)
    ap.add_argument("--outdir", type=Path, default=config.STEP0_DIR)
    ap.add_argument("--prefix", default="v1_213",
                    help="stem for the emitted fasta/roster (default v1_213)")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / args.prefix

    rows = [r for r in csv.DictReader(args.v1.open()) if normalize(r["aa_sequence"])]
    nodes = []
    for i, r in enumerate(sorted(rows, key=lambda r: int(r["rowid"])), start=1):
        seq = normalize(r["aa_sequence"])
        acc = rep_acc(r)
        nodes.append({
            "node_id": f"v{i:03d}",
            "sequence": seq,
            "sequence_md5": hashlib.md5(seq.encode()).hexdigest(),
            # identifier/pazy_id blank for v1 nodes -> graph script falls back to accession
            "identifier": "", "identifier_all": "", "pazy_id": "", "pazy_id_all": "",
            "accession": acc, "accession_all": acc,
            "enzyme_name": (r.get("enzyme_name") or "").strip(),
            "enzyme_name_all": (r.get("enzyme_name") or "").strip(),
            "organism": (r.get("host") or r.get("retrieved_host") or "").strip(),
            "sequence_length": len(seq),
            "n_source_rows": 1,
            "v1_component": (r.get("component") or "").strip(),
            "cath": (r.get("cath") or "").strip(),
        })

    write_fasta(f"{out}.fasta", nodes)
    write_fasta(f"{out}_rev.fasta", list(reversed(nodes)))  # orientation symmetrization

    cols = ["node_id", "sequence_md5", "identifier", "identifier_all", "accession",
            "accession_all", "enzyme_name", "enzyme_name_all", "organism",
            "pazy_id", "pazy_id_all", "sequence_length", "n_source_rows",
            "v1_component", "cath"]
    with Path(f"{out}_nodes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(nodes)

    n_md5 = len({n["sequence_md5"] for n in nodes})
    print(f"[step0] wrote {len(nodes)} v1 nodes ({n_md5} unique md5) -> {out}.fasta")


if __name__ == "__main__":
    main()
