#!/usr/bin/env python3
"""Step 1 (v3) - assemble the combined node set from the v1 reference graph plus a
fresh sequence file (pazy-260701_cleaned.csv: enzyme_name, organism, aa_sequence).

Node identity is the recomputed md5 of the normalized sequence (letters only,
uppercased) - the frozen normalize() from reconciliation v2.

Reference frame = the validated v1 42-component seed graph (205 unique md5 nodes),
reconstructed from v2/canonical_nodes.tsv + v2/canonical.fasta (the artifacts, since
plasticome.v1.tsv is not in data/). These anchor the partition and carry v1_component
+ cath as validation labels. The file's sequences are unioned in; any that are md5-
identical to a v1 node inherit that node (present_in_file=1), the rest are new nodes
(source="file", is_new=1) to be classified by the recomputed graph.

Output format matches step1_nodes.py exactly so v2/step23_graph.py runs unchanged.
"""
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # petadex-alignment/
V2 = ROOT / "v2"
FILE = ROOT / "data" / "pazy-260701_cleaned.csv"
OUT_PREFIX = ROOT / "v3" / "combined"


def normalize(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def md5_of(seq: str) -> str:
    return hashlib.md5(normalize(seq).encode()).hexdigest()


def load_v1_reference():
    """205 v1-unique md5 nodes with v1_component + cath, sequences from canonical.fasta."""
    fa, cur = {}, None
    for l in (V2 / "canonical.fasta").open():
        if l.startswith(">"):
            cur = l[1:].strip(); fa[cur] = []
        else:
            fa[cur].append(l.strip())
    fa = {k: "".join(v) for k, v in fa.items()}

    v1 = {}
    for r in csv.DictReader((V2 / "canonical_nodes.tsv").open(), delimiter="\t"):
        if r["source"] != "v1":
            continue
        seq = fa[r["node_id"]]
        h = r["sequence_md5"]
        assert md5_of(seq) == h, r["node_id"]
        v1[h] = {
            "sequence_md5": h, "sequence": seq, "accession": r["accession"],
            "sequence_length": len(seq), "v1_component": r["v1_component"],
            "cath": r["cath"], "source": "v1", "is_new": 0,
            "pazy_id": "", "activity": "",
        }
    assert len(v1) == 205, len(v1)
    return v1


def main():
    v1 = load_v1_reference()

    # ---- file: enzyme_name, organism, aa_sequence
    rows = list(csv.DictReader(FILE.open()))
    file_md5_rows = defaultdict(list)   # md5 -> list of (enzyme_name, organism)
    for r in rows:
        seq = normalize(r["aa_sequence"])
        if not seq:
            continue
        file_md5_rows[md5_of(seq)].append(
            ((r["enzyme_name"] or "").strip(), (r["organism"] or "").strip(), seq))

    delta = {}
    in_file_v1 = 0
    for h, recs in file_md5_rows.items():
        if h in v1:
            in_file_v1 += 1
            continue                     # already an anchor node; inherits its label
        names = sorted(n for n, _o, _s in recs if n)
        acc = names[0] if names else f"file:{h[:8]}"
        seq = recs[0][2]
        delta[h] = {
            "sequence_md5": h, "sequence": seq, "accession": acc,
            "sequence_length": len(seq), "v1_component": "", "cath": "",
            "source": "delta", "is_new": 1, "pazy_id": "",
            "activity": "",
        }

    nodes = {**v1, **delta}
    for i, h in enumerate(sorted(nodes), start=1):
        nodes[h]["node_id"] = f"m{i:04d}"
    ordered = [nodes[h] for h in sorted(nodes)]

    fa = Path(f"{OUT_PREFIX}.fasta")
    with fa.open("w") as fh:
        for n in ordered:
            fh.write(f">{n['node_id']}\n")
            for i in range(0, len(n["sequence"]), 60):
                fh.write(n["sequence"][i:i + 60] + "\n")

    # reverse-ordered fasta (B3: orientation symmetrization)
    with Path(f"{OUT_PREFIX}_rev.fasta").open("w") as fh:
        for n in reversed(ordered):
            fh.write(f">{n['node_id']}\n")
            for i in range(0, len(n["sequence"]), 60):
                fh.write(n["sequence"][i:i + 60] + "\n")

    cols = ["node_id", "sequence_md5", "accession", "sequence_length", "source",
            "is_new", "v1_component", "cath", "pazy_id", "activity"]
    with Path(f"{OUT_PREFIX}_nodes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for n in ordered:
            w.writerow(n)

    # md5 -> node_id map + full file-row roster for mapping assignments back to rows
    with Path(f"{OUT_PREFIX}_file_rows.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["enzyme_name", "organism", "sequence_md5", "node_id",
                    "sequence_length", "matches_v1_node"])
        for h, recs in sorted(file_md5_rows.items()):
            nid = nodes[h]["node_id"]
            for name, org, seq in recs:
                w.writerow([name, org, h, nid, len(seq), int(h in v1)])

    stats = {
        "n_nodes": len(nodes), "n_v1_reference": len(v1),
        "n_file_unique_md5": len(file_md5_rows),
        "n_file_rows": len(rows),
        "file_seqs_matching_v1_node": in_file_v1,
        "n_new_nodes": len(delta),
        "db_letters": sum(len(n["sequence"]) for n in ordered),
    }
    Path(f"{OUT_PREFIX}_stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
