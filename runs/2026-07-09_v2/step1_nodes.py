#!/usr/bin/env python3
"""Step 1 - assemble the combined node set (205 v1-unique + 287 delta = 492).

Node identity is the recomputed md5 of the normalized sequence (letters only,
uppercased) - the frozen normalize() from reconciliation v2.

Canonical run uses pull sequences only (D-seq: pull-only). The plasticome.2
curation overlay is applied only in the sensitivity variant (--overlay).
"""
import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent
V1_TSV = DATA / "data" / "plasticome.v1.tsv"
CUR_TSV = DATA / "data" / "plasticome.2.tsv"
PULL_CSV = DATA / "archive" / "pazy_pull_2026-06-30.csv"


def normalize(seq: str) -> str:
    """Frozen normalize() from reconciliation v2: letters only, uppercased."""
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def md5_of(seq: str) -> str:
    return hashlib.md5(normalize(seq).encode()).hexdigest()


def rep_accession(row):
    """Representative accession. 33 pull records carry a sequence but no primary
    accession, which would make D4's 'smallest member accession' rule resolve on
    blanks; fall back through the other identifier columns, then the pazy_id."""
    for col in ("accession", "uniprot_accession", "pdb_accession"):
        v = (row.get(col) or "").strip()
        if v:
            return v
    allacc = (row.get("all_accessions") or "").strip()
    if allacc:
        return sorted(a for a in allacc.split(";") if a.strip())[0].strip()
    return f"pazy:{row['pazy_id'].strip()}"


def load_curation():
    rows = list(csv.DictReader(CUR_TSV.open(), delimiter="\t"))
    seq_col = [k for k in rows[0] if k.startswith("aa_sequence (")][0]
    by_id = {}
    for r in rows:
        pid = r["pazy_id"].strip()
        if pid:
            by_id[pid] = {
                "seq": normalize(r[seq_col]),
                "activity": (r["direct_or_indirect_activity"] or "").strip(),
            }
    return by_id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlay", action="store_true",
                    help="apply plasticome.2 curated sequences (sensitivity variant)")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    curation = load_curation()

    # ---- v1 nodes: 213 rows -> unique md5, carrying the curated component label
    v1_rows = list(csv.DictReader(V1_TSV.open(), delimiter="\t"))
    assert len(v1_rows) == 213, len(v1_rows)
    v1 = {}  # md5 -> node dict
    for r in v1_rows:
        seq = normalize(r["aa_sequence"])
        h = md5_of(seq)
        acc = (r["pazy_accession"] or r["retrieved"] or "").strip()
        if h not in v1:
            v1[h] = {
                "sequence_md5": h, "sequence": seq, "accession": acc,
                "sequence_length": len(seq), "v1_component": r["component"].strip(),
                "cath": r["cath"].strip(), "source": "v1", "is_new": 0,
                "pazy_id": "", "activity": "", "v1_rowids": [r["rowid"]],
                "accessions": {acc},
            }
        else:
            n = v1[h]
            n["v1_rowids"].append(r["rowid"])
            n["accessions"].add(acc)
            # gate: redundant v1 rows must agree on their component label
            assert n["v1_component"] == r["component"].strip(), (h, r["rowid"])

    # ---- PAZy pull: status column already encodes the delta
    pull = list(csv.DictReader(PULL_CSV.open()))
    delta, rescued, revised = {}, [], []
    for r in pull:
        pid = r["pazy_id"].strip()
        seq = normalize(r["aa_sequence"])
        cur = curation.get(pid, {})
        if args.overlay and cur.get("seq"):
            if not seq:
                rescued.append(pid)
            elif cur["seq"] != seq:
                revised.append(pid)
            seq = cur["seq"]
        if not seq:
            continue  # sequence-less: cannot be aligned
        h = md5_of(seq)
        if h in v1:
            continue  # md5 already a v1 node; nothing new
        if h in delta:
            delta[h]["accessions"].add(r["accession"].strip())
            continue
        delta[h] = {
            "sequence_md5": h, "sequence": seq, "accession": rep_accession(r),
            "sequence_length": len(seq), "v1_component": "", "cath": "",
            "source": "delta", "is_new": 1, "pazy_id": pid,
            "activity": cur.get("activity", ""), "v1_rowids": [],
            "accessions": {rep_accession(r)},
        }

    nodes = {**v1, **delta}
    # deterministic node ids: sorted by md5
    for i, h in enumerate(sorted(nodes), start=1):
        nodes[h]["node_id"] = f"m{i:04d}"

    ordered = [nodes[h] for h in sorted(nodes)]

    fa = Path(f"{args.out_prefix}.fasta")
    with fa.open("w") as fh:
        for n in ordered:
            fh.write(f">{n['node_id']}\n")
            for i in range(0, len(n["sequence"]), 60):
                fh.write(n["sequence"][i:i + 60] + "\n")

    tsv = Path(f"{args.out_prefix}_nodes.tsv")
    cols = ["node_id", "sequence_md5", "accession", "sequence_length", "source",
            "is_new", "v1_component", "cath", "pazy_id", "activity"]
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for n in ordered:
            w.writerow(n)

    stats = {
        "overlay": args.overlay, "n_nodes": len(nodes), "n_v1": len(v1),
        "n_delta": len(delta), "rescued_sequence_less": sorted(rescued, key=int),
        "revised_sequences": sorted(revised, key=int),
        "v1_redundant_rows": sum(len(n["v1_rowids"]) - 1 for n in v1.values()),
    }
    Path(f"{args.out_prefix}_stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
