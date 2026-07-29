#!/usr/bin/env python3
"""Step 1 - assemble the md5-unique node set from the curated PAZy final TSV.

Flat de-novo run (same policy as v4-v8): no v1 anchor, no delta. Node identity is
the md5 of the frozen normalize() (letters only, uppercased), matching the paper's
"component" definition. The input's own `component_id`/`cath` and the v1 overlay are
attached only as descriptive columns for Step 3 interpretation -- never used to seed
or gate the partition.

  * Every row carries a stable unique `identifier` (PL<n>) and a `pazy_id`, so the
    canonical LABEL is the PL identifier.
  * md5-duplicate rows are COLLAPSED to one node; every contributing identifier /
    accession / enzyme_name / pazy_id is preserved in the *_all provenance columns.
  * Representative rule (G2): among collapsed rows, the smallest PL number wins
    (deterministic); full lists kept in *_all.

Outputs (to --outprefix, default outputs/combined):
  combined.fasta / combined_rev.fasta   node FASTA, both orientations (B3)
  combined_nodes.tsv                    node roster
  combined_stats.json                   input provenance + length stats
"""
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402


def normalize(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def md5_of(seq: str) -> str:
    return hashlib.md5(normalize(seq).encode()).hexdigest()


def pl_num(ident: str) -> int:
    m = re.search(r"(\d+)", ident or "")
    return int(m.group(1)) if m else 10**9


def load_v1_overlay(v1_csv: Path) -> dict:
    """md5 -> (v1_component, cath) from the v1 reference (overlay only)."""
    overlay = {}
    for r in csv.DictReader(v1_csv.open()):
        seq = normalize(r.get("aa_sequence", ""))
        if not seq:
            continue
        overlay.setdefault(md5_of(seq),
                           (r.get("component", "").strip(), r.get("cath", "").strip()))
    return overlay


def uniq_sorted(mem, col, key=None):
    return sorted({(m.get(col) or "").strip() for m in mem if (m.get(col) or "").strip()},
                  key=key)


def write_fasta(path: Path, seq_nodes: list) -> None:
    with Path(path).open("w") as fh:
        for n in seq_nodes:
            fh.write(f">{n['node_id']}\n")
            for i in range(0, len(n["sequence"]), 60):
                fh.write(n["sequence"][i:i + 60] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", type=Path, default=config.INPUT_TSV)
    ap.add_argument("--v1", type=Path, default=config.V1_CSV)
    ap.add_argument("--outprefix", type=Path, default=config.OUTPUTS / "combined")
    args = ap.parse_args()
    args.outprefix.parent.mkdir(parents=True, exist_ok=True)
    out = args.outprefix

    overlay = load_v1_overlay(args.v1)
    rows = list(csv.DictReader(args.tsv.open(), delimiter="\t"))

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
        mem = sorted(mem, key=lambda m: pl_num(m.get("identifier", "")))
        idents = uniq_sorted(mem, "identifier", key=pl_num)
        accs = uniq_sorted(mem, "accession")
        names = uniq_sorted(mem, "enzyme_name")
        orgs = uniq_sorted(mem, "organism")
        pazy = uniq_sorted(mem, "pazy_id")
        v1c, cath = overlay.get(h, ("", ""))
        nodes[h] = {
            "sequence_md5": h, "sequence": seq,
            "identifier": idents[0] if idents else "", "identifier_all": ";".join(idents),
            "accession": accs[0] if accs else "", "accession_all": ";".join(accs),
            "enzyme_name": names[0] if names else "", "enzyme_name_all": ";".join(names),
            "organism": orgs[0] if orgs else "",
            "pazy_id": pazy[0] if pazy else "", "pazy_id_all": ";".join(pazy),
            "sequence_length": len(seq), "n_source_rows": len(mem),
            "v1_component": v1c, "cath": cath,
        }

    # deterministic node ids by md5 order
    for i, h in enumerate(sorted(nodes), start=1):
        nodes[h]["node_id"] = f"m{i:04d}"
    ordered = [nodes[h] for h in sorted(nodes)]

    write_fasta(f"{out}.fasta", ordered)
    write_fasta(f"{out}_rev.fasta", list(reversed(ordered)))  # B3 orientation symmetrization

    cols = ["node_id", "sequence_md5", "identifier", "identifier_all", "accession",
            "accession_all", "enzyme_name", "enzyme_name_all", "organism",
            "pazy_id", "pazy_id_all", "sequence_length", "n_source_rows",
            "v1_component", "cath"]
    with Path(f"{out}_nodes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(ordered)

    lens = sorted(n["sequence_length"] for n in ordered)
    stats = {
        "input_tsv": str(args.tsv), "input_md5": hashlib.md5(args.tsv.read_bytes()).hexdigest(),
        "v1_overlay_csv": str(args.v1),
        "n_rows": len(rows), "n_blank_seq_dropped": n_blank,
        "n_rows_with_seq": len(rows) - n_blank, "n_nodes": len(nodes),
        "n_md5_collapsed_rows": (len(rows) - n_blank) - len(nodes),
        "n_multi_row_nodes": sum(1 for n in ordered if n["n_source_rows"] > 1),
        "n_unique_md5": len({n["sequence_md5"] for n in ordered}),
        "db_letters": sum(lens),
        "len_min": lens[0], "len_median": lens[len(lens) // 2], "len_max": lens[-1],
        "n_lt_200aa": sum(l < 200 for l in lens), "n_lt_100aa": sum(l < 100 for l in lens),
        "n_with_v1_overlay": sum(1 for n in ordered if n["v1_component"]),
    }
    Path(f"{out}_stats.json").write_text(json.dumps(stats, indent=2))
    print(f"[step1] {len(nodes)} nodes / {stats['n_unique_md5']} unique md5 "
          f"({n_blank} blank, {stats['n_md5_collapsed_rows']} md5-collapsed) -> {out}.fasta")


if __name__ == "__main__":
    main()
