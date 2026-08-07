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


def acc_key(acc: str) -> str:
    """Accession as an overlay key: version suffix stripped, upper-cased."""
    return (acc or "").strip().split(".")[0].upper()


def load_seq_by_accession(path) -> dict:
    """accession (version-stripped) -> aa_sequence, from a union-shaped TSV.

    Lets the overlay key on sequence md5 even when the v1 table itself no longer
    carries sequences: the union holds the sequence that was retrieved for each
    v1 accession, which is exactly what the md5 key needs.
    """
    seqs = {}
    if not path:
        return seqs
    path = Path(path)
    if not path.exists():
        print(f"[step1] no overlay sequence source at {path}")
        return seqs
    for r in csv.DictReader(path.open(), delimiter="\t"):
        a, s = acc_key(r.get("accession", "")), (r.get("aa_sequence") or "").strip()
        if a and s:
            seqs.setdefault(a, s)
    return seqs


def load_v1_overlay(v1_csv: Path, seq_source=None):
    """-> (md5 -> (v1_component, cath, reference), accession -> same).

    Two keys because neither alone is sufficient. The md5 of the sequence is
    exact, but it only works while the v1 table still carries sequences: the
    `.no-seq` re-export keeps 13 of 212, which cut the overlay from 161 nodes to
    10 without erroring. `retrieved` is the same accession the union joins on, so
    it covers every row that has one, sequence or not.

    md5 is tried first and accession second, so an exact sequence match always
    outranks an accession match on the same node.

    A missing file yields an empty overlay -- v1_component/cath/reference stay
    blank and nothing else changes, since the partition never reads them. The
    delimiter follows the extension, so either the v1 CSV or the tab-separated
    v1.1 table can be supplied.

    `reference` is the v1 table's literature citation, carried purely as
    annotation. Rows are still admitted on component/cath alone: on the v1.1
    table every one of the 211 rows with a reference also has a component, so
    gating on the citation as well would admit no extra row while risking a
    blank-overlay key shadowing a real hit.
    """
    by_md5, by_acc = {}, {}
    v1_csv = Path(v1_csv)
    if not v1_csv.exists():
        print(f"[step1] no v1 overlay at {v1_csv} -- v1_component/cath/reference left blank")
        return by_md5, by_acc
    borrowed = load_seq_by_accession(seq_source)
    n_borrowed = 0
    delim = "\t" if v1_csv.suffix.lower() in (".tsv", ".tab") else ","
    for r in csv.DictReader(v1_csv.open(), delimiter=delim):
        val = (r.get("component", "").strip(), r.get("cath", "").strip(),
               r.get("reference", "").strip())
        if not any(val[:2]):
            continue  # nothing to overlay; keying on it would only shadow a real hit
        acc = acc_key(r.get("retrieved", ""))
        seq = normalize(r.get("aa_sequence", ""))
        if not seq and acc in borrowed:
            seq = normalize(borrowed[acc])
            n_borrowed += 1
        if seq:
            by_md5.setdefault(md5_of(seq), val)
        if acc:
            by_acc.setdefault(acc, val)
    print(f"[step1] v1 overlay: {len(by_md5)} md5 keys "
          f"({n_borrowed} from the sequence source), {len(by_acc)} accession keys")
    return by_md5, by_acc


def load_doi_overlay(v2_tsv):
    """-> (pazy_id -> doi, accession -> doi) from the v260701 table.

    The DOI is the v260701 side of the same literature annotation the v1 table
    spells as `reference`, and it is the better-covered of the two: every one of
    the 473 rows carries one, against 167 of 413 clusters for the v1 overlay.
    It is read here rather than through the union because `01-union.tsv` has no
    doi column, and adding one would mean re-running step 1's NCBI fetches.

    pazy_id is the primary key and accession only the fallback. In the union,
    pazy_id is populated exactly on the rows that came from v260701 (468 of them,
    all distinct, blank on every v1.1-only row), so it cannot mis-join a v1.1 row
    onto a v260701 citation. Accession is looser -- 445 rows share 442 distinct
    version-stripped values -- so it is consulted only when pazy_id misses.

    A missing file yields an empty overlay and `doi` stays blank, exactly as the
    v1 overlay behaves; the partition never reads either.
    """
    by_pazy, by_acc = {}, {}
    if v2_tsv is None:
        return by_pazy, by_acc
    v2_tsv = Path(v2_tsv)
    if not v2_tsv.exists():
        print(f"[step1] no v260701 table at {v2_tsv} -- doi left blank")
        return by_pazy, by_acc
    delim = "\t" if v2_tsv.suffix.lower() in (".tsv", ".tab") else ","
    for r in csv.DictReader(v2_tsv.open(), delimiter=delim):
        doi = (r.get("doi") or "").strip()
        if not doi:
            continue
        pazy = (r.get("pazy_id") or "").strip()
        acc = acc_key(r.get("accession", ""))
        if pazy:
            by_pazy.setdefault(pazy, doi)
        if acc:
            by_acc.setdefault(acc, doi)
    print(f"[step1] doi overlay: {len(by_pazy)} pazy_id keys, {len(by_acc)} accession keys")
    return by_pazy, by_acc


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
    ap.add_argument("--v1-seqs", type=Path, default=None,
                    help="union TSV supplying sequences for v1 rows that carry none, "
                         "so the overlay can still key on sequence md5")
    ap.add_argument("--v2", type=Path, default=None,
                    help="v260701 TSV supplying the doi citation, keyed on pazy_id "
                         "then accession")
    ap.add_argument("--outprefix", type=Path, default=config.OUTPUTS / "combined")
    args = ap.parse_args()
    args.outprefix.parent.mkdir(parents=True, exist_ok=True)
    out = args.outprefix

    ov_md5, ov_acc = load_v1_overlay(args.v1, args.v1_seqs)
    doi_pazy, doi_acc = load_doi_overlay(args.v2)
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
        # md5 first (exact sequence), then any of the node's accessions
        v1c, cath, ref = ov_md5.get(h, ("", "", ""))
        via = "md5" if (v1c or cath) else ""
        if not via:
            for a in accs:
                hit = ov_acc.get(acc_key(a))
                if hit:
                    v1c, cath, ref = hit
                    via = "accession"
                    break
        # A node is md5-unique and can collapse several source rows, so every one
        # of its pazy_ids is consulted; accessions are the fallback for the rows
        # v260701 left without a pazy_id.
        dois = [doi_pazy[p] for p in pazy if p in doi_pazy]
        if not dois:
            dois = [doi_acc[k] for k in (acc_key(a) for a in accs) if k in doi_acc]
        nodes[h] = {
            "v1_overlay_via": via,
            "sequence_md5": h, "sequence": seq,
            "identifier": idents[0] if idents else "", "identifier_all": ";".join(idents),
            "accession": accs[0] if accs else "", "accession_all": ";".join(accs),
            "enzyme_name": names[0] if names else "", "enzyme_name_all": ";".join(names),
            "organism": orgs[0] if orgs else "",
            "pazy_id": pazy[0] if pazy else "", "pazy_id_all": ";".join(pazy),
            "sequence_length": len(seq), "n_source_rows": len(mem),
            "v1_component": v1c, "cath": cath, "reference": ref,
            "doi": ";".join(dict.fromkeys(dois)),
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
            "v1_component", "cath", "reference", "doi", "v1_overlay_via"]
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
        "n_v1_overlay_by_md5": sum(1 for n in ordered if n["v1_overlay_via"] == "md5"),
        "n_v1_overlay_by_accession": sum(1 for n in ordered if n["v1_overlay_via"] == "accession"),
    }
    Path(f"{out}_stats.json").write_text(json.dumps(stats, indent=2))
    print(f"[step1] {len(nodes)} nodes / {stats['n_unique_md5']} unique md5 "
          f"({n_blank} blank, {stats['n_md5_collapsed_rows']} md5-collapsed) -> {out}.fasta")


if __name__ == "__main__":
    main()
