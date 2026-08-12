#!/usr/bin/env python3
"""Emit the cluster-representative FASTA from an annotated clusters table.

Input is a Step-03 `alignment.tsv` (one row per cluster, carrying the centroid's
accession/pazy_id/sequence plus the connected-component id). Output is one record
per cluster with a five-field pipe header:

    >identifier|accession|genbank_accessions|pazy_id|component

  * `identifier`   -- the centroid's Step-01 union id, prefixed (default `PL`).
                      Every centroid *is* a union row, so it keeps that row's
                      `plasticome_id` and the FASTA stays joinable back to
                      `01-union.tsv`. The ids are therefore sparse: 411 centroids
                      drawn from 609 union rows means the numbering has gaps, and
                      the largest id is the largest union id, not 411. It is NOT a
                      1..N re-index over the emitted rows -- that is what this
                      script used to do, and it silently renamed every centroid
                      past the first renumbering point (`WP_015787089.1` was union
                      18 and was published as `PL17`).
  * `accession`    -- the representative centroid's accession, i.e. the one the
                      sequence was pulled from. Blank for centroids with no
                      accession (placeholders like `jmPE13`).
  * `genbank_accessions` -- other accessions mapping to the same sequence.
                      Always emitted empty here; the pipe is kept so the field
                      count is fixed and the column can be filled in later.
  * `pazy_id`      -- blank on centroids that have none.
  * `component`    -- the alignment-graph component id (`C001`...).

Resolving the union id, first hit wins:

    1. `rep_plasticome_id`, carried through steps 2 and 3 from the union row.
    2. the `U####` prefix of `rep_label`, which is the union row's 1-based
       position -- this is what a run built before step 2 started carrying
       `rep_plasticome_id` has to fall back on.
    3. the centroid's accession, when that accession names exactly one union row.

Routes 2 and 3 need the union table (`--union`, else `01-union.tsv` beside the
input). Whichever route resolves a row, the union row's sequence is checked
against the centroid's; a mismatch, an unresolvable row, or a repeated id is a
hard error rather than a quietly wrong identifier.

Ordering (mirrors `lib/01 union/build_union.py:sort_key`, with the amino-acid
sequence replacing input order as the final tiebreak):

    1. pazy_id, numerically; rows *without* a pazy_id sort after every row that
       has one -- same rule build_union.py applies to the union table.
    2. component, numerically (`C012` -> 12).
    3. aa_sequence, lexically.

Because pazy_id is unique across the rows that have one, keys 2 and 3 only
discriminate among the blank-pazy_id tail. Record order no longer determines the
numbering, so `--sort id` is available to emit in identifier order instead.

Example
-------
  python lib/fasta/clusters_to_fasta.py \\
      runs/2026-07-30.plasticome.v1.212.union-spec/03-alignment.tsv \\
      -o runs/2026-07-30.plasticome.v1.212.union-spec/04-plasticome.v1.212.union-spec.fasta
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tsv_to_fasta import duplicate_columns, normalize, wrap  # noqa: E402

INF = float("inf")

# `U0018|WP_015787089.1` -- the number is the union row's 1-based position, which is
# also its plasticome_id, and the accession after the pipe only keeps labels unique.
LABEL_RE = re.compile(r"^U0*(\d+)\|")

TSV_COLS = ["identifier", "accession", "genbank_accessions", "pazy_id",
            "component", "seq_len", "aa_sequence"]


def num(x: str) -> float:
    """Leading integer in a key, or +inf so blank/non-numeric sorts last.

    Handles both bare pazy ids (`16`) and prefixed component ids (`C012`); the
    union table's `num()` only ever saw the former.
    """
    m = re.search(r"\d+", x or "")
    return int(m.group()) if m else INF


def sort_key(row: dict, cols: dict) -> tuple:
    pazy = (row.get(cols["pazy"]) or "").strip()
    comp = (row.get(cols["component"]) or "").strip()
    seq = normalize(row.get(cols["seq"], ""))
    return (
        0 if pazy else 1,        # pazy_id rows first, blanks last
        num(pazy), pazy,
        num(comp), comp,
        seq,
    )


class Union:
    """The Step-01 union table, indexed the three ways the id resolver needs it."""

    def __init__(self, path: Path, rows: list[dict]):
        self.path, self.rows = path, rows
        self.by_id: dict[str, dict] = {}
        self.by_acc: dict[str, list[dict]] = defaultdict(list)
        for r in rows:
            pid = (r.get("plasticome_id") or "").strip()
            if pid:
                self.by_id.setdefault(pid, r)
            acc = (r.get("accession") or "").strip()
            if acc:
                self.by_acc[acc].append(r)

    def at(self, pos: int) -> dict | None:
        """The row at a 1-based `U####` position."""
        return self.rows[pos - 1] if 1 <= pos <= len(self.rows) else None

    @staticmethod
    def pid_of(row: dict, pos: int | None = None) -> str:
        """A union row's id, falling back to its position for a table without one."""
        return (row.get("plasticome_id") or "").strip() or (str(pos) if pos else "")


def find_union(tsv: Path) -> Path | None:
    """`01-union.tsv` beside the input, else the sole `01-*.tsv` in that directory."""
    canon = tsv.parent / "01-union.tsv"
    if canon.exists():
        return canon
    cands = sorted(tsv.parent.glob("01-*.tsv"))
    return cands[0] if len(cands) == 1 else None


def resolve_id(row: dict, cols: dict, union: Union | None) -> tuple[str, dict | None, str]:
    """(union id, union row, route) for one centroid; ("", None, "") if unresolvable.

    The union row comes back alongside the id so the caller can verify the two
    really describe the same sequence.
    """
    pid = (row.get(cols["id"]) or "").strip()
    if pid:
        return pid, (union.by_id.get(pid) if union else None), cols["id"]
    if union is None:
        return "", None, ""
    m = LABEL_RE.match((row.get(cols["label"]) or "").strip())
    if m:
        pos = int(m.group(1))
        u = union.at(pos)
        if u is not None:
            return Union.pid_of(u, pos), u, cols["label"]
    acc = (row.get(cols["accession"]) or "").strip()
    hits = union.by_acc.get(acc, ()) if acc else ()
    if len(hits) == 1:                   # an accession on two union rows names neither
        return Union.pid_of(hits[0]), hits[0], cols["accession"]
    return "", None, ""


def build_records(rows, cols, prefix, pad, union, order):
    """(records, errors); a record is (identifier, row, sequence), errors are strings.

    Blank-sequence rows are dropped, as they always were. Everything else must come
    out with exactly one union id: an identifier this script invents is an identifier
    that does not point back at `01-union.tsv`, which is the whole point of the file.
    """
    keep = [r for r in rows if normalize(r.get(cols["seq"], ""))]
    recs, errors, unresolved, mismatched = [], [], [], []
    for i, r in enumerate(keep, 1):
        pid, urow, route = resolve_id(r, cols, union)
        seq = normalize(r[cols["seq"]])
        # Name the row by whatever identifies it to a human; a centroid with neither
        # a label nor an accession is exactly the kind that fails to resolve, so the
        # row number is the last resort rather than an unhelpful placeholder.
        label = ((r.get(cols["label"]) or "").strip()
                 or (r.get(cols["accession"]) or "").strip()
                 or f"row {i}")
        if not pid:
            unresolved.append(label)
            continue
        if urow is not None and normalize(urow.get("aa_sequence", "")) != seq:
            mismatched.append(f"{label} -> union id {pid} (via {route})")
        ident = f"{prefix}{int(pid):0{pad}d}" if pad and pid.isdigit() else f"{prefix}{pid}"
        recs.append((ident, r, seq))

    if unresolved:
        hint = ("" if union is not None else
                "; no union table was found beside the input -- pass --union")
        errors.append(f"{len(unresolved)} centroid(s) have no union id{hint}: "
                      + ", ".join(unresolved[:10])
                      + (" ..." if len(unresolved) > 10 else ""))
    if mismatched:
        errors.append(f"{len(mismatched)} centroid(s) resolved to a union row with a "
                      f"different sequence: " + ", ".join(mismatched[:10])
                      + (" ..." if len(mismatched) > 10 else ""))
    seen: dict[str, int] = defaultdict(int)
    for ident, _, _ in recs:
        seen[ident] += 1
    dupes = sorted(i for i, n in seen.items() if n > 1)
    if dupes:
        errors.append(f"{len(dupes)} identifier(s) assigned to more than one centroid: "
                      + ", ".join(dupes[:10]) + (" ..." if len(dupes) > 10 else ""))

    if order == "id":
        recs.sort(key=lambda t: (num(t[0]), t[0]))
    else:
        recs.sort(key=lambda t: sort_key(t[1], cols))
    return recs, errors


def header_of(ident: str, row: dict, cols: dict) -> str:
    get = lambda c: (row.get(cols[c]) or "").strip()  # noqa: E731
    # genbank_accessions is intentionally empty -- the pipe is a placeholder.
    return f"{ident}|{get('accession')}||{get('pazy')}|{get('component')}"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", type=Path, help="Step-03 alignment/clusters TSV")
    ap.add_argument("-o", "--out", default=None,
                    help="output FASTA path, or '-' for stdout "
                         "(default: input name with a .fasta suffix)")
    ap.add_argument("--prefix", default="PL",
                    help="identifier prefix (default: PL)")
    ap.add_argument("--pad", type=int, default=0,
                    help="zero-pad the identifier to this width (default: 0, unpadded)")
    ap.add_argument("--width", type=int, default=0,
                    help="sequence line-wrap width; 0 = one line per record "
                         "(default: 0, matching the prior centroid FASTA)")
    ap.add_argument("--union", type=Path, default=None,
                    help="Step-01 union TSV, the source of the identifiers "
                         "(default: 01-union.tsv beside the input)")
    ap.add_argument("--tsv-out", type=Path, default=None,
                    help="also write the same records as a TSV "
                         f"({', '.join(TSV_COLS)})")
    ap.add_argument("--sort", choices=("spec", "id"), default="spec",
                    help="record order: 'spec' is the union sort spec "
                         "(pazy_id, component, sequence; default), 'id' is "
                         "ascending identifier")
    ap.add_argument("--accession-col", default="rep_accession")
    ap.add_argument("--pazy-col", default="rep_pazy_id")
    ap.add_argument("--component-col", default="component_id")
    ap.add_argument("--seq-col", default="rep_aa_sequence")
    ap.add_argument("--id-col", default="rep_plasticome_id",
                    help="column carrying the centroid's union id (default: "
                         "rep_plasticome_id)")
    ap.add_argument("--label-col", default="rep_label",
                    help="column carrying the centroid's U#### union label, used "
                         "when --id-col is blank (default: rep_label)")
    args = ap.parse_args()

    if not args.tsv.exists():
        print(f"error: no such file: {args.tsv}", file=sys.stderr)
        return 1

    cols = {"accession": args.accession_col, "pazy": args.pazy_col,
            "component": args.component_col, "seq": args.seq_col,
            "id": args.id_col, "label": args.label_col}
    # id/label are the fallback chain, not data this script emits: a clusters table
    # predating either column still resolves through whatever it does carry.
    required = ["accession", "pazy", "component", "seq"]

    with args.tsv.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        have = reader.fieldnames or []
        dupes = duplicate_columns(have)
        if dupes:
            print(f"error: duplicate column name(s) {dupes} in {args.tsv.name}; "
                  f"only the last of each is readable, which silently drops data.",
                  file=sys.stderr)
            return 1
        missing = [cols[c] for c in required if cols[c] not in have]
        if missing:
            print(f"error: column(s) {missing} not in {args.tsv.name}; "
                  f"available: {have}", file=sys.stderr)
            return 1
        rows = list(reader)

    union_path = args.union or find_union(args.tsv)
    union = None
    if union_path is not None:
        if not union_path.exists():
            print(f"error: no such union table: {union_path}", file=sys.stderr)
            return 1
        with union_path.open(newline="") as fh:
            union = Union(union_path, list(csv.DictReader(fh, delimiter="\t")))

    recs, errors = build_records(rows, cols, args.prefix, args.pad, union, args.sort)
    if errors:
        for e in errors:
            print(f"error: {e}", file=sys.stderr)
        print("       identifiers must come from the union table; refusing to "
              "invent them.", file=sys.stderr)
        return 1

    if args.out == "-":
        out_fh, out_name, close = sys.stdout, "<stdout>", False
    else:
        out_path = Path(args.out) if args.out else args.tsv.with_suffix(".fasta")
        out_fh, out_name, close = out_path.open("w"), str(out_path), True

    n = 0
    try:
        for ident, row, seq in recs:
            out_fh.write(f">{header_of(ident, row, cols)}\n")
            if args.width > 0:
                for line in wrap(seq, args.width):
                    out_fh.write(line + "\n")
            else:
                out_fh.write(seq + "\n")
            n += 1
    finally:
        if close:
            out_fh.close()

    print(f"{out_name}: wrote {n} records from {len(rows)} rows "
          f"({len(rows) - n} blank seq)", file=sys.stderr)

    if args.tsv_out:
        # Same records, same order, exploded out of the pipe header, so the two
        # deliverables cannot drift apart in identifier or content.
        with args.tsv_out.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(TSV_COLS)
            for ident, row, seq in recs:
                get = lambda c: (row.get(cols[c]) or "").strip()  # noqa: E731
                w.writerow([ident, get("accession"), "", get("pazy"),
                            get("component"), len(seq), seq])
        print(f"{args.tsv_out}: wrote {n} rows", file=sys.stderr)

    if union is not None:
        print(f"identifiers from {union_path.name} ({len(union.rows)} rows)",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(0)
