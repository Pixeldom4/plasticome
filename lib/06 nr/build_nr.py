#!/usr/bin/env python3
"""Step 6 - the 100% non-redundant view of the union, one row per distinct sequence.

Input is `01-union.tsv` and nothing else. Output is set C of the three deliverable
sets:

    A  01-union.tsv     609 rows  the starting set of record          (step 1)
    B  04-<run>.fasta   411 rows  90% centroids, the node set         (steps 2-4)
    C  06-nr.tsv        493 rows  100% non-redundant, this script     (step 6)

C exists so a downstream `nr` search can ask, once per distinct amino-acid
sequence, which other GenBank/UniProt/PDB accessions retrieve it. Collapsing the
union first turns 609 queries into 493 and gives every alternate accession
already present *inside* the union for free.

B and C are parallel branches off step 1, not a chain
-----------------------------------------------------
C never feeds steps 2 or 3. Building B on top of C would shrink the database
steps 2 and 3 search, and a smaller database means smaller e-values for the same
alignment score, so the fixed `evalue < 1e-5` post-filter would become
effectively more permissive -- more edges, fewer and larger components. That is
the same database-size confound that made `runs/2026-08-05.align-before-cluster`
uninterpretable. B is the delivered artifact and its lineage stays frozen; C is
built beside it.

Consequently this script computes no alignment, no clustering and no component
partition, and it mints no identifier. `--only 6` runs in a directory holding
nothing but `01-union.tsv`.

Hashing, not clustering at 100%
-------------------------------
`usearch cluster_fast -id 1.0` and `cd-hit -c 1.0` are not exact-identity
operators: both measure identity over an alignment with a coverage model, so a
300 aa protein wholly contained in a 320 aa protein collapses at 1.0 under
defaults. That is a domain match, not "the same amino acid sequence", and it is
the wrong relation for an accession lookup. A hash groupby is O(n),
deterministic, needs no binary, and drops the Docker dependency for this step.

`normalize()` is imported from `lib/03 alignment/scripts/step1_nodes.py` rather
than redefined. Two identity spaces would diverge on exactly the rows that
differ by a trailing `*`, lowercase masking or stray whitespace, and the join
between C and everything else would then fail silently. The emitted
`aa_sequence` is the *normalized* one, so what was hashed and what gets queried
are the same bytes.

Rules
-----
  * Representative is the smallest `plasticome_id` in the md5 group, matching
    step 3's smallest-PL-represents convention.
  * `member_plasticome_ids` and `member_accessions` are positionally parallel
    `;`-delimited lists in ascending `plasticome_id`, so the representative is
    always first and `zip()` recovers which accession belongs to which row. A
    union row with no accession leaves an empty slot rather than shortening the
    accession list -- 41 of the 609 rows on the reference union have none.
  * **Accession versions are kept.** Step 2 folds seeds version-insensitively,
    which is right for seeds and wrong here: `WP_012345.1` and `.2` can be
    different sequences, and 06 is precisely about accessions that retrieve an
    identical one.
  * `;` inside an accession is a hard error, not something to escape.

Fixture
-------
609 union rows -> 493 distinct sequences, 109 duplicate groups covering 225 rows
(225 - 109 = 116 eliminated, 609 - 116 = 493). Any future run against this union
that does not produce 493 rows means the union or `normalize()` changed.

Deliverables (beside the input, unless overridden)
--------------------------------------------------
  06-nr.tsv                             one row per distinct sequence
  06-nr.fasta                           5-field pipe header, normalized sequence
  06-nr.intermediates/provenance.json   input md5, normalize() hash, counts

Example
-------
  python "lib/06 nr/build_nr.py" runs/2026-08-06.final-usearch.2/01-union.tsv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

# normalize() is the one shared definition of sequence identity in this repo; see
# the module docstring for why 06 imports it instead of carrying its own copy.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "03 alignment" / "scripts"))
try:
    from step1_nodes import normalize  # noqa: E402
except ImportError as exc:  # pragma: no cover - a broken checkout, not a run-time path
    raise SystemExit(
        f"error: cannot import normalize() from 'lib/03 alignment/scripts/step1_nodes.py': {exc}\n"
        "       step 6 shares step 3's sequence identity rather than defining a second one."
    )

csv.field_size_limit(10 ** 7)  # aa_sequence makes for long fields

DELIM = ";"

# Order matters to a human reading the file, not to the code: rows are written
# with DictWriter and read back with DictReader, so this list is the single place
# the layout is decided. `seq_md5` sits last, after the long `aa_sequence`, so the
# identifying columns are the ones visible without scrolling; it is still the key
# every other file joins on.
COLS = ["rep_plasticome_id", "rep_accession", "n_members",
        "member_plasticome_ids", "member_accessions", "seq_len", "aa_sequence",
        "seq_md5"]

REQUIRED_IN = ["plasticome_id", "accession", "aa_sequence"]


def pl_num(pid: str) -> int:
    """Leading integer of a plasticome_id, or +inf so a non-numeric id sorts last.

    Mirrors `step1_nodes.pl_num`. The union's ids are bare integers today; the
    fallback only keeps the sort total if that ever stops being true.
    """
    m = re.search(r"\d+", pid or "")
    return int(m.group()) if m else 10 ** 9


def sort_key(pid: str) -> tuple:
    return (pl_num(pid), pid)


def read_union(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        have = reader.fieldnames or []
        missing = [c for c in REQUIRED_IN if c not in have]
        if missing:
            raise SystemExit(f"error: column(s) {missing} not in {path.name}; available: {have}")
        return list(reader)


def check(rows: list[dict], path: Path) -> None:
    """Preconditions. Each of these means step 1 did not finish, so none are warnings."""
    blank = [(r.get("plasticome_id") or "?") for r in rows
             if not normalize(r.get("aa_sequence", ""))]
    if blank:
        raise SystemExit(
            f"error: {len(blank)} row(s) in {path.name} have a blank aa_sequence: "
            + ", ".join(blank[:10]) + (" ..." if len(blank) > 10 else "")
            + "\n       step 1 backfills every sequence from its accession; a blank means "
              "fetch_sequences.py did not complete."
        )

    seen: dict[str, int] = defaultdict(int)
    for r in rows:
        seen[(r.get("plasticome_id") or "").strip()] += 1
    dupes = sorted((p for p, n in seen.items() if n > 1), key=sort_key)
    if dupes:
        raise SystemExit(
            f"error: {len(dupes)} plasticome_id(s) appear on more than one row of "
            f"{path.name}: " + ", ".join(dupes[:10]) + (" ..." if len(dupes) > 10 else "")
            + "\n       the id is the join key for every downstream set; it has to be unique."
        )

    # `;` separates accessions inside a field; `|` separates the fields. Either one
    # inside an accession corrupts the header, so both are fatal.
    bad = [f"{(r.get('plasticome_id') or '?')} ({(r.get('accession') or '').strip()})"
           for r in rows
           if DELIM in (r.get("accession") or "") or "|" in (r.get("accession") or "")]
    if bad:
        raise SystemExit(
            f"error: {len(bad)} accession(s) contain '{DELIM}' or '|', which separate "
            f"accessions and header fields respectively: "
            + ", ".join(bad[:10]) + (" ..." if len(bad) > 10 else "")
            + "\n       an accession carrying either is a hard error rather than "
              "something to escape."
        )


def group(rows: list[dict]) -> list[dict]:
    """md5 -> one output row, groups ordered by their representative's plasticome_id."""
    by_md5: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        seq = normalize(r.get("aa_sequence", ""))
        by_md5[hashlib.md5(seq.encode()).hexdigest()].append(r)

    out = []
    for md5, members in by_md5.items():
        members.sort(key=lambda r: sort_key((r.get("plasticome_id") or "").strip()))
        rep = members[0]
        seq = normalize(rep.get("aa_sequence", ""))
        ids = [(m.get("plasticome_id") or "").strip() for m in members]
        accs = [(m.get("accession") or "").strip() for m in members]
        out.append({
            # Header field 4. Deliberately not a TSV column: it is the union's
            # own value for the representative row, one join away on
            # rep_plasticome_id, so storing it here would be a second copy.
            "rep_pazy_id": (rep.get("pazy_id") or "").strip(),
            "seq_md5": md5,
            "rep_plasticome_id": ids[0],
            "rep_accession": accs[0],
            "n_members": len(members),
            "member_plasticome_ids": DELIM.join(ids),
            "member_accessions": DELIM.join(accs),
            "seq_len": len(seq),
            "aa_sequence": seq,
        })
    out.sort(key=lambda g: sort_key(g["rep_plasticome_id"]))
    return out


def write_tsv(path: Path, groups: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        # extrasaction: groups carry rep_pazy_id for the FASTA header, which is
        # not a column of this file.
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(groups)


def alternates_of(g: dict) -> str:
    """Header field 3: the accessions of a group's *other* members.

    The complement of field 2, so the representative's own accession never
    repeats in the list. Blanks are dropped and order is preserved, which is
    ascending `plasticome_id` because that is how the group was sorted. Duplicates
    are removed: two union rows can carry the same accession at the same version.
    """
    accs = [a.strip() for a in (g["member_accessions"] or "").split(DELIM)]
    rep = (g["rep_accession"] or "").strip()
    out, seen = [], {rep} if rep else set()
    for a in accs[1:]:
        if a and a not in seen:
            seen.add(a)
            out.append(a)
    return DELIM.join(out)


def header_of(g: dict, component: str = "") -> str:
    """`PL<n>|accession|alt_accessions|pazy_id|component`, five positional fields.

    Same shape and arity as the step-4 and step-5 headers, so one `split("|")`
    parser reads all three. Field 5 is a branch-B fact this script cannot know:
    `build_nr.py` leaves it empty and `crosswalk.py` re-emits the file with it
    filled once the component partition is on disk. See the module docstring.
    """
    return (f"PL{g['rep_plasticome_id']}|{g['rep_accession']}|{alternates_of(g)}"
            f"|{g.get('rep_pazy_id', '')}|{component}")


def write_fasta(path: Path, groups: list[dict], width: int = 0,
                components: dict | None = None) -> int:
    """One record per distinct sequence. Returns how many carry a component.

    `components` maps seq_md5 -> component_id; anything absent leaves field 5
    empty, with both pipes in place so the arity stays five.
    """
    components = components or {}
    n = 0
    with path.open("w") as fh:
        for g in groups:
            comp = components.get(g["seq_md5"], "")
            n += bool(comp)
            fh.write(f">{header_of(g, comp)}\n")
            seq = g["aa_sequence"]
            if width > 0:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")
    return n


def file_md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(src: Path, rows: list[dict], groups: list[dict]) -> dict:
    """What a rebuild has to reproduce to count as the same build.

    `normalize_md5` hashes the *source text* of the imported normalizer, so a
    silent edit to it shows up here rather than as an unexplained change in the
    group count three steps downstream.
    """
    dup = [g for g in groups if g["n_members"] > 1]
    return {
        "step": 6,
        "input": src.name,
        "input_md5": file_md5(src),
        "normalize_source": "lib/03 alignment/scripts/step1_nodes.py::normalize",
        "normalize_md5": hashlib.md5(inspect.getsource(normalize).encode()).hexdigest(),
        "n_union_rows": len(rows),
        "n_distinct_sequences": len(groups),
        "n_duplicate_groups": len(dup),
        "n_rows_in_duplicate_groups": sum(g["n_members"] for g in dup),
        "largest_group": max((g["n_members"] for g in groups), default=0),
        "keeps_accession_versions": True,
        "delimiter": DELIM,
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("union", type=Path, help="Step-01 union TSV")
    ap.add_argument("-o", "--out", type=Path, default=None,
                    help="output TSV path (default: 06-nr.tsv beside the input)")
    ap.add_argument("--fasta-out", type=Path, default=None,
                    help="output FASTA path (default: the TSV path with a .fasta suffix)")
    ap.add_argument("--intermediates", type=Path, default=None,
                    help="sidecar directory for provenance.json "
                         "(default: <out stem>.intermediates beside the output)")
    ap.add_argument("--width", type=int, default=0,
                    help="FASTA line-wrap width; 0 = one line per record (default: 0)")
    ap.add_argument("--expect", type=int, default=0,
                    help="fail unless exactly this many distinct sequences come out "
                         "(the reference union gives 493)")
    args = ap.parse_args()

    if not args.union.exists():
        print(f"error: no such file: {args.union}", file=sys.stderr)
        return 1

    rows = read_union(args.union)
    check(rows, args.union)
    groups = group(rows)

    # The universe is conserved: every union row lands in exactly one group.
    total = sum(g["n_members"] for g in groups)
    if total != len(rows):
        print(f"error: members sum to {total} but the union has {len(rows)} rows",
              file=sys.stderr)
        return 1
    if args.expect and len(groups) != args.expect:
        print(f"error: {len(groups)} distinct sequences, expected {args.expect}; "
              f"the union or normalize() changed", file=sys.stderr)
        return 1

    out = args.out or args.union.parent / "06-nr.tsv"
    fasta = args.fasta_out or out.with_suffix(".fasta")
    inter = args.intermediates or out.parent / f"{out.name.split('.')[0]}.intermediates"

    out.parent.mkdir(parents=True, exist_ok=True)
    inter.mkdir(parents=True, exist_ok=True)
    write_tsv(out, groups)
    n_comp = write_fasta(fasta, groups, args.width)
    prov = provenance(args.union, rows, groups)
    with (inter / "provenance.json").open("w") as fh:
        json.dump(prov, fh, indent=2)
        fh.write("\n")

    print(f"{out}: {len(groups)} distinct sequences from {len(rows)} union rows "
          f"({prov['n_duplicate_groups']} duplicate groups covering "
          f"{prov['n_rows_in_duplicate_groups']} rows, largest {prov['largest_group']})",
          file=sys.stderr)
    print(f"{fasta}: {len(groups)} records, {n_comp} with a component "
          f"(crosswalk.py fills field 5 once branch B is on disk)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(0)
