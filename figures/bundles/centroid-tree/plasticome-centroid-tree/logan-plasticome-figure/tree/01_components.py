#!/usr/bin/env python3
"""
Step 1: turn the SSN components into the units of the tree.

A component is a connected component of the 30 percent identity, 50 percent
coverage graph from ../ssn (the family cluster definition of panel a-i). Each
one gets its own alignment and tree. Components are then binned into fold
groups by the CATH superfamily of the catalytic domain. Only components in the
same fold group can share an ancestor, so only those are grafted together.

Two cleanups before inference:
  identical sequences  the union set carries the same protein twice when both
                       PAZy pulls list it. They are collapsed to one tip, named
                       for the smallest plasticome_id, and the tip keeps every
                       record it stands for.
  fragments            a sequence shorter than half its component's median
                       length is left out of inference and listed as such.

Writes work/tips.tsv, work/components.tsv and work/components/cNN/seqs.faa.
"""

import collections
import os
import statistics

from common import (HITS, NODES, SET, WORK, cdir, family_of, in_set, load_rows,
                    read_tsv, shape_of, substrate_of, write_fasta, write_tsv)

FRAGMENT_FRACTION = 0.5


def fold_key(rows_in_comp):
    """CATH superfamily of the catalytic domain, or the domain name when the
    record has no CATH assignment."""
    cnt = collections.Counter()
    for r in rows_in_comp:
        cath = (r["catalytic-cath"] or "NA").strip()
        if cath in ("NA", "-", ""):
            cath = "domain:" + (r["domain"] or "NA").split(";")[0].strip()
        cnt[cath] += 1
    return sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def main():
    rows = load_rows()
    nodes = {r["plasticome_id"]: r for r in read_tsv(NODES)}
    comp_members = collections.defaultdict(list)
    for pid, n in nodes.items():
        if in_set(rows[pid]):
            comp_members[int(n["component"])].append(pid)
    os.makedirs(WORK, exist_ok=True)

    # mean identity to the rest of the component, for the medoid ranking
    ident = collections.defaultdict(dict)
    with open(HITS) as fh:
        for line in fh:
            q, s, pid, *_ = line.split("\t")
            if q != s:
                ident[q][s] = max(float(pid), ident[q].get(s, 0.0))

    tips, comps = [], []
    for c in sorted(comp_members):
        members = sorted(comp_members[c], key=int)
        med = statistics.median(len(rows[m]["aa_sequence"].strip()) for m in members)

        by_seq = collections.OrderedDict()
        for m in members:
            by_seq.setdefault(rows[m]["aa_sequence"].strip().upper(), []).append(m)

        kept = {}
        for seq, recs in by_seq.items():
            tip = recs[0]
            status = "fragment" if len(seq) < FRAGMENT_FRACTION * med else "kept"
            subs = sorted({substrate_of(rows[r]["substrate"]) for r in recs})
            # records of one protein that disagree on substrate are a
            # combination by definition; NR beside a named plastic is not
            named = [s for s in subs if s != "NR"]
            sub = "NR" if not named else named[0] if len(named) == 1 else "Combi."
            if status == "kept":
                kept[tip] = seq
            r0 = rows[tip]
            tips.append([tip, ";".join(recs), r0["enzyme_name"], r0["accession"],
                         c, len(seq), status, sub, shape_of(r0),
                         family_of(r0["domain"]).replace("\n", " "),
                         r0["uniprot_accession"], r0["pdb_accession"]])

        os.makedirs(cdir(c), exist_ok=True)
        write_fasta(os.path.join(cdir(c), "seqs.faa"), kept)

        # Medoid ranking over the kept tips: highest mean identity to the other
        # tips first. Step 3 walks down this list to the first tip that has a
        # structure, so the representative is as central as the data allows.
        def centrality(t):
            others = [o for o in kept if o != t]
            if not others:
                return 0.0
            return sum(ident[t].get(o, 0.0) for o in others) / len(others)
        ranking = sorted(kept, key=lambda t: (-centrality(t), int(t)))

        fam = collections.Counter(family_of(rows[m]["domain"]).replace("\n", " ")
                                  for m in members).most_common(1)[0][0]
        comps.append([c, len(members), len(kept), fold_key([rows[m] for m in members]),
                      fam, int(med), ";".join(ranking)])

    # fold groups, largest first, numbered so the order is stable
    by_fold = collections.defaultdict(list)
    for row in comps:
        by_fold[row[3]].append(row)
    order = sorted(by_fold, key=lambda k: (-sum(r[2] for r in by_fold[k]), k))
    gid = {k: i for i, k in enumerate(order)}
    comps = [[r[0], gid[r[3]]] + r[1:] for r in comps]

    write_tsv(os.path.join(WORK, "tips.tsv"),
              ["tip", "records", "enzyme_name", "accession", "component", "length",
               "status", "substrate_group", "shape", "family_group",
               "uniprot_accession", "pdb_accession"], tips)
    write_tsv(os.path.join(WORK, "components.tsv"),
              ["component", "fold_group", "n_records", "n_tips", "fold_key",
               "family", "median_length", "medoid_ranking"], comps)

    n_kept = sum(1 for t in tips if t[6] == "kept")
    n_rec = sum(1 for r in rows.values() if in_set(r))
    print(f"set: {SET}. {n_rec} of {len(rows)} records -> {len(tips)} distinct sequences, "
          f"{n_kept} kept, {len(tips) - n_kept} fragments left out")
    for t in tips:
        if t[6] == "fragment":
            print(f"  fragment: {t[0]} {t[2]} ({t[5]} aa, component {t[4]})")
    print(f"{len(comps)} components in {len(order)} fold groups")
    for k in order:
        cs = by_fold[k]
        print(f"  g{gid[k]:02d} {k:<48} {len(cs):>2} components, "
              f"{sum(r[2] for r in cs):>3} tips")


if __name__ == "__main__":
    main()
