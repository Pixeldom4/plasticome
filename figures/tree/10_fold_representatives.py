#!/usr/bin/env python3
"""
Step 10: one representative enzyme per fold, and the FASTA of those sequences.

A fold here is a slice of the pie from step 8: the components whose catalytic
CATH superfamily is the same and whose homology step 4 confirmed. Singleton
folds, the ones that are a single distinct sequence and are drawn in the
shared grey, are skipped. The sequence is its own representative, so a table
row for it would say nothing the tree does not already say.

How the representative is picked. Candidates are the distinct sequences of the
fold that the trees kept, so fragments are out. They are ranked by:

  1. curated seed      an enzyme from source-data/plasticome-curated-seeds
  2. solved structure  has_structure = yes in the union table
  3. largest component the fold's biggest SSN component first
  4. medoid centrality the medoid ranking from step 1, highest mean identity
                       to the rest of its component first
  5. plasticome_id     so the choice is stable

Seed beats structure beats size beats centrality, in that order, and the
`chosen_for` column records which rungs the winner actually stood on.

Categorising the folds. The CATH class and architecture come from the
accession. The mechanism column is hand-assigned, one label per fold, and is
the only column here that is not read off the table or computed.

Writes results/fold_representatives.tsv and results/fold_representatives.faa.
"""

import collections
import os
import re

from common import (NODES, OUT, SET, SUBSTRATE_ORDER, WORK, final_groups,
                    in_set, load_rows, read_tsv, seed_nodes, write_tsv)

CLASS = {"1": "Mainly alpha", "2": "Mainly beta", "3": "Alpha beta",
         "4": "Few secondary structures", "6": "Special"}
ARCH = {"1.10": "Orthogonal bundle", "2.40": "Beta barrel", "2.60": "Sandwich",
        "3.20": "Alpha-beta barrel", "3.40": "3-layer (aba) sandwich",
        "3.60": "4-layer sandwich", "3.90": "Alpha-beta complex"}
ORIGIN = [("D", "Bacterial"), ("h", "Eukaryotic"), ("o", "Metagenomic")]
MISSING = ("", "NA", "-")

# Hand-assigned, one per catalytic CATH superfamily of a non-singleton fold.
# Keyed on the catalytic accession, "none" for the fold with no accession.
MECHANISM = {
    "3.40.50.1820": "Serine hydrolase, Ser-His-Asp triad on the alpha/beta hydrolase fold",
    "3.60.70.12":   "Ntn amidohydrolase, autoprocessed N-terminal Thr nucleophile",
    "3.90.1300.10": "Amidase signature family, Ser-cisSer-Lys triad",
    "none":         "Rubber oxygenase, heme b MpaB/Lcp family",
    "3.40.50.200":  "Subtilisin-like serine protease, Asp-His-Ser triad",
    "2.40.10.10":   "Chymotrypsin-like serine protease, His-Asp-Ser triad",
    "3.40.710.10":  "Penicillin-recognising serine hydrolase, Ser-Lys dyad",
    "1.10.760.10":  "Rubber oxygenase, c-type diheme",
    "2.60.40.420":  "Multicopper oxidase, laccase-like",
    "3.40.50.1110": "SGNH esterase, Ser-Gly-Asn-His catalytic tetrad",
}


def ids(cell):
    return [x.strip() for x in (cell or "").split(";") if x.strip() not in MISSING]


def names(cell):
    return [x.strip() for x in re.split(r";(?! )", cell or "") if x.strip() not in MISSING]


def clean(name):
    name = " ".join(name.replace("[topology]", "").split())
    return name[0].upper() + name[1:]


def wrap(fh, header, seq, width=60):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def main():
    rows = {pid: r for pid, r in load_rows().items() if in_set(r)}
    comp_of = {r["plasticome_id"]: int(r["component"]) for r in read_tsv(NODES)}
    group_of, verdict = final_groups()
    comps = {int(r["component"]): r for r in read_tsv(os.path.join(WORK, "components.tsv"))}
    tips = read_tsv(os.path.join(OUT, "tree_tips.tsv"))
    tip_of = {rec: t for t in tips for rec in t["records"].split(";")}
    seeds = {pid: name for name, pid in seed_nodes(rows).items()}

    # --- group the records, and the domain names behind the CATH ids --------
    by_g, name_of = collections.defaultdict(list), {}
    for pid, r in rows.items():
        cat, acc, nm = ids(r["catalytic-cath"]), ids(r["accessory-cath"]), names(r["domain"])
        if not cat and nm:
            cat = ["none"]
        for c, n in zip(cat + acc, nm):
            name_of.setdefault(c if c != "none" else f"none:{n}", clean(n))
        by_g[f"g{group_of[comp_of[pid]]:02d}"].append(
            {"pid": pid, "comp": comp_of[pid], "cat": cat, "acc": acc,
             "tip": tip_of.get(pid)})

    order = sorted(by_g, key=lambda g: (-len({r["tip"]["tip"] for r in by_g[g] if r["tip"]}), g))

    out_rows, chosen = [], []
    for g in order:
        rs = by_g[g]
        seqs = {r["tip"]["tip"] for r in rs if r["tip"]}
        if len(seqs) < 2:
            continue                       # singleton fold, skip

        # --- candidates: one per distinct sequence the trees kept -----------
        recs_of = collections.defaultdict(list)
        for r in rs:
            if r["tip"]:
                recs_of[r["tip"]["tip"]].append(r)
        cand = []
        for tip, group in recs_of.items():
            comp = group[0]["comp"]
            ranking = comps[comp]["medoid_ranking"].split(";")
            if tip not in ranking:
                continue                   # fragment, left out of the trees
            seed = next((seeds[r["pid"]] for r in group if r["pid"] in seeds), "")
            struct = any(rows[r["pid"]]["has_structure"].strip().lower() == "yes"
                         for r in group)
            cand.append({"tip": tip, "group": group, "comp": comp, "seed": seed,
                         "struct": struct, "rank": ranking.index(tip) + 1,
                         "of": len(ranking), "size": int(comps[comp]["n_tips"])})
        best = min(cand, key=lambda c: (not c["seed"], not c["struct"], -c["size"],
                                        c["rank"], int(c["tip"])))

        # the record of the winner that carries the most metadata
        rep = max((rows[r["pid"]] for r in best["group"]),
                  key=lambda r: sum(1 for v in r.values() if v.strip() not in MISSING))
        why = ", ".join(filter(None, [
            "curated seed" if best["seed"] else "",
            "solved structure" if best["struct"] else "",
            f"largest component ({best['size']} seqs)" if best["size"] > 1 else "",
            f"medoid rank {best['rank']} of {best['of']}"]))

        cats = [c for c, _ in collections.Counter(c for r in rs for c in r["cat"]).most_common()]
        if not cats:
            cats = ["none"]
        acc = collections.Counter(a for r in rs for a in r["acc"])
        key = cats[0]
        # a fold with no CATH accession is named by its domain name instead
        fold = (name_of[key] if key != "none"
                else clean(names(rows[rs[0]["pid"]]["domain"])[0]))
        tipinfo = [r["tip"] for r in rs if r["tip"]]
        sub = collections.Counter(t["substrate_group"] for t in {t["tip"]: t for t in tipinfo}.values())
        org = collections.Counter(t["shape"] for t in {t["tip"]: t for t in tipinfo}.values())
        v = {verdict[r["comp"]] for r in rs}

        out_rows.append([
            g, fold.split(";")[0].strip(), "; ".join(cats),
            CLASS.get(key.split(".")[0], "") if key[:1].isdigit() else "",
            ARCH.get(".".join(key.split(".")[:2]), "") if key[:1].isdigit() else "",
            MECHANISM.get(key, ""),
            len({r["comp"] for r in rs}), len(seqs), len(rs),
            *[sub.get(s, 0) for s in SUBSTRATE_ORDER],
            *[org.get(k, 0) for k, _ in ORIGIN],
            "confirmed by structure" if "confirmed" in v else "single component",
            "; ".join(f"{a} {name_of[a]} ({n})" for a, n in acc.most_common()),
            best["seed"] or rep["enzyme_name"], best["tip"], rep["enzyme_name"],
            rep["accession"], rep["uniprot_accession"], rep["pdb_accession"],
            rep["has_structure"], len(rep["aa_sequence"].strip()),
            rep["substrate"], rep["organism"], rep["superkingdom"],
            f"c{best['comp']:02d}", why])
        chosen.append((g, best, rep))

    head = (["slice", "fold", "catalytic_cath", "cath_class", "cath_architecture",
             "mechanism", "components", "sequences", "records"]
            + [s.lower().replace(".", "") for s in SUBSTRATE_ORDER]
            + [lab.lower() for _, lab in ORIGIN]
            + ["homology", "accessory_domains", "representative", "rep_plasticome_id",
               "rep_enzyme_name", "rep_accession", "rep_uniprot", "rep_pdb",
               "rep_has_structure", "rep_length_aa", "rep_substrate", "rep_organism",
               "rep_superkingdom", "rep_component", "chosen_for"])
    write_tsv(os.path.join(OUT, "fold_representatives.tsv"), head, out_rows)
    REP, SEQ = head.index("representative"), head.index("sequences")

    path = os.path.join(OUT, "fold_representatives.faa")
    with open(path, "w") as fh:
        for row, (g, best, rep) in zip(out_rows, chosen):
            wrap(fh, f"{g}|{row[REP]}|{rep['accession']}|{row[2]}|{row[1]}",
                 rep["aa_sequence"].strip())

    print(f"set: {SET}. {len(out_rows)} non-singleton folds of {len(order)}")
    for row in out_rows:
        print(f"  {row[0]}  {row[2]:<24} {row[REP]:<22} {row[1]:<42} "
              f"{row[SEQ]:>3} seqs  [{row[-1]}]")


if __name__ == "__main__":
    main()
