#!/usr/bin/env python3
"""
Step 9: the CATH domains of the plasticome as a table, coloured by fold.

Two views of the same annotation columns (catalytic-cath, accessory-cath,
domain) in the union table:

  Folds     one row per slice of the pie (step 8): the catalytic domain that
            defines the fold, with counts, substrates, origin, seeds, the
            accessory domains found beside it and the homology verdict.
  Domains   one row per CATH domain per fold it occurs in. A catalytic domain
            has one row. An accessory domain has a row for every fold whose
            enzymes carry it, so a module shared by two folds shows up twice.

Every row is filled with its fold's tint from 08_draw_pie.py, so the table and
the pie read as one figure.

Counting. "Sequences" are distinct sequences, as in the trees. "Records" are
rows of the union table, which lists a protein twice when both PAZy pulls
carry it. Substrate and origin columns count distinct sequences.

Domain names. The `domain` column joins names with ";" and some names contain
"; " themselves ("Peroxidase; domain 1"). A separator is a ";" with no space
after it, which lines the names up one to one with the CATH ids.

Writes out/cath_domains.xlsx (sheets Folds, Domains, Notes),
out/cath_folds.tsv and out/cath_domains.tsv.
"""

import collections
import importlib
import os
import re

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

from common import (NODES, OUT, SET, SUBSTRATE_ORDER, WORK, final_groups,
                    in_set, load_rows, read_tsv, seed_nodes, write_tsv)

pie = importlib.import_module("08_draw_pie")        # the fold tints live there

CLASS = {"1": "Mainly alpha", "2": "Mainly beta", "3": "Alpha beta",
         "4": "Few secondary structures", "6": "Special"}
ORIGIN = [("D", "Bacterial"), ("h", "Eukaryotic"), ("o", "Metagenomic")]
MISSING = ("", "NA", "-")


def ids(cell):
    return [x.strip() for x in (cell or "").split(";") if x.strip() not in MISSING]


def names(cell):
    return [x.strip() for x in re.split(r";(?! )", cell or "") if x.strip() not in MISSING]


def clean(name):
    name = " ".join(name.replace("[topology]", "").split())
    return name[0].upper() + name[1:]


def cath_class(cath):
    return CLASS.get(cath.split(".")[0], "") if cath[:1].isdigit() else ""


def main():
    rows = {pid: r for pid, r in load_rows().items() if in_set(r)}
    comp_of = {r["plasticome_id"]: int(r["component"]) for r in read_tsv(NODES)}
    group_of, verdict = final_groups()
    tips = read_tsv(os.path.join(OUT, "tree_tips.tsv"))
    tip_of = {rec: t for t in tips for rec in t["records"].split(";")}
    seeds = {pid: name for name, pid in seed_nodes(rows).items()}
    hom = {r["fold_group"]: r for r in read_tsv(os.path.join(OUT, "homology_check.tsv"))}
    bb = {r["fold_group"]: r for r in read_tsv(os.path.join(OUT, "backbone_summary.tsv"))}
    kept = collections.Counter()
    for r in read_tsv(os.path.join(OUT, "backbone_branches.tsv")):
        if int(r["ufboot"]) >= 95 and (r["in_aa_only_tree"] == "1" or r["in_tm_nj_tree"] == "1"):
            kept[r["fold_group"]] += 1

    # --- per record: fold group, catalytic and accessory domains -------------
    recs, name_of = [], {}
    for pid, r in rows.items():
        cat, acc, nm = ids(r["catalytic-cath"]), ids(r["accessory-cath"]), names(r["domain"])
        if not cat and nm:                       # a named domain with no CATH id
            cat = ["none"]
        if len(nm) != len(cat) + len(acc):
            raise SystemExit(f"record {pid}: {len(nm)} domain names for "
                             f"{len(cat) + len(acc)} CATH ids")
        for c, n in zip(cat + acc, nm):
            key = c if c != "none" else f"none:{n}"
            if name_of.setdefault(key, clean(n)) != clean(n):
                raise SystemExit(f"{c} has two names: {name_of[key]} / {clean(n)}")
        g = f"g{group_of[comp_of[pid]]:02d}"
        recs.append({"pid": pid, "g": g, "comp": comp_of[pid],
                     "cat": [c if c != "none" else f"none:{nm[0]}" for c in cat],
                     "acc": acc, "tip": tip_of.get(pid)})

    by_g = collections.defaultdict(list)
    for r in recs:
        by_g[r["g"]].append(r)
    order = sorted(by_g, key=lambda g: (-len({r["tip"]["tip"] for r in by_g[g] if r["tip"]}), g))

    def stats(rs):
        """Counts over a set of records; substrates and origin by sequence."""
        seqs = {r["tip"]["tip"]: r["tip"] for r in rs if r["tip"]}
        sub = collections.Counter(t["substrate_group"] for t in seqs.values())
        org = collections.Counter(t["shape"] for t in seqs.values())
        return ([len({r["comp"] for r in rs}), len(seqs), len(rs)]
                + [sub.get(s) for s in SUBSTRATE_ORDER]
                + [org.get(k) for k, _ in ORIGIN])

    def tint(g):
        key = next(t for t in tips if t["fold_group"] == g)["fold_key"]
        n = sum(1 for t in tips if t["fold_group"] == g)
        return pie.SINGLE_TINT if n == 1 else pie.FOLD_TINT.get(key, pie.SINGLE_TINT)

    def fold_label(g):
        cat = collections.Counter(c for r in by_g[g] for c in r["cat"]).most_common(1)[0][0]
        note = " (unconfirmed)" if any(verdict[r["comp"]] == "unconfirmed" for r in by_g[g]) else ""
        # the short form used on the pie: up to the first ";" of the name
        return name_of[cat].split(";")[0].strip() + note

    def show(cath):
        return "none" if cath.startswith("none:") else cath

    count_head = (["Components", "Sequences", "Records"] + SUBSTRATE_ORDER
                  + [lab for _, lab in ORIGIN])

    # --- sheet 1: folds ------------------------------------------------------
    fold_head = (["Slice", "Fold (catalytic domain)", "CATH", "CATH class", "Fill"]
                 + count_head + ["With PDB structure", "Curated seeds",
                                 "Accessory domains (records)", "Homology check",
                                 "Weakest TM link", "Backbone branches kept"])
    fold_rows = []
    for g in order:
        rs = by_g[g]
        cats = [c for c, _ in collections.Counter(c for r in rs for c in r["cat"]).most_common()]
        acc = collections.Counter(a for r in rs for a in r["acc"])
        pdb = len({r["tip"]["tip"] for r in rs if r["tip"]
                   and rows[r["pid"]]["has_structure"].strip().lower() == "yes"})
        v = {verdict[r["comp"]] for r in rs}
        check = ("not confirmed, stands alone" if "unconfirmed" in v else
                 "confirmed by structure" if "confirmed" in v else
                 "single component, nothing to check")
        h, b = hom.get(g, {}), bb.get(g)
        fold_rows.append(
            [g, fold_label(g), "; ".join(show(c) for c in cats),
             cath_class(cats[0]), tint(g)] + stats(rs)
            + [pdb or None,
               ", ".join(sorted((seeds[r["pid"]] for r in rs if r["pid"] in seeds),
                                key=str.lower)) or None,
               "; ".join(f"{a} {name_of[a]} ({n})" for a, n in acc.most_common()) or None,
               check,
               float(h["weakest_best_tm_in_group"])
               if h.get("weakest_best_tm_in_group", "NA") != "NA" and "unconfirmed" not in v
               and len({r["comp"] for r in rs}) > 1 else None,
               f"{kept[g]} of {b['branches']}" if b and int(b["branches"]) else None])

    # --- sheet 2: domains, one row per domain per fold -----------------------
    dom_head = (["CATH", "Domain", "CATH class", "Role", "Fold it occurs in", "Slice",
                 "Fill"] + count_head + ["Share of the fold's sequences"])
    dom_rows = []
    for g in order:
        rs = by_g[g]
        seen = collections.OrderedDict()
        for role in ("cat", "acc"):
            cnt = collections.Counter(d for r in rs for d in r[role])
            for d, _ in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])):
                seen[(d, role)] = [r for r in rs if d in r[role]]
        for (d, role), sub in seen.items():
            dom_rows.append([show(d), name_of[d], cath_class(d),
                             "catalytic" if role == "cat" else "accessory",
                             fold_label(g), g, tint(g)] + stats(sub) + [None])
        unassigned = [r for r in rs if not r["cat"]]
        if unassigned:
            dom_rows.append(["none", "No domain annotated", "", "unassigned",
                             fold_label(g), g, tint(g)] + stats(unassigned) + [None])

    write_tsv(os.path.join(OUT, "cath_folds.tsv"),
              [h.lower().replace(" ", "_") for h in fold_head],
              [["" if v is None else v for v in r] for r in fold_rows])
    write_tsv(os.path.join(OUT, "cath_domains.tsv"),
              [h.lower().replace(" ", "_") for h in dom_head[:-1]],
              [["" if v is None else v for v in r[:-1]] for r in dom_rows])

    # --- workbook ------------------------------------------------------------
    wb = Workbook()
    font = Font(name="Arial", size=10)
    bold = Font(name="Arial", size=10, bold=True)
    thin = Side(style="thin", color="8A8A8A")
    hair = Side(style="hair", color="C8C8C4")

    def sheet(ws, head, data, widths, wrap_cols=()):
        ws.append(head)
        for r in data:
            ws.append(r)
        fill_col = head.index("Fill") + 1
        for j in range(1, len(head) + 1):
            c = ws.cell(row=1, column=j)
            c.font, c.border = bold, Border(bottom=thin)
            c.alignment = Alignment(wrap_text=True, vertical="bottom",
                                    horizontal="left" if j <= fill_col or j in wrap_cols
                                    else "right")
            ws.column_dimensions[get_column_letter(j)].width = widths.get(head[j - 1], 9)
        for i, r in enumerate(data, start=2):
            fill = PatternFill("solid", fgColor=r[fill_col - 1].lstrip("#").upper())
            for j in range(1, len(head) + 1):
                c = ws.cell(row=i, column=j)
                c.font, c.fill, c.border = font, fill, Border(bottom=hair)
                c.alignment = Alignment(vertical="top", wrap_text=j in wrap_cols)
        ws.row_dimensions[1].height = 30
        ws.freeze_panes = ws.cell(row=2, column=3)
        ws.auto_filter.ref = f"A1:{get_column_letter(len(head))}{len(data) + 1}"

    ws = wb.active
    ws.title = "Folds"
    wide = {"Fold (catalytic domain)": 46, "CATH": 22, "CATH class": 13, "Fill": 9,
            "Components": 12, "Sequences": 11, "With PDB structure": 11,
            "Curated seeds": 38, "Accessory domains (records)": 70,
            "Homology check": 30, "Weakest TM link": 10, "Backbone branches kept": 12,
            "Bacterial": 10, "Eukaryotic": 10, "Metagenomic": 12}
    sheet(ws, fold_head, fold_rows, wide,
          wrap_cols={fold_head.index(h) + 1 for h in
                     ("Curated seeds", "Accessory domains (records)", "CATH")})
    last = len(fold_rows) + 1
    tot = last + 1
    ws.cell(row=tot, column=2, value="Total").font = bold
    first_n = fold_head.index("Components") + 1
    last_n = fold_head.index("With PDB structure") + 1
    for j in range(first_n, last_n + 1):
        L = get_column_letter(j)
        c = ws.cell(row=tot, column=j, value=f"=SUM({L}2:{L}{last})")
        c.font, c.border = bold, Border(top=thin)
    ws.cell(row=tot, column=2).border = Border(top=thin)

    ws2 = wb.create_sheet("Domains")
    wide2 = {"CATH": 14, "Domain": 52, "CATH class": 13, "Role": 11,
             "Fold it occurs in": 46, "Slice": 7, "Fill": 9, "Components": 12,
             "Sequences": 11, "Share of the fold's sequences": 14,
             "Bacterial": 10, "Eukaryotic": 10, "Metagenomic": 12}
    sheet(ws2, dom_head, dom_rows, wide2)
    seq_col = get_column_letter(dom_head.index("Sequences") + 1)
    slice_col = get_column_letter(dom_head.index("Slice") + 1)
    share = dom_head.index("Share of the fold's sequences") + 1
    f_seq = get_column_letter(fold_head.index("Sequences") + 1)
    for i in range(2, len(dom_rows) + 2):
        c = ws2.cell(row=i, column=share,
                     value=f"=IFERROR({seq_col}{i}/INDEX(Folds!${f_seq}$2:${f_seq}${last},"
                           f"MATCH({slice_col}{i},Folds!$A$2:$A${last},0)),\"\")")
        c.number_format = "0%"

    notes = wb.create_sheet("Notes")
    lines = [
        ("What this is", "CATH domains of the plasticome v2 union set, grouped by the "
                         "fold groups of the tree figure (out/figure_pie.svg). Built by "
                         "tree/09_cath_table.py."),
        ("Source", "catalytic-cath, accessory-cath and domain columns of "
                   "plasticome.v2 - plasticome-v1-v260701-union-annotated (1).csv, "
                   f"609 records. This workbook covers the set '{SET}': "
                   f"{len(rows)} records."),
        ("Fold / slice", "One slice of the pie. Components (connected components of the "
                         "30 percent identity SSN) whose catalytic domain is the same CATH "
                         "superfamily and whose homology was confirmed by structure."),
        ("Row colour", "The fold's background tint in the pie. The hex value is in the "
                       "Fill column. Folds that are a single sequence share one grey."),
        ("Sequences", f"Distinct sequences, as in the trees. {len(tips)} in total."
                      + (" One 28 residue fragment (Arylacylamidase P80008.1) is in "
                         "Records but not in Sequences, because it was left out of "
                         "the trees." if SET == "all" else
                         " Only the 90 percent centroids are included.")),
        ("Records", "Rows of the union table. A protein is listed twice when both PAZy "
                    "pulls carry it, so Records is larger than Sequences."),
        ("Substrate columns", "Distinct sequences per substrate. Combi. is a sequence "
                              "with more than one substrate. NR is not reported."),
        ("Origin columns", "Distinct sequences: bacterial (includes archaeal), eukaryotic, "
                           "metagenomic."),
        ("With PDB structure", "Distinct sequences with has_structure = yes in the union "
                               "table."),
        ("Accessory domains", "Domains beside the catalytic one, with the number of "
                              "records carrying each. On the Domains sheet an accessory "
                              "domain has one row per fold it occurs in."),
        ("Homology check", "From step 4. A fold with several components is confirmed "
                           "when every component reaches TM-score 0.5 with another member. "
                           "Weakest TM link is the lowest of those best scores."),
        ("Backbone branches kept", "From steps 5 and 6. Internal branches between "
                                   "components with UFBoot 95 or more that an independent "
                                   "tree also recovers, out of all internal branches."),
        ("No CATH id", "The oxygenase MpaB family has a domain name but no CATH "
                       "accession in the union table, shown as none."),
        ("CATH class", "First number of the accession: 1 mainly alpha, 2 mainly beta, "
                       "3 alpha beta, 4 few secondary structures, 6 special."),
        ("Totals", "The Total row on Folds and the share column on Domains are formulas."),
    ]
    for i, (k, v) in enumerate(lines, start=1):
        a, b = notes.cell(row=i, column=1, value=k), notes.cell(row=i, column=2, value=v)
        a.font, b.font = bold, font
        a.alignment = Alignment(vertical="top")
        b.alignment = Alignment(wrap_text=True, vertical="top")
    notes.column_dimensions["A"].width = 24
    notes.column_dimensions["B"].width = 110

    path = os.path.join(OUT, "cath_domains.xlsx")
    wb.save(path)
    n_dom = len({(r[0], r[1]) for r in dom_rows if r[3] != "unassigned"})
    print(f"{len(fold_rows)} folds; {n_dom} distinct CATH domains in "
          f"{len(dom_rows)} domain by fold rows "
          f"({sum(1 for r in dom_rows if r[3] == 'accessory')} accessory)")
    print(f"records {sum(r[7] for r in fold_rows)}, sequences {sum(r[6] for r in fold_rows)}, "
          f"components {sum(r[5] for r in fold_rows)}")


if __name__ == "__main__":
    main()
