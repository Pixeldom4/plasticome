#!/usr/bin/env python3
"""
Step 2: do accessory domains or source metadata track measured activity?

A first, simple probe, as the feedback of 2026-10-27 asked: "using activity
datasets from large studies, see if correlation exists between accessory
domains and activity, and the same for metadata".

Activity. PETai labels (petadex-adjacency-landscape/data/petai_labels.csv),
13,092 measurements on 1,219 sequences from 696 references. Units differ from
study to study, so every measurement is turned into a percentile inside its
own group: same reference, variable, substrate and unit, groups of 5 or more
sequences only. Only variables that rise with activity are used (degradation
proxies, kcat, kcat/Km, specific and relative activity, halo size); Km and
melting temperature are left out of the activity score. A sequence's score on
a substrate class is the mean of its percentiles there.

Substrate class. "plastic" is a polymer (PET and gfPET, PBAT, PLA, PBSA, PBS,
PCL, PHB, PU and the like); "model" is a small ester or monomer (pNP esters,
BHET, MHET, tributyrin and the like). This is the known substrate against
plastic split the story asks for, measured instead of assumed.

Accessory domains. The step 1 Pfam scan of the measured sequences themselves.
The highest scoring hit is the catalytic domain; any other hit that overlaps it
by 30 percent of its own length or less is accessory.

Metadata. Temperature, pH and salinity of the source BioSamples of the atlas
family the sequence maps to (petai_families.csv to family_metadata.csv), and
the family's ecosystem. Family level, not sequence level.

Writes out/figure_probe.{png,pdf,svg}, out/probe_summary.tsv and
out/sequence_table.tsv.
"""

import collections
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
LAND = os.path.join(HERE, "..", "..", "..", "petadex-adjacency-landscape", "data")
OUT = os.path.join(HERE, "out")
INK, MUTED, SURFACE = "#1a1a1a", "#6b6b6b", "#fcfcfb"
WITH, WITHOUT = "#c0632b", "#4f7fb3"

ACTIVITY = re.compile(r"proxy for plastic degradation|^kcat$|kcat/km|specific activity|"
                      r"relative activity|activity|halo|weight loss|released|conversion", re.I)
NOT_ACTIVITY = re.compile(r"^km$|melting|tm\b|temperature|half-life|expression|yield", re.I)
PLASTIC = re.compile(r"^(gf)?pet|pbat|pla|pbs|pcl|phb|pha|pu\b|pur|polyure|nylon|pa6|pe\b|"
                     r"film|powder|polymer|impranil|pes|ppc|pc\b|pbsa", re.I)
MODEL = re.compile(r"pnp|p-np|nitrophenyl|bhet|mhet|tributyrin|ester|acetate|butyrate|"
                   r"tpa|eg\b|acetylcholine|olive|tween|trioleate", re.I)


def substrate_class(s):
    s = str(s).strip()
    if MODEL.search(s):
        return "model"
    if PLASTIC.search(s):
        return "plastic"
    return None


def accessory(pf):
    """seq_md5 -> (catalytic pfam name, [accessory pfam names])."""
    out = {}
    for sid, g in pf.groupby("seq_md5"):
        g = g.sort_values("score", ascending=False)
        top = g.iloc[0]
        acc = []
        for _, h in g.iloc[1:].iterrows():
            ov = max(0, min(h.end, top.end) - max(h.start, top.start) + 1)
            if ov <= 0.3 * (h.end - h.start + 1):
                acc.append(h.pfam_name)
        out[sid] = (top.pfam_name, sorted(set(acc)))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    lab = pd.read_csv(os.path.join(LAND, "petai_labels.csv"))
    pf = pd.read_csv(os.path.join(HERE, "work", "petai_pfam.tsv"), sep="\t")
    fam = pd.read_csv(os.path.join(LAND, "petai_families.csv"))
    meta = pd.read_csv(os.path.join(LAND, "family_metadata.csv"))
    ann = pd.read_csv(os.path.join(LAND, "family_annotations.csv"))
    acc = accessory(pf)

    # --- activity score ---------------------------------------------------
    q = lab[(lab.task == "quantitative") & lab.variable.str.contains(ACTIVITY)
            & ~lab.variable.str.contains(NOT_ACTIVITY)].copy()
    q["cls"] = q.substrate.map(substrate_class)
    q = q[q.cls.notna()]
    q["grp"] = q.ref.astype(str) + "|" + q.variable + "|" + q.substrate.astype(str) + "|" + q.unit.astype(str)
    per = q.groupby(["grp", "seq_md5"]).value.mean().reset_index()
    size = per.groupby("grp").seq_md5.nunique()
    per = per[per.grp.map(size) >= 5]
    per["pct"] = per.groupby("grp").value.rank(pct=True)
    per = per.merge(q[["grp", "cls", "ref"]].drop_duplicates("grp"), on="grp")
    score = per.groupby(["seq_md5", "cls"]).pct.mean().reset_index()

    # --- per sequence table -------------------------------------------------
    seqs = pd.DataFrame({"seq_md5": sorted(set(lab.seq_md5))})
    seqs["catalytic"] = seqs.seq_md5.map(lambda s: acc.get(s, (None, []))[0])
    seqs["accessory"] = seqs.seq_md5.map(lambda s: ";".join(acc.get(s, (None, []))[1]))
    seqs["has_accessory"] = seqs.accessory != ""
    seqs["length"] = seqs.seq_md5.map(pf.drop_duplicates("seq_md5").set_index("seq_md5").length)
    f = fam.dropna(subset=["family_id"]).drop_duplicates("seq_md5").set_index("seq_md5").family_id
    seqs["family_id"] = seqs.seq_md5.map(f)
    m = meta.drop_duplicates("family_id").set_index("family_id")
    a = ann.drop_duplicates("family_id").set_index("family_id")
    for col in ("temperature_c", "ph", "salinity_psu"):
        seqs[col] = seqs.family_id.map(m[col])
    seqs["ecosystem"] = seqs.family_id.map(a["ecosystem"])
    tm = lab[lab.variable.str.contains("Melting temperature", case=False)].groupby("seq_md5").value.median()
    seqs["tm_measured"] = seqs.seq_md5.map(tm)
    for cls in ("plastic", "model"):
        s = score[score.cls == cls].set_index("seq_md5").pct
        seqs[f"activity_{cls}"] = seqs.seq_md5.map(s)
    halo = lab[lab.ref.astype(str).str.startswith("halo:")]
    seqs["halo_active"] = seqs.seq_md5.map(~halo.groupby("seq_md5").is_negative.max().astype(bool))
    seqs.to_csv(os.path.join(OUT, "sequence_table.tsv"), sep="\t", index=False)

    # --- tests ----------------------------------------------------------------
    rows = []

    def mw(cls):
        x = seqs[seqs.has_accessory][f"activity_{cls}"].dropna()
        y = seqs[~seqs.has_accessory][f"activity_{cls}"].dropna()
        p = stats.mannwhitneyu(x, y).pvalue if len(x) > 2 and len(y) > 2 else np.nan
        rows.append([f"accessory vs none, {cls} substrates", len(x), len(y),
                     f"{x.median():.2f}", f"{y.median():.2f}", f"{p:.3g}"])
        return x, y, p

    # paired inside references: a reference that measured both kinds
    def within_ref(cls):
        d = per[per.cls == cls].copy()
        d["has"] = d.seq_md5.map(seqs.set_index("seq_md5").has_accessory)
        diffs = []
        for ref, g in d.groupby("ref"):
            if g.has.nunique() == 2 and g.seq_md5.nunique() >= 5:
                diffs.append(g[g.has].pct.mean() - g[~g.has].pct.mean())
        if len(diffs) >= 3:
            p = stats.wilcoxon(diffs).pvalue
        else:
            p = np.nan
        rows.append([f"within reference, {cls}: references with both", len(diffs), "",
                     f"{np.median(diffs):.2f}" if diffs else "", "", f"{p:.3g}"])
        return diffs, p

    res = {cls: mw(cls) for cls in ("plastic", "model")}
    wr = {cls: within_ref(cls) for cls in ("plastic", "model")}
    corr = []
    for col, label in (("temperature_c", "source temperature"), ("ph", "source pH"),
                       ("salinity_psu", "source salinity")):
        for cls in ("plastic", "model"):
            d = seqs[[col, f"activity_{cls}"]].dropna()
            r, p = stats.spearmanr(d[col], d[f"activity_{cls}"]) if len(d) > 5 else (np.nan, np.nan)
            corr.append((label, cls, len(d), r, p))
            rows.append([f"spearman {label} vs activity, {cls}", len(d), "", f"{r:.2f}", "", f"{p:.3g}"])
    d = seqs[["temperature_c", "tm_measured"]].dropna()
    rt, pt = stats.spearmanr(d.temperature_c, d.tm_measured)
    rows.append(["spearman source temperature vs measured Tm", len(d), "", f"{rt:.2f}", "", f"{pt:.3g}"])
    h = seqs.dropna(subset=["halo_active"])
    ct = pd.crosstab(h.has_accessory, h.halo_active)
    if ct.shape == (2, 2):
        o, ph = stats.fisher_exact(ct.values)
        rows.append(["halo assay active, accessory vs none (Fisher)", int(ct.loc[True].sum()),
                     int(ct.loc[False].sum()), f"{ct.loc[True, True] / ct.loc[True].sum():.2f}",
                     f"{ct.loc[False, True] / ct.loc[False].sum():.2f}", f"{ph:.3g}"])
    pd.DataFrame(rows, columns=["test", "n_a", "n_b", "stat_a", "stat_b", "p"]).to_csv(
        os.path.join(OUT, "probe_summary.tsv"), sep="\t", index=False)

    # --- figure -----------------------------------------------------------------
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.6))
    fig.patch.set_facecolor(SURFACE)
    a0 = ax[0, 0]
    for i, cls in enumerate(("model", "plastic")):
        x, y, p = res[cls]
        bp = a0.boxplot([y, x], positions=[i * 3, i * 3 + 1], widths=0.7, patch_artist=True,
                        showfliers=False)
        for patch, c in zip(bp["boxes"], (WITHOUT, WITH)):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        for j, v in enumerate((y, x)):
            a0.scatter(np.random.default_rng(1).normal(i * 3 + j, 0.08, len(v)), v, s=4,
                       c=INK, alpha=0.35, lw=0)
        a0.text(i * 3 + 0.5, 1.05, f"p = {p:.2g}\nn = {len(y)} / {len(x)}", ha="center",
                fontsize=7, color=MUTED)
    a0.set_xticks([0, 1, 3, 4])
    a0.set_xticklabels(["no acc.", "accessory", "no acc.", "accessory"], fontsize=7)
    a0.text(0.5, -0.17, "model ester substrates", ha="center", fontsize=8)
    a0.text(3.5, -0.17, "plastic substrates", ha="center", fontsize=8)
    a0.set_ylim(-0.05, 1.22)
    a0.set_ylabel("activity, percentile inside its study", fontsize=8)
    a0.set_title("A  Accessory domain present or not", loc="left", fontsize=9.5)

    a1 = ax[0, 1]
    cnt = collections.Counter(d for s in seqs.accessory if s for d in s.split(";"))
    top = [d for d, _ in cnt.most_common(10)]
    ys = []
    for i, dname in enumerate(top):
        has = seqs.accessory.str.split(";").map(lambda L: dname in L)
        for k, cls in enumerate(("model", "plastic")):
            v = seqs[has][f"activity_{cls}"].dropna()
            if len(v):
                a1.errorbar(v.median(), i + (k - 0.5) * 0.3, xerr=[[v.median() - v.quantile(.25)],
                            [v.quantile(.75) - v.median()]], fmt="o", ms=4,
                            color=("#8a8a8a", INK)[k], lw=0.8)
                a1.text(1.02, i + (k - 0.5) * 0.3, f"{len(v)}", fontsize=6, va="center",
                        color=("#8a8a8a", INK)[k], transform=a1.get_yaxis_transform())
    base = {cls: seqs[~seqs.has_accessory][f"activity_{cls}"].median() for cls in ("model", "plastic")}
    a1.axvline(base["plastic"], color=INK, lw=0.6, ls=(0, (3, 2)))
    a1.axvline(base["model"], color="#8a8a8a", lw=0.6, ls=(0, (3, 2)))
    a1.set_yticks(range(len(top)))
    a1.set_yticklabels([f"{d} ({cnt[d]})" for d in top], fontsize=6.8)
    a1.invert_yaxis()
    a1.set_xlim(0, 1)
    a1.set_xlabel("activity percentile, median and IQR", fontsize=8)
    a1.set_title("B  Each accessory domain (black plastic, grey model ester)\n"
                 "dashed: sequences without accessory domains", loc="left", fontsize=9)

    a2 = ax[0, 2]
    for cls, c in (("model", "#8a8a8a"), ("plastic", INK)):
        diffs, p = wr[cls]
        a2.scatter(diffs, np.full(len(diffs), 0 if cls == "model" else 1) +
                   np.random.default_rng(2).normal(0, 0.05, len(diffs)), s=14, c=c, alpha=0.7, lw=0)
        a2.text(0.98, 0 if cls == "model" else 1, f"{len(diffs)} studies, p = {p:.2g}",
                transform=a2.get_yaxis_transform(), ha="right", va="bottom", fontsize=7, color=c)
    a2.axvline(0, color=MUTED, lw=0.6)
    a2.set_yticks([0, 1])
    a2.set_yticklabels(["model ester", "plastic"], fontsize=8)
    a2.set_ylim(-0.6, 1.6)
    a2.set_xlabel("accessory minus none, mean percentile inside one study", fontsize=8)
    a2.set_title("C  Same comparison inside each study", loc="left", fontsize=9.5)

    a3 = ax[1, 0]
    a3.scatter(d.temperature_c, d.tm_measured, s=10, c=INK, alpha=0.5, lw=0)
    a3.set_xlabel("source BioSample temperature of the family (C)", fontsize=8)
    a3.set_ylabel("measured melting temperature (C)", fontsize=8)
    a3.set_title(f"D  Source temperature vs measured Tm\nSpearman {rt:.2f}, p = {pt:.2g}, n = {len(d)}",
                 loc="left", fontsize=9)

    a4 = ax[1, 1]
    labs = [f"{l}\n{c}" for l, c, *_ in corr]
    rs = [r for *_, r, p in corr]
    cols = [INK if c == "plastic" else "#8a8a8a" for _, c, *_ in corr]
    a4.barh(range(len(corr)), rs, color=cols)
    for i, (l, c, n, r, p) in enumerate(corr):
        a4.text(r + (0.01 if r >= 0 else -0.01), i, f"n={n}, p={p:.2g}", va="center",
                ha="left" if r >= 0 else "right", fontsize=6.5, color=MUTED)
    a4.set_yticks(range(len(corr)))
    a4.set_yticklabels(labs, fontsize=6.8)
    a4.axvline(0, color=MUTED, lw=0.6)
    a4.set_xlim(-0.5, 0.5)
    a4.invert_yaxis()
    a4.set_xlabel("Spearman correlation with activity percentile", fontsize=8)
    a4.set_title("E  Source metadata against activity", loc="left", fontsize=9.5)

    a5 = ax[1, 2]
    eco = seqs.dropna(subset=["ecosystem", "activity_plastic"])
    top_e = eco.ecosystem.value_counts()
    top_e = top_e[top_e >= 8].index[:8]
    data = [eco[eco.ecosystem == e].activity_plastic for e in top_e]
    if data:
        a5.boxplot(data, orientation="horizontal", widths=0.6, showfliers=False)
        a5.set_yticks(range(1, len(top_e) + 1))
        a5.set_yticklabels([f"{e} ({len(x)})" for e, x in zip(top_e, data)], fontsize=7)
        if len(data) > 1:
            pk = stats.kruskal(*data).pvalue
            a5.set_title(f"F  Ecosystem of the family, plastic activity\nKruskal-Wallis p = {pk:.2g}",
                         loc="left", fontsize=9)
            rows.append(["kruskal ecosystem vs plastic activity", sum(len(x) for x in data), "",
                         "", "", f"{pk:.3g}"])
    a5.set_xlabel("activity percentile on plastic", fontsize=8)
    for a in ax.ravel():
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(labelsize=7)
        a.set_facecolor(SURFACE)
    fig.suptitle("Probe: accessory domains and source metadata against measured activity (PETai labels)",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(OUT, f"figure_probe.{ext}"), dpi=250, facecolor=SURFACE)
    pd.DataFrame(rows, columns=["test", "n_a", "n_b", "stat_a", "stat_b", "p"]).to_csv(
        os.path.join(OUT, "probe_summary.tsv"), sep="\t", index=False)

    print(f"{len(seqs)} sequences; {seqs.has_accessory.sum()} with an accessory domain; "
          f"scored on plastic {seqs.activity_plastic.notna().sum()}, model {seqs.activity_model.notna().sum()}")
    print("catalytic:", collections.Counter(seqs.catalytic).most_common(8))
    print("accessory:", cnt.most_common(12))
    for r in rows:
        print("  ", r)


if __name__ == "__main__":
    main()
