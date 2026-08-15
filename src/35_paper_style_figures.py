#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
35 - Paper-style figure gallery (Nature-Metabolism layout, our immune data)
===========================================================================

Reproduces the *figure TYPES* of You et al., Nature Metabolism 2025
("Mapping the plasma metabolome to human health and disease", 10 main figures)
-- but built entirely from OUR plasma-immunome causal atlas, using ALL immune
proteins we can legitimately instrument (673 cis-eQTL-instrumented immune genes
x 28 FinnGen R12 diseases = 18,844 Mendelian-randomization tests).

It is purely additive: it reads existing result tables, writes only new PNGs
into 08_figures/paper_style/, and never modifies any existing data or figure.

Panel  <- their analogue
  P1  study-overview infographic          (their Fig 1)
  P2  protein-disease MR mirrored Manhattan by disease chapter (their Fig 2)
  P6  atlas discrimination / evidence     (their Fig 6 MetRS)
  P7  causal proteins per disease category, stacked by immune class (their Fig 7)
  P8  colocalization bubble map           (their Fig 8)
  P9  prioritised-target summary          (their Fig 9)

Run:
    python src/35_paper_style_figures.py
"""

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GC   = os.path.join(ROOT, "06_genetic_causality")
OUT  = os.path.join(ROOT, "08_figures", "paper_style")
os.makedirs(OUT, exist_ok=True)

FDR_SIG = 0.05

# ---- Nature-ish styling ---------------------------------------------------- #
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333",
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "legend.frameon": False,
})

# disease-category palette (mirrors the ICD-chapter colouring of the paper)
CAT_COLORS = {
    "Autoimmune":     "#D1495B",
    "Cardiovascular": "#2E6F95",
    "Metabolic":      "#E9963A",
    "Neuro/Aging":    "#6A4C93",
    "Renal":          "#2A9D8F",
}
CLASS_COLORS = {
    "Immune-cell-enriched":     "#4C72B0",
    "Cytokine/Interleukin":     "#DD8452",
    "CD/Leukocyte surface":     "#55A868",
    "TNF superfamily":          "#C44E52",
    "Complement/Coagulation":   "#8172B3",
    "Immune checkpoint":        "#937860",
    "Chemokine axis":           "#DA8BC3",
    "Ig/B-cell/Fc":             "#8C8C8C",
    "Interferon axis":          "#CCB974",
    "HLA/Antigen presentation": "#64B5CD",
    "Acute-phase":              "#B07AA1",
}
def class_color(c):
    return CLASS_COLORS.get(c, "#BBBBBB")


def load():
    mr = pd.read_csv(os.path.join(GC, "cis_MR_phenome_results.tsv"), sep="\t")
    co = pd.read_csv(os.path.join(GC, "coloc_phenome_results.tsv"), sep="\t")
    nov = pd.read_csv(os.path.join(GC, "novelty_engine_ranked.tsv"), sep="\t")
    return mr, co, nov


# --------------------------------------------------------------------------- #
def panel1_overview(mr):
    """Study-overview infographic (their Fig 1)."""
    n_prot = 1007
    n_inst = mr["gene_symbol"].nunique()
    n_dis  = mr["disease"].nunique()
    n_test = len(mr)
    n_hit  = int((mr["FDR"] < FDR_SIG).sum())
    n_cat  = mr["disease_category"].nunique()

    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.axis("off")
    ax.set_xlim(0, 100); ax.set_ylim(0, 46)
    ax.text(50, 43, "An open, genetics-anchored plasma immunome\u2013phenome atlas",
            ha="center", va="center", fontsize=15, fontweight="bold")
    ax.text(50, 38.6, "the same immune cis-eQTL instruments tested across the disease phenome",
            ha="center", va="center", fontsize=10, style="italic", color="#555555")

    stages = [
        ("Plasma immune\nproteome", f"{n_prot:,} curated\nOlink immune proteins", "#4C72B0"),
        ("Genetic\ninstruments", f"{n_inst} proteins with a\ncis-eQTL instrument\n(eQTLGen)", "#55A868"),
        ("Disease\nphenome", f"{n_dis} FinnGen R12\ndiseases in {n_cat} categories", "#E9963A"),
        ("Mendelian\nrandomization", f"{n_test:,} cis-MR tests\n{n_hit} hits at FDR<0.05", "#D1495B"),
        ("Colocalization\n+ pQTL + replication", "transcript & protein\ncausal evidence tiers", "#6A4C93"),
    ]
    w, gap = 15.5, 3.0
    x0 = (100 - (len(stages) * w + (len(stages) - 1) * gap)) / 2
    cy = 20
    for i, (title, sub, col) in enumerate(stages):
        x = x0 + i * (w + gap)
        box = FancyBboxPatch((x, cy - 8), w, 16, boxstyle="round,pad=0.4,rounding_size=1.4",
                             linewidth=1.4, edgecolor=col, facecolor=col + "22")
        ax.add_patch(box)
        ax.text(x + w / 2, cy + 4.6, title, ha="center", va="center",
                fontsize=10.5, fontweight="bold", color=col)
        ax.text(x + w / 2, cy - 2.6, sub, ha="center", va="center", fontsize=8.4, color="#333333")
        if i < len(stages) - 1:
            ax.add_patch(FancyArrowPatch((x + w + 0.2, cy), (x + w + gap - 0.2, cy),
                         arrowstyle="-|>", mutation_scale=16, linewidth=1.6, color="#888888"))
    ax.text(50, 3.2,
            "Every claim bound to an evidence tier (T2 MHC-caution \u2192 T5 protein-level causal); "
            "no individual-level or gated data used.",
            ha="center", va="center", fontsize=8.6, style="italic", color="#666666")
    fig.savefig(os.path.join(OUT, "P1_study_overview.png"))
    plt.close(fig)
    print("  P1_study_overview.png")


# --------------------------------------------------------------------------- #
def panel2_manhattan(mr):
    """Mirrored Manhattan of protein-disease MR, grouped/coloured by category
    (their Fig 2a,b). Risk (OR>1) above axis, protective (OR<1) below."""
    d = mr.copy()
    d = d[np.isfinite(d["MR_p"]) & (d["MR_p"] > 0)]
    d["signed_logp"] = -np.log10(d["MR_p"]) * np.where(d["OR"] >= 1, 1, -1)

    cats = ["Autoimmune", "Cardiovascular", "Metabolic", "Neuro/Aging", "Renal"]
    cats = [c for c in cats if c in set(d["disease_category"])]
    # x layout: diseases grouped by category
    order, xpos, xticks, xlabels, boundaries = [], {}, [], [], []
    x = 0
    rng = np.random.default_rng(7)
    for c in cats:
        dis = sorted(d.loc[d["disease_category"] == c, "disease"].unique())
        for dd in dis:
            xpos[dd] = x; x += 1
        boundaries.append((c, x - len(dis), x - 1))
    d["x"] = d["disease"].map(xpos) + rng.uniform(-0.32, 0.32, len(d))

    thr = -np.log10(FDR_SIG)
    fig, ax = plt.subplots(figsize=(13, 6))
    for c in cats:
        sub = d[d["disease_category"] == c]
        sig = sub["FDR"] < FDR_SIG
        ax.scatter(sub.loc[~sig, "x"], sub.loc[~sig, "signed_logp"], s=5,
                   color=CAT_COLORS[c], alpha=0.20, linewidths=0, rasterized=True)
        ax.scatter(sub.loc[sig, "x"], sub.loc[sig, "signed_logp"], s=20,
                   color=CAT_COLORS[c], alpha=0.9, linewidths=0.3, edgecolor="white")
    ax.axhline(0, color="#333333", linewidth=0.9)
    ax.axhline(thr, color="grey", linewidth=0.8, linestyle="--")
    ax.axhline(-thr, color="grey", linewidth=0.8, linestyle="--")

    # label the strongest hits
    hits = d[d["FDR"] < FDR_SIG].copy()
    hits["abs"] = hits["signed_logp"].abs()
    for _, r in hits.sort_values("abs", ascending=False).head(22).iterrows():
        ax.annotate(r["gene_symbol"], (r["x"], r["signed_logp"]),
                    fontsize=6.8, ha="center",
                    va="bottom" if r["signed_logp"] > 0 else "top",
                    color="#222222")
    # category bands + labels (reserve headroom so labels never touch the title)
    ymin, ymax = ax.get_ylim()
    band = ymax + (ymax - ymin) * 0.16
    ax.set_ylim(ymin, band)
    lab_y = ymax + (ymax - ymin) * 0.05
    for c, a, b in boundaries:
        ax.axvspan(a - 0.5, b + 0.5, color=CAT_COLORS[c], alpha=0.05, zorder=0)
        ax.text((a + b) / 2, lab_y, c, ha="center", va="bottom",
                fontsize=9.5, fontweight="bold", color=CAT_COLORS[c])
    ax.set_xticks([xpos[dd] for dd in xpos])
    ax.set_xticklabels(list(xpos.keys()), rotation=90, fontsize=6.2)
    ax.set_ylabel("signed \u2212log$_{10}$($P_{MR}$)   \u2190 protective | risk \u2192")
    ax.set_xlim(-1, x)
    ax.set_title("Plasma-immune protein \u2192 disease Mendelian randomization across the phenome",
                 loc="left", pad=14)
    ax.text(0.005, 0.03, f"{len(d):,} cis-MR tests \u00b7 {int((d['FDR']<FDR_SIG).sum())} at FDR<0.05",
            transform=ax.transAxes, fontsize=8, color="#555555")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "P2_mr_manhattan.png"))
    plt.close(fig)
    print("  P2_mr_manhattan.png")


# --------------------------------------------------------------------------- #
def panel7_category_stack(mr):
    """Causal proteins per disease category, stacked by immune class (Fig 7a)."""
    hits = mr[mr["FDR"] < FDR_SIG].copy()
    cats = ["Autoimmune", "Cardiovascular", "Metabolic", "Neuro/Aging", "Renal"]
    cats = [c for c in cats if c in set(hits["disease_category"])]
    classes = [c for c in CLASS_COLORS if c in set(hits["immune_class"])]
    tab = (hits.groupby(["disease_category", "immune_class"]).size()
                .unstack(fill_value=0).reindex(index=cats, columns=classes, fill_value=0))

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    bottom = np.zeros(len(cats))
    for cl in classes:
        vals = tab[cl].values
        ax.bar(range(len(cats)), vals, bottom=bottom, width=0.68,
               color=class_color(cl), edgecolor="white", linewidth=0.5, label=cl)
        bottom += vals
    for i, tot in enumerate(bottom):
        ax.text(i, tot + 0.6, int(tot), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, fontsize=9)
    ax.set_ylabel("causal gene\u2013disease hits (FDR<0.05)")
    ax.set_title("Causal plasma-immune proteins by disease category", loc="left")
    ax.legend(fontsize=7.2, ncol=1, loc="upper right", title="immune class",
              title_fontsize=7.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "P7_category_stack.png"))
    plt.close(fig)
    print("  P7_category_stack.png")


# --------------------------------------------------------------------------- #
def panel8_coloc_bubble(co):
    """Colocalization bubble map: gene x disease, size/colour by PP.H4 (Fig 8)."""
    d = co.copy()
    d = d[d["PP_H4"] >= 0.5]
    if not len(d):
        print("  (no coloc PP.H4>=0.5 -> P8 skipped)"); return
    genes = (d.groupby("gene")["PP_H4"].max().sort_values(ascending=False).index.tolist())
    dis = sorted(d["disease"].unique())
    gi = {g: i for i, g in enumerate(genes)}
    di = {s: j for j, s in enumerate(dis)}

    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(dis) + 3), max(5, 0.32 * len(genes) + 2)))
    sc = ax.scatter([di[s] for s in d["disease"]], [gi[g] for g in d["gene"]],
                    s=30 + (d["PP_H4"] ** 3) * 320, c=d["PP_H4"], cmap="viridis",
                    vmin=0.5, vmax=1.0, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(dis))); ax.set_xticklabels(dis, rotation=90, fontsize=7)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=7)
    ax.set_ylim(-1, len(genes)); ax.set_xlim(-1, len(dis))
    ax.set_title("Transcript colocalization (coloc.abf) of immune protein \u2194 disease", loc="left")
    cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("PP.H4 (shared causal variant)")
    ax.grid(True, color="#EEEEEE", linewidth=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "P8_coloc_bubble.png"))
    plt.close(fig)
    print("  P8_coloc_bubble.png")


# --------------------------------------------------------------------------- #
def panel9_priority(nov):
    """Prioritised-target summary (their Fig 9): top novelty-priority targets."""
    d = nov.copy().sort_values("novelty_priority", ascending=False).head(25).iloc[::-1]
    labels = d["gene_symbol"] + "  \u2192  " + d["disease"]
    colors = [CAT_COLORS.get(c, "#999999") for c in d["disease_category"]]
    fig, ax = plt.subplots(figsize=(8.8, 7.2))
    y = np.arange(len(d))
    ax.barh(y, d["novelty_priority"], color=colors, edgecolor="white", linewidth=0.5)
    for i, (_, r) in enumerate(d.iterrows()):
        marks = []
        if r.get("PP_H4", 0) >= 0.8: marks.append("coloc")
        if r.get("druggability", 0) and r["druggability"] > 0: marks.append("druggable")
        if marks:
            ax.text(r["novelty_priority"] + 0.03, i, " \u00b7 ".join(marks),
                    va="center", fontsize=6.6, color="#444444")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("novelty-priority score")
    ax.set_title("Top prioritised novel plasma-immune targets", loc="left")
    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                      markerfacecolor=CAT_COLORS[c], markeredgecolor="white", label=c)
               for c in CAT_COLORS if c in set(d["disease_category"])]
    ax.legend(handles=handles, fontsize=7.4, loc="lower right", title="disease category",
              title_fontsize=7.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "P9_priority_summary.png"))
    plt.close(fig)
    print("  P9_priority_summary.png")


# --------------------------------------------------------------------------- #
def panel6_atlas_perf(mr, nov):
    """Atlas discrimination panel (their Fig 6 MetRS analogue): per-disease
    causal yield + directional split. PIRS is unshipped, so this shows the
    causal atlas signal honestly rather than a fabricated risk-score AUC."""
    hits = mr[mr["FDR"] < FDR_SIG].copy()
    g = hits.groupby("disease")
    tab = pd.DataFrame({
        "risk": g.apply(lambda x: int((x["OR"] > 1).sum())),
        "protective": g.apply(lambda x: int((x["OR"] < 1).sum())),
    })
    tab["cat"] = hits.groupby("disease")["disease_category"].first()
    tab["tot"] = tab["risk"] + tab["protective"]
    tab = tab.sort_values(["cat", "tot"], ascending=[True, False])

    fig, ax = plt.subplots(figsize=(9, 6.4))
    y = np.arange(len(tab))
    ax.barh(y, tab["risk"], color="#C44E52", edgecolor="white", linewidth=0.4, label="risk (OR>1)")
    ax.barh(y, -tab["protective"], color="#4C72B0", edgecolor="white", linewidth=0.4,
            label="protective (OR<1)")
    ax.axvline(0, color="#333333", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d}" for d in tab.index], fontsize=7.4)
    # colour disease tick labels by category
    for tick, c in zip(ax.get_yticklabels(), tab["cat"]):
        tick.set_color(CAT_COLORS.get(c, "#333333"))
    ax.set_xlabel("protective  \u2190   causal immune proteins   \u2192  risk")
    ax.set_title("Per-disease causal immune-protein yield and direction", loc="left")
    ax.legend(fontsize=8, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "P6_atlas_direction.png"))
    plt.close(fig)
    print("  P6_atlas_direction.png")


def main():
    print("[35] Paper-style figure gallery (all immune data) ...")
    mr, co, nov = load()
    print(f"    loaded {len(mr):,} MR tests, {mr['gene_symbol'].nunique()} genes, "
          f"{mr['disease'].nunique()} diseases, {len(co)} coloc rows, {len(nov)} ranked targets")
    panel1_overview(mr)
    panel2_manhattan(mr)
    panel6_atlas_perf(mr, nov)
    panel7_category_stack(mr)
    panel8_coloc_bubble(co)
    panel9_priority(nov)
    print(f"[35] Done -> {OUT}")


if __name__ == "__main__":
    main()
