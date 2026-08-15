#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
38 - Extra paper-style panels (addable-now analogues of the metabolome atlas)
=============================================================================
Adds the remaining reference-paper figure TYPES that CAN be built honestly from
our existing summary tables (no individual-level, imaging or longitudinal data,
which we deliberately never fabricate):

  X1  diseases-per-immune-protein lollipop        (their Fig 3c analogue)
  X2  shared-disease-category distribution         (their Fig 3d analogue)
  X3  prioritised-target score decomposition       (their Fig 6d importance analogue)
  X4  immune cell-source enrichment of causal hits  (their Fig 3e/cell-origin analogue)

Purely additive: reads existing tables, writes new PNGs into
08_figures/paper_style/, modifies nothing.

Run:
    python src/38_extra_paper_panels.py
"""

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GC   = os.path.join(ROOT, "06_genetic_causality")
OUT  = os.path.join(ROOT, "08_figures", "paper_style")
os.makedirs(OUT, exist_ok=True)

FDR_SIG = 0.05

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

CAT_COLORS = {
    "Autoimmune":     "#D1495B",
    "Cardiovascular": "#2E6F95",
    "Metabolic":      "#E9963A",
    "Neuro/Aging":    "#6A4C93",
    "Renal":          "#2A9D8F",
}


# --------------------------------------------------------------------------- #
def panel_x1_diseases_per_protein(mr):
    """X1 (their Fig 3c): how many diseases each immune protein causally drives."""
    hits = mr[mr["FDR"] < FDR_SIG].copy()
    g = hits.groupby("gene_symbol")
    tab = pd.DataFrame({
        "n_dis": g["disease"].nunique(),
        "cat":   g["disease_category"].agg(lambda s: s.value_counts().index[0]),
    }).sort_values("n_dis", ascending=True)
    tab = tab[tab["n_dis"] >= 2]          # show pleiotropic proteins (>=2 diseases)

    fig, ax = plt.subplots(figsize=(7.6, max(4.5, 0.30 * len(tab) + 1.2)))
    y = np.arange(len(tab))
    cols = [CAT_COLORS.get(c, "#999999") for c in tab["cat"]]
    ax.hlines(y, 0, tab["n_dis"], color="#CCCCCC", linewidth=1.3, zorder=1)
    ax.scatter(tab["n_dis"], y, s=70, c=cols, edgecolor="white", linewidth=0.8, zorder=2)
    for i, v in enumerate(tab["n_dis"]):
        ax.text(v + 0.12, i, str(int(v)), va="center", fontsize=8, color="#333333")
    ax.set_yticks(y); ax.set_yticklabels(tab.index, fontsize=8)
    ax.set_xlabel("distinct diseases causally affected (FDR<0.05)")
    ax.set_xlim(0, tab["n_dis"].max() + 1)
    ax.set_title("Pleiotropic immune proteins \u2014 diseases per protein", loc="left")
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=8,
                      markerfacecolor=CAT_COLORS[c], markeredgecolor="white", label=c)
               for c in CAT_COLORS if c in set(tab["cat"])]
    ax.legend(handles=handles, fontsize=7.4, loc="lower right",
              title="dominant category", title_fontsize=7.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "X1_diseases_per_protein.png"))
    plt.close(fig)
    print("  X1_diseases_per_protein.png")


# --------------------------------------------------------------------------- #
def panel_x2_shared_categories(mr):
    """X2 (their Fig 3d): distribution of how many disease CATEGORIES each causal
    immune protein spans (cross-chapter pleiotropy)."""
    hits = mr[mr["FDR"] < FDR_SIG].copy()
    per_gene = hits.groupby("gene_symbol")["disease_category"].nunique()
    dist = per_gene.value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    x = dist.index.values
    palette = ["#BFD7EA", "#7FA8C9", "#2E6F95", "#1B4965"]
    cols = [palette[min(i, len(palette) - 1)] for i in range(len(x))]
    bars = ax.bar(x, dist.values, width=0.62, color=cols, edgecolor="white", linewidth=0.8)
    for b, v in zip(bars, dist.values):
        ax.text(b.get_x() + b.get_width() / 2, v + max(dist.values) * 0.01, str(int(v)),
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xlabel("number of disease categories a protein causally spans")
    ax.set_ylabel("causal immune proteins")
    ax.set_title("Cross-chapter pleiotropy of causal immune proteins", loc="left")
    n_multi = int((per_gene >= 2).sum()); n_tot = int(per_gene.shape[0])
    ax.text(0.98, 0.95, f"{n_multi}/{n_tot} proteins act across \u22652 categories",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#555555")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "X2_shared_category_distribution.png"))
    plt.close(fig)
    print("  X2_shared_category_distribution.png")


# --------------------------------------------------------------------------- #
def panel_x3_score_decomposition(nov):
    """X3 (their Fig 6d): decomposition of the prioritisation score into its
    evidence components for the top targets."""
    comps = [("s_causal", "causal MR", "#D1495B"),
             ("s_coloc",  "colocalization", "#2E6F95"),
             ("s_pleio",  "pleiotropy", "#E9963A"),
             ("s_drug",   "druggability", "#6A4C93"),
             ("s_cell",   "cell-source", "#2A9D8F")]
    d = nov.copy().sort_values("novelty_priority", ascending=False).head(20).iloc[::-1]
    labels = d["gene_symbol"] + "  \u2192  " + d["disease"]

    fig, ax = plt.subplots(figsize=(8.8, 6.6))
    y = np.arange(len(d))
    left = np.zeros(len(d))
    for col, lab, c in comps:
        vals = d[col].values
        ax.barh(y, vals, left=left, color=c, edgecolor="white", linewidth=0.4, label=lab)
        left += vals
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("cumulative evidence score (summed components)")
    ax.set_title("Prioritised-target score decomposition (top 20)", loc="left")
    ax.legend(fontsize=7.6, loc="lower right", ncol=1, title="evidence layer",
              title_fontsize=7.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "X3_score_decomposition.png"))
    plt.close(fig)
    print("  X3_score_decomposition.png")


# --------------------------------------------------------------------------- #
def panel_x4_cellsource(cs):
    """X4 (their Fig 3e / cell-origin): enrichment of causal hits by the immune
    cell that expresses each protein."""
    d = cs.copy().sort_values("odds_ratio", ascending=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    y = np.arange(len(d))
    cols = ["#C44E52" if o > 1 else "#4C72B0" for o in d["odds_ratio"]]
    ax.barh(y, d["odds_ratio"], color=cols, edgecolor="white", linewidth=0.6)
    ax.axvline(1.0, color="#333333", linewidth=0.9, linestyle="--")
    for i, r in enumerate(d.itertuples()):
        ax.text(r.odds_ratio + 0.02, i,
                f"{int(r.n_causal)}/{int(r.n_source)}  (OR={r.odds_ratio:.2f})",
                va="center", fontsize=7.4, color="#333333")
    ax.set_yticks(y); ax.set_yticklabels(d["source_cell"], fontsize=8.4)
    ax.set_xlabel("enrichment of causal hits (odds ratio vs. all instrumented)")
    ax.set_title("Which immune cell of origin yields more causal proteins", loc="left")
    ax.set_xlim(0, max(1.6, d["odds_ratio"].max() * 1.25))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "X4_cellsource_enrichment.png"))
    plt.close(fig)
    print("  X4_cellsource_enrichment.png")


def main():
    print("[38] Extra paper-style panels ...")
    mr  = pd.read_csv(os.path.join(GC, "cis_MR_phenome_results.tsv"), sep="\t")
    nov = pd.read_csv(os.path.join(GC, "novelty_engine_ranked.tsv"), sep="\t")
    cs  = pd.read_csv(os.path.join(GC, "phenome_cellsource_map.tsv"), sep="\t")
    panel_x1_diseases_per_protein(mr)
    panel_x2_shared_categories(mr)
    panel_x3_score_decomposition(nov)
    panel_x4_cellsource(cs)
    print(f"[38] Done -> {OUT}")


if __name__ == "__main__":
    main()
