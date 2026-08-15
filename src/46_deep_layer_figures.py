#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
46 - DEEP-LAYER figures for the whole-phenome extension
=======================================================
Visualises the three real deep-analysis layers built on top of the pan-phenome
cis-MR (src/41): additive, nothing existing is modified.

  DL1  Phenome-wide colocalization (src/43, coloc_ALL_finngen_results.tsv):
       PP.H4 distribution across all pan-phenome hits + top colocalised pairs.
  DL2  Two-population replication (src/44, uk_panphenome_concordance.tsv):
       FinnGen(Finland) vs UKB(England) effect concordance scatter + summary.
  DL2b Widened two-population replication (src/52, uk_panphenome_concordance_ALL.tsv):
       same design plus normalised name matches, UK Biobank-only cohorts.
  DL3  Extended blood immune-trait MR (src/45, extended_cell_crp_MR_results.tsv):
       causal immune-protein count per blood trait (counts / fractions / CRP).

Outputs -> 08_figures/paper_style/DL1_*, DL2_*, DL3_*
Run:  python src/46_deep_layer_figures.py
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

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
    "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.labelsize": 9.5, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "legend.frameon": False,
})


# --------------------------------------------------------------------------- #
def panel_dl1_coloc():
    f = os.path.join(GC, "coloc_ALL_finngen_results.tsv")
    if not os.path.exists(f):
        print("DL1 skipped (no coloc file)"); return
    d = pd.read_csv(f, sep="\t")
    d = d[d["nsnp"] >= 5].copy()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4),
                                   gridspec_kw={"width_ratios": [1, 1.25]})
    # left: PP.H4 histogram
    ax1.hist(d["PP_H4"].dropna(), bins=30, color="#55A868", edgecolor="white")
    ax1.axvline(0.8, color="#C44E52", ls="--", lw=1.0)
    n8 = int((d["PP_H4"] >= 0.8).sum()); n5 = int((d["PP_H4"] >= 0.5).sum())
    ax1.set_xlabel("PP.H4 (shared causal variant)")
    ax1.set_ylabel("gene-disease pairs")
    ax1.set_title(f"DL1  Phenome-wide colocalization\n{len(d)} tested pairs | "
                  f"PP.H4>=0.8: {n8} | >=0.5: {n5}")
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)
    # right: top colocalised pairs
    top = d.sort_values("PP_H4", ascending=False).head(22).iloc[::-1]
    lab = (top["gene"] + "  ->  " + top["disease"].str.slice(0, 40))
    y = np.arange(len(top))
    ax2.hlines(y, 0, top["PP_H4"], color="#bbb", lw=1.2, zorder=1)
    ax2.scatter(top["PP_H4"], y, c="#55A868", s=40, zorder=2,
                edgecolors="white", linewidths=0.6)
    ax2.set_yticks(y); ax2.set_yticklabels(lab, fontsize=6.5)
    ax2.set_xlim(0, 1.02); ax2.axvline(0.8, color="#C44E52", ls="--", lw=0.9)
    ax2.set_xlabel("PP.H4")
    ax2.set_title("Top colocalised immune-protein -> disease pairs")
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "DL1_phenome_coloc.png"))
    plt.close(fig)
    print("wrote DL1")


# --------------------------------------------------------------------------- #
def panel_dl2_replication():
    f = os.path.join(GC, "uk_panphenome_concordance.tsv")
    if not os.path.exists(f):
        print("DL2 skipped (no UKB file)"); return
    d = pd.read_csv(f, sep="\t")
    d = d[(d["fin_OR"] > 0) & (d["uk_OR"] > 0)].copy()
    lx = np.log(d["fin_OR"]); ly = np.log(d["uk_OR"])
    val = d["two_population_validated"].astype(bool)
    conc = d["concordant"].astype(bool)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4),
                                   gridspec_kw={"width_ratios": [1.3, 1]})
    # left: log-OR concordance scatter
    ax1.axhline(0, color="#999", lw=0.7); ax1.axvline(0, color="#999", lw=0.7)
    lim = max(abs(lx).max(), abs(ly).max()) * 1.1
    ax1.plot([-lim, lim], [-lim, lim], color="#ccc", ls=":", lw=0.9)
    ax1.scatter(lx[~conc], ly[~conc], s=26, c="#bbb", edgecolors="none",
                alpha=0.7, label="discordant")
    ax1.scatter(lx[conc & ~val], ly[conc & ~val], s=30, c="#4C72B0",
                edgecolors="none", alpha=0.8, label="concordant")
    ax1.scatter(lx[val], ly[val], s=42, c="#C44E52", edgecolors="white",
                linewidths=0.5, label="two-population validated")
    ax1.set_xlim(-lim, lim); ax1.set_ylim(-lim, lim)
    ax1.set_xlabel("FinnGen (Finland)  ln(OR)")
    ax1.set_ylabel("UK Biobank (England)  ln(OR)")
    ax1.set_title("DL2  FinnGen vs UK Biobank causal-effect replication")
    ax1.legend(loc="upper left", fontsize=7.5)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)
    # right: summary bars
    n_pair = len(d); n_conc = int(conc.sum()); n_val = int(val.sum())
    cats = ["tested\npairs", "same\ndirection", "two-pop.\nvalidated"]
    vals = [n_pair, n_conc, n_val]
    cols = ["#8C8C8C", "#4C72B0", "#C44E52"]
    ax2.bar(np.arange(3), vals, color=cols, edgecolor="white")
    for i, v in enumerate(vals):
        ax2.text(i, v + max(vals) * 0.01, str(v), ha="center", va="bottom", fontsize=9)
    ax2.set_xticks(np.arange(3)); ax2.set_xticklabels(cats, fontsize=8)
    ax2.set_ylabel("gene-disease pairs")
    ax2.set_title(f"{d['phenocode'].nunique()} exact-code diseases | "
                  f"{d['gene_symbol'].nunique()} genes")
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "DL2_uk_replication.png"))
    plt.close(fig)
    print("wrote DL2")


# --------------------------------------------------------------------------- #
def panel_dl3_bloodtraits():
    f = os.path.join(GC, "extended_cell_crp_MR_results.tsv")
    if not os.path.exists(f):
        print("DL3 skipped (no extended cell/CRP file)"); return
    d = pd.read_csv(f, sep="\t")
    sig = d[d["trait_FDR"] < 0.05]
    per = (sig.groupby(["trait", "trait_group"]).size()
           .reset_index(name="n").sort_values("n"))
    gcol = {"count": "#4C72B0", "fraction": "#DD8452", "inflammation": "#C44E52"}
    colors = [gcol.get(g, "#8C8C8C") for g in per["trait_group"]]
    fig, ax = plt.subplots(figsize=(9, 6.5))
    y = np.arange(len(per))
    ax.barh(y, per["n"].values, color=colors, edgecolor="white")
    ax.set_yticks(y); ax.set_yticklabels(per["trait"].values, fontsize=8)
    for i, v in enumerate(per["n"].values):
        ax.text(v + 0.3, i, str(int(v)), va="center", fontsize=7.5)
    ax.set_xlabel("causal immune proteins (FDR<0.05)")
    ax.set_title("DL3  Immune proteome -> blood immune traits "
                 f"({d['trait'].nunique()} traits, {len(d)} tests)")
    leg = [Line2D([0], [0], marker="s", color="w", markerfacecolor=gcol["count"],
                  markersize=9, label="cell count"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor=gcol["fraction"],
                  markersize=9, label="differential fraction"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor=gcol["inflammation"],
                  markersize=9, label="inflammation (CRP)")]
    ax.legend(handles=leg, loc="lower right", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "DL3_blood_traits.png"))
    plt.close(fig)
    print("wrote DL3")


def panel_dl2b_replication_expanded():
    """Widened two-population arm (src/52): exact-code + normalised name matches."""
    f = os.path.join(GC, "uk_panphenome_concordance_ALL.tsv")
    if not os.path.exists(f):
        print("DL2b skipped (no expanded UKB file)"); return
    d = pd.read_csv(f, sep="\t")
    d = d[(d["fin_OR"] > 0) & (d["uk_OR"] > 0)].copy()
    lx = np.log(d["fin_OR"]); ly = np.log(d["uk_OR"])
    val = d["two_population_validated"].astype(bool)
    conc = d["concordant"].astype(bool)
    nm = d["match_method"].eq("name-match")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4),
                                   gridspec_kw={"width_ratios": [1.3, 1]})
    # FinnGen log-odds and UKB (largely linear-model) betas are on different
    # scales, so each axis is scaled independently; what is being read off this
    # panel is the SIGN agreement (quadrant), not the slope.
    ax1.axhline(0, color="#999", lw=0.7); ax1.axvline(0, color="#999", lw=0.7)
    xl = abs(lx).max() * 1.12; yl = abs(ly).max() * 1.12
    ax1.axhspan(0, yl, xmin=0.5, color="#eef4ea", zorder=0)
    ax1.axhspan(-yl, 0, xmax=0.5, color="#eef4ea", zorder=0)
    ax1.scatter(lx[~conc], ly[~conc], s=26, c="#bbb", edgecolors="none",
                alpha=0.7, label="discordant")
    ax1.scatter(lx[conc & ~val], ly[conc & ~val], s=30, c="#4C72B0",
                edgecolors="none", alpha=0.8, label="concordant")
    ax1.scatter(lx[val], ly[val], s=42, c="#C44E52", edgecolors="white",
                linewidths=0.5, label="two-population validated")
    ax1.scatter(lx[nm], ly[nm], s=95, facecolors="none", edgecolors="#1f6f3f",
                linewidths=0.9, label="newly reachable (name-matched)")
    ax1.set_xlim(-xl, xl); ax1.set_ylim(-yl, yl)
    ax1.set_xlabel("FinnGen (Finland)  ln(OR)")
    ax1.set_ylabel("UK Biobank (England)  ln(OR)   [note: axis scaled independently]")
    ax1.set_title("DL2b  Widened FinnGen vs UK Biobank replication\n"
                  "shaded quadrants = same causal direction", fontsize=10)
    ax1.legend(loc="upper left", fontsize=7.5)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)
    # right: exact-code vs name-match contribution, stacked
    stages = ["tested\npairs", "same\ndirection", "two-pop.\nvalidated"]
    ex = [int((~nm).sum()), int((conc & ~nm).sum()), int((val & ~nm).sum())]
    na = [int(nm.sum()), int((conc & nm).sum()), int((val & nm).sum())]
    x = np.arange(3)
    ax2.bar(x, ex, color="#4C72B0", edgecolor="white", label="exact-code (src/44)")
    ax2.bar(x, na, bottom=ex, color="#1f6f3f", edgecolor="white",
            label="name-matched (new)")
    for i, (a, b) in enumerate(zip(ex, na)):
        ax2.text(i, a + b + max(np.add(ex, na)) * 0.012, str(a + b),
                 ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax2.set_xticks(x); ax2.set_xticklabels(stages, fontsize=8)
    ax2.set_ylabel("gene-disease pairs")
    ax2.set_title(f"{d['phenocode'].nunique()} diseases | "
                  f"{d['gene_symbol'].nunique()} genes | UK Biobank-only cohorts")
    ax2.legend(fontsize=7.5, frameon=False)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "DL2b_uk_replication_expanded.png"))
    plt.close(fig)
    print("wrote DL2b")


def main():
    panel_dl1_coloc()
    panel_dl2_replication()
    panel_dl2b_replication_expanded()
    panel_dl3_bloodtraits()
    print("[46] done ->", OUT)


if __name__ == "__main__":
    main()
