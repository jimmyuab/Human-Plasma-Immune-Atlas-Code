#!/usr/bin/env python
"""HDDM Layer 54 - Step 8 : figures for the cis-MR results."""
import os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="whitegrid", context="talk")

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality")
FIG  = os.path.join(ROOT, "08_figures", "main_figures")
res  = pd.read_csv(os.path.join(GEN, "cis_MR_immune_results.tsv"), sep="\t")
res["nlog10p"] = -np.log10(res["MR_p"].clip(lower=1e-300))

def save(fig, n):
    fig.savefig(os.path.join(FIG, n), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig); print(" wrote", n)

# Fig13 volcano
def volcano():
    fig, ax = plt.subplots(figsize=(11, 8))
    sig = res.FDR < 0.05
    ax.scatter(res.loc[~sig, "MR_beta"], res.loc[~sig, "nlog10p"], s=14, c="#cccccc", label="ns")
    up = sig & (res.MR_beta > 0); dn = sig & (res.MR_beta < 0)
    ax.scatter(res.loc[up, "MR_beta"], res.loc[up, "nlog10p"], s=40, c="#c0392b", label="FDR<0.05 risk")
    ax.scatter(res.loc[dn, "MR_beta"], res.loc[dn, "nlog10p"], s=40, c="#2471a3", label="FDR<0.05 protective")
    for _, r in res[sig].head(14).iterrows():
        ax.annotate(f"{r.gene_symbol}\u2192{r.disease.split()[0]}",
                    (r.MR_beta, r.nlog10p), fontsize=8,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=1)
    ax.set_xlabel("cis-MR effect (log-OR per SD cis-expression)")
    ax.set_ylabel(r"$-\log_{10}$ MR p-value")
    ax.set_title("Fig.13  cis-MR: immune-gene expression \u2192 immune-disease risk", fontweight="bold")
    ax.legend(fontsize=11)
    save(fig, "Fig13_MR_volcano.png")

# Fig14 forest of top 18 by FDR
def forest():
    top = res.sort_values("FDR").head(18).iloc[::-1]
    fig, ax = plt.subplots(figsize=(10, 9))
    y = np.arange(len(top))
    ax.errorbar(top.OR, y, xerr=[top.OR - top.OR_l95, top.OR_u95 - top.OR],
                fmt="o", color="#1f496e", ecolor="#888", capsize=3, ms=7)
    ax.axvline(1, ls="--", c="red", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.gene_symbol} \u2192 {r.disease}" for _, r in top.iterrows()], fontsize=10)
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio per SD cis-expression (log scale)")
    ax.set_title("Fig.14  Top causal immune gene\u2013disease pairs (FDR<0.05)", fontweight="bold")
    save(fig, "Fig14_MR_forest.png")

# Fig15 heatmap MR-Z for genes significant in >=1 disease
def heatmap():
    siggenes = res.loc[res.FDR < 0.05, "gene_symbol"].unique()
    sub = res[res.gene_symbol.isin(siggenes)]
    piv = sub.pivot_table(index="gene_symbol", columns="disease", values="MR_z", aggfunc="first")
    piv = piv.loc[piv.abs().max(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(13, max(6, .45 * len(piv))))
    sns.heatmap(piv, cmap="RdBu_r", center=0, vmin=-8, vmax=8, ax=ax,
                cbar_kws={"label": "MR Z (risk + / protective \u2212)"}, linewidths=.4)
    ax.set_title("Fig.15  Causal immune genes \u00d7 disease (MR Z-score)", fontweight="bold")
    ax.set_xlabel(""); ax.set_ylabel("immune gene")
    plt.xticks(rotation=40, ha="right")
    save(fig, "Fig15_MR_heatmap.png")

# Fig16 significant counts per disease
def perdisease():
    c = res[res.FDR < 0.05].groupby("disease").size().sort_values()
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(c.index, c.values, color=sns.color_palette("rocket", len(c)))
    for i, v in enumerate(c.values):
        ax.text(v + .05, i, str(v), va="center", fontweight="bold")
    ax.set_xlabel("number of causal immune genes (FDR<0.05)")
    ax.set_title("Fig.16  Causal immune genes per disease", fontweight="bold")
    save(fig, "Fig16_MR_per_disease.png")

if __name__ == "__main__":
    print("MR figures ->", FIG)
    volcano(); forest(); heatmap(); perdisease()
    print("done")
