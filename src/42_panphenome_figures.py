#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
42 - PAN-PHENOME figures (all 2,469 FinnGen R12 diseases)
=========================================================
Visualises the full immune-proteome x whole-phenome cis-MR produced by
src/41_panphenome_mr.py (cis_MR_ALL_finngen_results.tsv): 672 immune proteins
tested against every FinnGen R12 endpoint (~2,466 with data, 1.66M tests).

Panels (all new, additive; nothing existing is modified):
  PP1  phenome-wide Manhattan  -- diseases (grouped by ICD chapter) vs -log10(FDR)
  PP2  top 30 immune->disease causal hits (lollipop, coloured by direction)
  PP3  immune-causal disease count per ICD chapter (how far the immune
       proteome reaches across the medical phenome)

Outputs -> 08_figures/paper_style/PP1_*, PP2_*, PP3_*
Run:  python src/42_panphenome_figures.py
"""
import os
import re
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GC = os.path.join(ROOT, "06_genetic_causality")
OUT = os.path.join(ROOT, "08_figures", "paper_style")
os.makedirs(OUT, exist_ok=True)

FDR_SIG = 0.05

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
    "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.labelsize": 9.5, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "legend.frameon": False,
})


def short_chapter(cat):
    """Collapse the verbose FinnGen category string to a short chapter tag."""
    if not isinstance(cat, str):
        return "Other"
    m = re.match(r"\s*([IVXLC]+)\s", cat)
    roman = m.group(1) if m else ""
    # readable keyword
    kw = cat
    for junk in ["Diseases of the", "Diseases of", "from hospital discharges",
                 "from cancer register", "(ICD-O-3)"]:
        kw = kw.replace(junk, "")
    kw = re.sub(r"\([^)]*\)", "", kw)
    kw = re.sub(r"^\s*[IVXLC]+\s*", "", kw).strip(" ,")
    kw = kw[:22]
    return f"{roman} {kw}".strip()


CHAP_COLORS = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860",
    "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD", "#4E79A7", "#F28E2B",
    "#59A14F", "#E15759", "#B07AA1", "#76B7B2", "#FF9DA7", "#9C755F",
    "#BAB0AC", "#86BCB6", "#D37295", "#A0CBE8", "#FFBE7D", "#8CD17D",
]


def load():
    r = pd.read_csv(os.path.join(GC, "cis_MR_ALL_finngen_results.tsv"), sep="\t")
    r["chapter"] = r["category"].apply(short_chapter)
    r["nlfdr"] = -np.log10(np.clip(r["FDR"], 1e-320, None))
    return r


# --------------------------------------------------------------------------- #
def panel_pp1_manhattan(r):
    sig = r[r["FDR"] < FDR_SIG].copy()
    # order chapters by number of hits
    order = (sig.groupby("chapter").size().sort_values(ascending=False).index.tolist())
    order = [c for c in order if c and c != " "]
    cmap = {c: CHAP_COLORS[i % len(CHAP_COLORS)] for i, c in enumerate(order)}

    fig, ax = plt.subplots(figsize=(13, 5.2))
    x = 0
    xticks, xlabels = [], []
    rng = np.random.default_rng(7)
    for c in order:
        d = sig[sig["chapter"] == c]
        n = len(d)
        xs = x + rng.uniform(0, 1, n) * max(n, 8)
        ax.scatter(xs, d["nlfdr"], s=14, c=cmap[c], alpha=0.75,
                   edgecolors="none", rasterized=True)
        xticks.append(x + max(n, 8) / 2)
        xlabels.append(c)
        x += max(n, 8) + 4
    ax.axhline(-np.log10(FDR_SIG), color="#888", lw=0.8, ls="--")
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=55, ha="right", fontsize=6.5)
    ax.set_ylabel("-log10 FDR")
    ax.set_title("PP1  Immune-proteome causal map across the entire FinnGen R12 phenome "
                 f"({sig['phenocode'].nunique()} diseases with FDR<0.05 hits, of 2,466 tested)")
    ax.set_xlim(-4, x)
    # label a few extreme hits
    top = sig.sort_values("nlfdr", ascending=False).head(8)
    for _, row in top.iterrows():
        ax.annotate(f"{row['gene_symbol']}", (xticks[order.index(row['chapter'])], row["nlfdr"]),
                    fontsize=6.5, ha="center", va="bottom", color="#222")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "PP1_phenome_manhattan.png"))
    plt.close(fig)
    print("wrote PP1")


def panel_pp2_tophits(r):
    sig = r[r["FDR"] < FDR_SIG].copy()
    sig["lab"] = sig["gene_symbol"] + "  ->  " + sig["phenotype"].str.slice(0, 46)
    top = sig.sort_values("nlfdr", ascending=False).drop_duplicates("lab").head(30).iloc[::-1]
    col = ["#C44E52" if o > 1 else "#4C72B0" for o in top["OR"]]
    fig, ax = plt.subplots(figsize=(9, 8.5))
    y = np.arange(len(top))
    ax.hlines(y, 0, top["nlfdr"], color="#bbb", lw=1.2, zorder=1)
    ax.scatter(top["nlfdr"], y, c=col, s=42, zorder=2, edgecolors="white", linewidths=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(top["lab"], fontsize=7)
    ax.set_xlabel("-log10 FDR")
    ax.set_title("PP2  Strongest immune -> disease causal effects (all 2,466 diseases)")
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#C44E52",
                  markersize=8, label="risk-increasing (OR>1)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="#4C72B0",
                  markersize=8, label="protective (OR<1)")]
    ax.legend(handles=leg, loc="lower right", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "PP2_top_hits.png"))
    plt.close(fig)
    print("wrote PP2")


def panel_pp3_chapters(r):
    sig = r[r["FDR"] < FDR_SIG]
    per = (sig.groupby("chapter")["phenocode"].nunique()
           .sort_values(ascending=False))
    per = per[[c for c in per.index if c and c != " "]].head(18).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    colors = [CHAP_COLORS[i % len(CHAP_COLORS)] for i in range(len(per))]
    ax.barh(np.arange(len(per)), per.values, color=colors, edgecolor="white")
    ax.set_yticks(np.arange(len(per)))
    ax.set_yticklabels(per.index, fontsize=7.5)
    ax.set_xlabel("distinct diseases with >=1 immune causal protein (FDR<0.05)")
    ax.set_title("PP3  Reach of the plasma immune proteome across the medical phenome")
    for i, v in enumerate(per.values):
        ax.text(v + 0.3, i, str(int(v)), va="center", fontsize=7.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "PP3_chapter_reach.png"))
    plt.close(fig)
    print("wrote PP3")


def main():
    r = load()
    n_sig = int((r["FDR"] < FDR_SIG).sum())
    print(f"[42] {len(r):,} tests | {r['phenocode'].nunique()} diseases | "
          f"{n_sig} hits FDR<0.05 | {r[r.FDR<FDR_SIG]['phenocode'].nunique()} diseases with a hit")
    panel_pp1_manhattan(r)
    panel_pp2_tophits(r)
    panel_pp3_chapters(r)
    print("[42] done ->", OUT)


if __name__ == "__main__":
    main()
