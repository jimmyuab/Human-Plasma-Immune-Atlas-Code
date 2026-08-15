#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
40 - Complete the reference-atlas figure TYPESET (remaining distinctive panels)
===============================================================================
Fills in every remaining reference-paper (You et al., Nat Metab 2025) panel TYPE
that our plasma-immune data can honestly support, so the gallery covers the full
~20-panel typeset rather than a subset.

New panel types added here (with their reference analogue):
  C1  circular per-disease contribution plot  (their Fig 6g,h) x2 diseases
  V1  cross-category overlap / set map         (their Fig 2f Venn)
  E1  evidence-source stacked bar              (their Fig 2e source split)
  S1  eQTL-MR vs pQTL-MR effect scatter        (their Fig 2g,h effect comparison)
  T1  evidence-tier ladder                     (their Fig 6 validation summary)
  D1  per-disease directional volcano grid     (their Fig 2 multi-disease)

Purely additive: reads existing tables, writes new PNGs into
08_figures/paper_style/, modifies nothing.

Run:
    python src/40_full_paper_typeset.py
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
    "font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333", "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.labelsize": 9.5, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "legend.frameon": False,
})

CAT_COLORS = {
    "Autoimmune": "#D1495B", "Cardiovascular": "#2E6F95", "Metabolic": "#E9963A",
    "Neuro/Aging": "#6A4C93", "Renal": "#2A9D8F",
}
CLASS_COLORS = {
    "Immune-cell-enriched": "#4C72B0", "Cytokine/Interleukin": "#DD8452",
    "CD/Leukocyte surface": "#55A868", "TNF superfamily": "#C44E52",
    "Complement/Coagulation": "#8172B3", "Immune checkpoint": "#937860",
    "Chemokine axis": "#DA8BC3", "Ig/B-cell/Fc": "#8C8C8C",
    "Interferon axis": "#CCB974", "HLA/Antigen presentation": "#64B5CD",
    "Acute-phase": "#B07AA1",
}
def cls_col(c): return CLASS_COLORS.get(c, "#BBBBBB")


# --------------------------------------------------------------------------- #
def panel_c1_circular(mr, disease, fname):
    """C1 (their Fig 6g,h): circular contribution plot of every causal immune
    protein for one disease; bar length = -log10(P), direction = risk/protective,
    colour = immune class."""
    d = mr[(mr["disease"] == disease) & (mr["FDR"] < FDR_SIG)].copy()
    if not len(d):
        print(f"  (no hits for {disease} -> {fname} skipped)"); return
    d["logp_raw"] = -np.log10(np.clip(d["MR_p"], 1e-300, None))
    # cap extreme -log10(P) so the wheel stays legible (a handful of hits reach
    # P~0); capped bars are flagged with a '+' so nothing is hidden
    cap = float(np.percentile(d["logp_raw"], 80))
    cap = max(cap, 8.0)
    d["logp"] = np.clip(d["logp_raw"], None, cap)
    d["capped"] = d["logp_raw"] > cap
    d = d.sort_values("immune_class")
    n = len(d)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    width = 2 * np.pi / n * 0.9
    rmax = d["logp"].max()

    fig = plt.figure(figsize=(8.4, 8.4))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(np.pi / 2); ax.set_theta_direction(-1)
    for a, (_, r) in zip(ang, d.iterrows()):
        h = r["logp"]
        col = cls_col(r["immune_class"])
        risk = r["OR"] >= 1
        ax.bar(a, h, width=width, bottom=2, color=col,
               edgecolor="white", linewidth=0.4,
               alpha=0.95 if risk else 0.55, hatch=None if risk else "///")
        # label gene
        rot = np.degrees(a)
        ha = "left"
        if 90 < rot < 270:
            rot += 180; ha = "right"
        lab = r["gene_symbol"] + ("\u207a" if r["capped"] else "")
        ax.text(a, 2 + h + rmax * 0.06, lab, rotation=rot,
                rotation_mode="anchor", ha=ha, va="center", fontsize=6.0,
                color="#222222")
    ax.set_ylim(0, 2 + rmax * 1.25)
    ax.set_yticklabels([]); ax.set_xticks([]); ax.spines["polar"].set_visible(False)
    ax.text(0, 0, disease.replace(" ", "\n"), ha="center", va="center",
            fontsize=11, fontweight="bold")
    # legends
    classes = [c for c in CLASS_COLORS if c in set(d["immune_class"])]
    h1 = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                 markerfacecolor=cls_col(c), markeredgecolor="white", label=c)
          for c in classes]
    leg1 = ax.legend(handles=h1, fontsize=7, loc="upper left",
                     bbox_to_anchor=(-0.14, 1.12), title="immune class",
                     title_fontsize=7.6)
    ax.add_artist(leg1)
    h2 = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                 markerfacecolor="#888", markeredgecolor="white", label="risk (OR>1)"),
          Line2D([0], [0], marker="s", linestyle="", markersize=8, alpha=0.55,
                 markerfacecolor="#888", markeredgecolor="white", label="protective (OR<1)")]
    ax.legend(handles=h2, fontsize=7.2, loc="upper right", bbox_to_anchor=(1.14, 1.12))
    ax.set_title(f"Causal immune-protein contribution wheel \u2014 {disease}\n"
                 f"bar = \u2212log10(P$_{{MR}}$) capped at {cap:.0f} (\u207a=exceeds), "
                 f"colour = immune class ({n} proteins)",
                 loc="center", pad=26, fontsize=10.5)
    fig.savefig(os.path.join(OUT, fname))
    plt.close(fig)
    print(f"  {fname}")


# --------------------------------------------------------------------------- #
def panel_v1_overlap(mr):
    """V1 (their Fig 2f): overlap map of causal immune proteins shared between
    the largest disease categories (upset-style presence matrix)."""
    hits = mr[mr["FDR"] < FDR_SIG].copy()
    cats = ["Autoimmune", "Cardiovascular", "Metabolic", "Neuro/Aging", "Renal"]
    cats = [c for c in cats if c in set(hits["disease_category"])]
    gene_cats = hits.groupby("gene_symbol")["disease_category"].apply(lambda s: set(s))
    shared = gene_cats[gene_cats.apply(len) >= 2]
    if not len(shared):
        print("  (no cross-category genes -> V1 skipped)"); return
    # build presence matrix for shared genes
    genes = sorted(shared.index, key=lambda g: (-len(shared[g]), g))
    mat = np.zeros((len(cats), len(genes)))
    for j, g in enumerate(genes):
        for i, c in enumerate(cats):
            if c in shared[g]:
                mat[i, j] = 1

    fig, ax = plt.subplots(figsize=(max(7, 0.32 * len(genes) + 2), 3.6))
    for i, c in enumerate(cats):
        for j, g in enumerate(genes):
            if mat[i, j]:
                ax.scatter(j, i, s=90, color=CAT_COLORS[c], edgecolor="white", linewidth=0.6, zorder=3)
        # connect dots within a gene column
    for j, g in enumerate(genes):
        ys = [i for i in range(len(cats)) if mat[i, j]]
        if len(ys) >= 2:
            ax.plot([j, j], [min(ys), max(ys)], color="#999999", linewidth=1.2, zorder=2)
    ax.set_yticks(range(len(cats)))
    ax.set_yticklabels(cats, fontsize=8.6)
    for tick, c in zip(ax.get_yticklabels(), cats):
        tick.set_color(CAT_COLORS[c]); tick.set_fontweight("bold")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_xlim(-1, len(genes)); ax.set_ylim(-0.6, len(cats) - 0.4)
    ax.set_title("Immune proteins causal across \u22652 disease categories "
                 "(cross-chapter set map)", loc="left")
    ax.grid(True, axis="y", color="#EEEEEE", linewidth=0.5); ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "V1_crosscategory_setmap.png"))
    plt.close(fig)
    print("  V1_crosscategory_setmap.png")


# --------------------------------------------------------------------------- #
def panel_e1_evidence_source(mr, nov, ft):
    """E1 (their Fig 2e): per-category hit counts split by strongest evidence
    SOURCE reached (MR only / +coloc / +pQTL protein-level)."""
    hits = mr[mr["FDR"] < FDR_SIG][["gene_symbol", "disease", "disease_category"]].copy()
    hits["key"] = hits["gene_symbol"] + "|" + hits["disease"]
    # coloc level from novelty (PP_H4) and pQTL from tiers
    nv = nov.set_index(nov["gene_symbol"] + "|" + nov["disease"])
    ftk = ft.set_index(ft["gene_symbol"] + "|" + ft["disease"])
    def level(k):
        pph4 = nv["PP_H4"].get(k, np.nan) if k in nv.index else np.nan
        pqtl = ftk["pQTL_OR"].get(k, np.nan) if k in ftk.index else np.nan
        if isinstance(pqtl, pd.Series): pqtl = pqtl.iloc[0]
        if isinstance(pph4, pd.Series): pph4 = pph4.iloc[0]
        if pd.notna(pqtl):
            return "MR + coloc + pQTL (protein-level)"
        if pd.notna(pph4) and pph4 >= 0.8:
            return "MR + colocalization"
        return "MR only"
    hits["level"] = hits["key"].map(level)
    levels = ["MR only", "MR + colocalization", "MR + coloc + pQTL (protein-level)"]
    lvl_col = {"MR only": "#BFD7EA",
               "MR + colocalization": "#7FA8C9",
               "MR + coloc + pQTL (protein-level)": "#1B4965"}
    cats = ["Autoimmune", "Cardiovascular", "Metabolic", "Neuro/Aging", "Renal"]
    cats = [c for c in cats if c in set(hits["disease_category"])]
    tab = (hits.groupby(["disease_category", "level"]).size()
                .unstack(fill_value=0).reindex(index=cats, columns=levels, fill_value=0))

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    bottom = np.zeros(len(cats))
    for lv in levels:
        vals = tab[lv].values
        ax.bar(range(len(cats)), vals, bottom=bottom, width=0.66,
               color=lvl_col[lv], edgecolor="white", linewidth=0.5, label=lv)
        bottom += vals
    for i, t in enumerate(bottom):
        ax.text(i, t + 0.5, int(t), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(range(len(cats))); ax.set_xticklabels(cats, fontsize=9)
    for tick, c in zip(ax.get_xticklabels(), cats):
        tick.set_color(CAT_COLORS[c]); tick.set_fontweight("bold")
    ax.set_ylabel("causal gene\u2013disease hits (FDR<0.05)")
    ax.set_title("Depth of evidence per disease category", loc="left")
    ax.legend(fontsize=7.6, loc="upper right", title="strongest evidence reached",
              title_fontsize=7.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "E1_evidence_source_stack.png"))
    plt.close(fig)
    print("  E1_evidence_source_stack.png")


# --------------------------------------------------------------------------- #
def panel_s1_pqtl_scatter(ft):
    """S1 (their Fig 2g,h): effect-comparison scatter, transcript (eQTL) MR OR
    vs protein (pQTL) MR OR, for targets with both."""
    d = ft.copy()
    d = d[d["pQTL_OR"].notna() & d["OR"].notna()]
    d["e"] = np.log(d["OR"]); d["p"] = np.log(d["pQTL_OR"])
    d = d[np.isfinite(d["e"]) & np.isfinite(d["p"])]
    n = len(d)
    conc = int((np.sign(d["e"]) == np.sign(d["p"])).sum())

    lim = max(d["e"].abs().max(), d["p"].abs().max()) * 1.2
    fig, ax = plt.subplots(figsize=(6.6, 6.4))
    ax.axhspan(0, lim, xmin=0.5, xmax=1.0, color="#2A9D8F", alpha=0.05, zorder=0)
    ax.axhspan(-lim, 0, xmin=0.0, xmax=0.5, color="#2A9D8F", alpha=0.05, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color="#CCCCCC", linewidth=0.8, linestyle="--", zorder=1)
    ax.axhline(0, color="#999999", linewidth=0.8); ax.axvline(0, color="#999999", linewidth=0.8)
    cols = [CLASS_COLORS.get(c, "#888") for c in d["immune_class"]]
    ax.scatter(d["e"], d["p"], s=70, c=cols, edgecolor="white", linewidth=0.7, zorder=3)
    for _, r in d.iterrows():
        ax.annotate(r["gene_symbol"], (r["e"], r["p"]), fontsize=6.6,
                    ha="center", va="bottom", color="#222222")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("transcript (eQTL) MR  ln(OR)")
    ax.set_ylabel("protein (pQTL) MR  ln(OR)")
    ax.set_title("Transcript vs protein causal-effect concordance", loc="left")
    ax.text(0.02, 0.02, f"{conc}/{n} same-direction ({conc/n*100:.0f}%)",
            transform=ax.transAxes, fontsize=8, color="#555555", va="bottom")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "S1_eqtl_vs_pqtl_scatter.png"))
    plt.close(fig)
    print("  S1_eqtl_vs_pqtl_scatter.png")


# --------------------------------------------------------------------------- #
def panel_t1_tier_ladder(ft):
    """T1 (their Fig 6 validation summary): evidence-tier ladder of the
    validated targets, coloured by tier, annotated with claim strength."""
    d = ft.copy().sort_values(["final_tier", "PP_H4"], ascending=[True, True])
    tier_col = {2: "#BFD7EA", 3: "#7FA8C9", 4: "#2E6F95", 5: "#1B4965"}
    labels = d["gene_symbol"] + "  \u2192  " + d["disease"]
    fig, ax = plt.subplots(figsize=(8.8, max(5, 0.34 * len(d) + 1.2)))
    y = np.arange(len(d))
    ax.barh(y, -np.log10(np.clip(d["FDR"], 1e-300, None)),
            color=[tier_col.get(int(t), "#999") for t in d["final_tier"]],
            edgecolor="white", linewidth=0.5)
    for i, r in enumerate(d.itertuples()):
        note = f"T{int(r.final_tier)}"
        if str(getattr(r, "rep_status", "")) == "replicated":
            note += " \u00b7 replicated"
        ax.text(0.15, i, note, va="center", fontsize=6.8, color="#222222")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("causal strength  \u2212log10(FDR)")
    ax.set_title("Evidence-tier ladder of validated immune targets", loc="left")
    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                      markerfacecolor=tier_col[t], markeredgecolor="white",
                      label={2: "T2 MHC-caution", 3: "T3 transcript-causal",
                             4: "T4 +colocalized", 5: "T5 protein-level"}[t])
               for t in sorted(tier_col) if t in set(d["final_tier"])]
    ax.legend(handles=handles, fontsize=7.2, loc="lower right", title="evidence tier",
              title_fontsize=7.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "T1_evidence_tier_ladder.png"))
    plt.close(fig)
    print("  T1_evidence_tier_ladder.png")


# --------------------------------------------------------------------------- #
def panel_d1_volcano_grid(mr):
    """D1 (their Fig 2 multi-disease): compact directional volcano grid for the
    six diseases with the most causal immune proteins."""
    top = mr[mr["FDR"] < FDR_SIG]["disease"].value_counts().head(6).index.tolist()
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2))
    for ax, dis in zip(axes.ravel(), top):
        d = mr[mr["disease"] == dis].copy()
        d = d[np.isfinite(d["MR_p"]) & (d["MR_p"] > 0)]
        d["logp"] = -np.log10(d["MR_p"]); d["lnor"] = np.log(d["OR"])
        sig = d["FDR"] < FDR_SIG
        cat = d["disease_category"].iloc[0]
        ax.scatter(d.loc[~sig, "lnor"], d.loc[~sig, "logp"], s=6, color="#CCCCCC",
                   alpha=0.5, linewidths=0, rasterized=True)
        ax.scatter(d.loc[sig, "lnor"], d.loc[sig, "logp"], s=18,
                   color=CAT_COLORS.get(cat, "#888"), edgecolor="white", linewidth=0.3)
        ax.axhline(-np.log10(FDR_SIG), color="grey", linewidth=0.7, linestyle="--")
        ax.axvline(0, color="#999999", linewidth=0.7)
        for _, r in d[sig].sort_values("logp", ascending=False).head(5).iterrows():
            ax.annotate(r["gene_symbol"], (r["lnor"], r["logp"]), fontsize=6,
                        ha="center", va="bottom")
        ax.set_title(f"{dis}  ({int(sig.sum())} hits)", fontsize=9, color=CAT_COLORS.get(cat, "#333"))
        ax.set_xlabel("ln(OR)"); ax.set_ylabel("\u2212log10(P)")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("Per-disease causal volcanoes \u2014 six most immune-driven diseases",
                 fontsize=12, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(OUT, "D1_volcano_grid.png"))
    plt.close(fig)
    print("  D1_volcano_grid.png")


def main():
    print("[40] Completing the reference figure typeset ...")
    mr  = pd.read_csv(os.path.join(GC, "cis_MR_phenome_results.tsv"), sep="\t")
    nov = pd.read_csv(os.path.join(GC, "novelty_engine_ranked.tsv"), sep="\t")
    ft  = pd.read_csv(os.path.join(GC, "FINAL_evidence_tiers_repl.tsv"), sep="\t")
    panel_c1_circular(mr, "Hypertension", "C1_wheel_hypertension.png")
    panel_c1_circular(mr, "Type 2 diabetes", "C1_wheel_type2diabetes.png")
    panel_v1_overlap(mr)
    panel_e1_evidence_source(mr, nov, ft)
    panel_s1_pqtl_scatter(ft)
    panel_t1_tier_ladder(ft)
    panel_d1_volcano_grid(mr)
    print(f"[40] Done -> {OUT}")


if __name__ == "__main__":
    main()
