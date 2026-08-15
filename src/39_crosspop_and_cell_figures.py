#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
39 - Cross-population replication + immune-cell-count figures
=============================================================
Two data-gated reference-paper analogues, now buildable because the throttled
OpenGWAS jobs have finished:

  R1  Finland (FinnGen) vs England (UK Biobank) two-population replication
      (their Fig 3i,j replication analogue) -- from finngen_vs_uk_concordance.tsv
  R2  immune protein -> circulating immune-CELL-COUNT MR
      (their Fig 3 metabolite->trait arm analogue) -- from
      immune_cell_count_MR_results.tsv

Purely additive: reads existing tables, writes new PNGs into
08_figures/paper_style/, modifies nothing.

Run:
    python src/39_crosspop_and_cell_figures.py
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

CELL_COLORS = {
    "Lymphocyte count": "#4C72B0",
    "Monocyte count":   "#DD8452",
    "Neutrophil count": "#C44E52",
    "Eosinophil count": "#55A868",
    "Basophil count":   "#8172B3",
    "Leukocyte count":  "#937860",
}


# --------------------------------------------------------------------------- #
def panel_r1_crosspop(c):
    """R1 (their Fig 3i,j): Finland vs England effect concordance scatter."""
    d = c.copy()
    d["fin_beta"] = np.log(d["fin_OR"])
    d["uk_beta"]  = np.log(d["uk_OR"])
    d = d[np.isfinite(d["fin_beta"]) & np.isfinite(d["uk_beta"])]

    n = len(d)
    n_conc = int(d["concordant"].sum())
    n_val  = int(d["two_population_validated"].sum())

    fig, ax = plt.subplots(figsize=(6.8, 6.4))
    # independent per-axis scales (FinnGen ln(OR) and Neale-UKB betas differ in
    # magnitude; only the SIGN defines direction concordance)
    xlim = d["fin_beta"].abs().max() * 1.18
    ylim = d["uk_beta"].abs().max() * 1.25
    # shade the two concordant (same-sign) quadrants
    ax.axhspan(0, ylim, xmin=0.5, xmax=1.0, color="#2A9D8F", alpha=0.05, zorder=0)
    ax.axhspan(-ylim, 0, xmin=0.0, xmax=0.5, color="#2A9D8F", alpha=0.05, zorder=0)
    ax.axhline(0, color="#999999", linewidth=0.8)
    ax.axvline(0, color="#999999", linewidth=0.8)

    disc = ~d["concordant"]
    val  = d["two_population_validated"]
    ax.scatter(d.loc[disc, "fin_beta"], d.loc[disc, "uk_beta"], s=26,
               facecolor="#DDDDDD", edgecolor="#999999", linewidth=0.5,
               label=f"discordant ({int(disc.sum())})", zorder=2)
    conc_only = d["concordant"] & ~val
    ax.scatter(d.loc[conc_only, "fin_beta"], d.loc[conc_only, "uk_beta"], s=32,
               facecolor="#7FA8C9", edgecolor="white", linewidth=0.5,
               label=f"concordant ({int(conc_only.sum())})", zorder=3)
    ax.scatter(d.loc[val, "fin_beta"], d.loc[val, "uk_beta"], s=52,
               facecolor="#C44E52", edgecolor="white", linewidth=0.6,
               label=f"two-population validated ({n_val})", zorder=4)

    # label a few strongest validated pairs
    for _, r in d[val].reindex(d[val]["fin_beta"].abs().sort_values(ascending=False).index).head(8).iterrows():
        ax.annotate(f"{r['gene_symbol']}", (r["fin_beta"], r["uk_beta"]),
                    fontsize=6.6, ha="center", va="bottom", color="#222222")

    ax.set_xlim(-xlim, xlim); ax.set_ylim(-ylim, ylim)
    ax.set_xlabel("FinnGen effect  ln(OR)   \u2190 protective | risk \u2192")
    ax.set_ylabel("UK Biobank effect  (SD units)   \u2190 protective | risk \u2192")
    ax.set_title("Two-population replication of immune\u2192disease causal effects\n"
                 "Finland (FinnGen) vs England (UK Biobank)", loc="left", fontsize=10.3)
    ax.text(0.02, 0.02, f"{n_conc}/{n} directionally concordant "
                        f"({n_conc/n*100:.0f}%)\n{n_val} validated at UK FDR<0.05",
            transform=ax.transAxes, fontsize=8, color="#555555", va="bottom")
    ax.legend(fontsize=7.6, loc="upper left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "R1_crosspop_replication.png"))
    plt.close(fig)
    print("  R1_crosspop_replication.png")


# --------------------------------------------------------------------------- #
def panel_r2_cellcount_manhattan(m):
    """R2 (their Fig 3 trait-arm): immune protein -> immune cell count MR,
    mirrored by cell type."""
    d = m.copy()
    d = d[np.isfinite(d["cell_p"]) & (d["cell_p"] > 0)]
    d["signed_logp"] = -np.log10(d["cell_p"]) * np.where(d["cell_beta"] >= 0, 1, -1)

    traits = [t for t in CELL_COLORS if t in set(d["cell_trait"])]
    order, xpos, boundaries = {}, {}, []
    x = 0
    rng = np.random.default_rng(11)
    for t in traits:
        genes = sorted(d.loc[d["cell_trait"] == t, "gene_symbol"].unique())
        a = x
        for _ in genes:
            x += 1
        boundaries.append((t, a, x - 1))
    # x = jittered position within each trait band
    xvals = []
    tstart = {t: b[1] for t, b in zip(traits, boundaries)}
    tspan  = {t: (b[2] - b[1] + 1) for t, b in zip(traits, boundaries)}
    for _, r in d.iterrows():
        t = r["cell_trait"]
        xvals.append(tstart[t] + rng.uniform(0, tspan[t] - 1))
    d = d.assign(x=xvals)

    thr = -np.log10(FDR_SIG)
    fig, ax = plt.subplots(figsize=(13, 6))
    for t in traits:
        sub = d[d["cell_trait"] == t]
        sig = sub["cell_FDR"] < FDR_SIG
        ax.scatter(sub.loc[~sig, "x"], sub.loc[~sig, "signed_logp"], s=6,
                   color=CELL_COLORS[t], alpha=0.18, linewidths=0, rasterized=True)
        ax.scatter(sub.loc[sig, "x"], sub.loc[sig, "signed_logp"], s=20,
                   color=CELL_COLORS[t], alpha=0.9, linewidths=0.3, edgecolor="white")
    ax.axhline(0, color="#333333", linewidth=0.9)
    ax.axhline(thr, color="grey", linewidth=0.8, linestyle="--")
    ax.axhline(-thr, color="grey", linewidth=0.8, linestyle="--")

    # label strongest hits
    hits = d[d["cell_FDR"] < FDR_SIG].copy()
    hits["abs"] = hits["signed_logp"].abs()
    for _, r in hits.sort_values("abs", ascending=False).head(20).iterrows():
        ax.annotate(r["gene_symbol"], (r["x"], r["signed_logp"]), fontsize=6.4,
                    ha="center", va="bottom" if r["signed_logp"] > 0 else "top",
                    color="#222222")
    ymin, ymax = ax.get_ylim()
    band = ymax + (ymax - ymin) * 0.14
    ax.set_ylim(ymin, band)
    lab_y = ymax + (ymax - ymin) * 0.04
    for t, a, b in boundaries:
        ax.axvspan(a - 0.5, b + 0.5, color=CELL_COLORS[t], alpha=0.05, zorder=0)
        ax.text((a + b) / 2, lab_y, t.replace(" count", ""), ha="center", va="bottom",
                fontsize=9.5, fontweight="bold", color=CELL_COLORS[t])
    ax.set_xticks([])
    ax.set_xlim(-1, x)
    ax.set_ylabel("signed \u2212log$_{10}$($P_{MR}$)   \u2190 decreases | increases \u2192")
    ax.set_title("Plasma-immune protein \u2192 circulating immune-cell-count Mendelian randomization",
                 loc="left", pad=14)
    n_hit = int((d["cell_FDR"] < FDR_SIG).sum())
    ax.text(0.005, 0.03, f"{len(d):,} cis-MR tests \u00b7 {n_hit} at FDR<0.05 "
                         f"(Neale/UKB blood counts)", transform=ax.transAxes,
            fontsize=8, color="#555555")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "R2_cellcount_manhattan.png"))
    plt.close(fig)
    print("  R2_cellcount_manhattan.png")


# --------------------------------------------------------------------------- #
def panel_r3_cellcount_bars(m):
    """R3: causal immune proteins per immune cell count, split by direction."""
    hits = m[m["cell_FDR"] < FDR_SIG].copy()
    traits = [t for t in CELL_COLORS if t in set(hits["cell_trait"])]
    up = [int(((hits["cell_trait"] == t) & (hits["cell_beta"] > 0)).sum()) for t in traits]
    dn = [int(((hits["cell_trait"] == t) & (hits["cell_beta"] < 0)).sum()) for t in traits]

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    y = np.arange(len(traits))
    ax.barh(y, up, color="#C44E52", edgecolor="white", linewidth=0.5, label="increases count")
    ax.barh(y, [-v for v in dn], color="#4C72B0", edgecolor="white", linewidth=0.5,
            label="decreases count")
    ax.axvline(0, color="#333333", linewidth=0.9)
    for i, (u, d_) in enumerate(zip(up, dn)):
        ax.text(u + 1.2, i, str(u), va="center", fontsize=8, color="#C44E52")
        ax.text(-d_ - 1.2, i, str(d_), va="center", ha="right", fontsize=8, color="#4C72B0")
    ax.set_yticks(y); ax.set_yticklabels([t.replace(" count", "") for t in traits], fontsize=9)
    ax.set_xlabel("decreases  \u2190   causal immune proteins   \u2192  increases")
    ax.set_title("Causal immune proteins per circulating cell count, by direction", loc="left")
    ax.legend(fontsize=8, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "R3_cellcount_direction.png"))
    plt.close(fig)
    print("  R3_cellcount_direction.png")


def main():
    print("[39] Cross-population + cell-count figures ...")
    c = pd.read_csv(os.path.join(GC, "finngen_vs_uk_concordance.tsv"), sep="\t")
    m = pd.read_csv(os.path.join(GC, "immune_cell_count_MR_results.tsv"), sep="\t")
    print(f"    concordance rows: {len(c)}  |  cell-count tests: {len(m)}")
    panel_r1_crosspop(c)
    panel_r2_cellcount_manhattan(m)
    panel_r3_cellcount_bars(m)
    print(f"[39] Done -> {OUT}")


if __name__ == "__main__":
    main()
