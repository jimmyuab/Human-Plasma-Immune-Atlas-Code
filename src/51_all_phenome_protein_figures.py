#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
51 - Figures for the whole-phenome protein layer
================================================
Six panels documenting that every evidence layer now runs across the FinnGen
R12 phenome rather than a curated 28-disease list, and what the newly
phenome-wide protein (INTERVAL cis-pQTL) layer found.

  A  disease coverage per evidence layer (before/after)
  B  pan-phenome protein-level MR volcano
  C  transcript (eQTL) vs protein (pQTL) effect concordance
  D  novelty tiers across the 1,016 causal pairs, by organ system
  E  top novelty-priority targets, protein-confirmed highlighted
  F  strongest protein-level colocalizations

Outputs: 08_figures/nature/Figure17..Figure22 *.png
Run: python src/51_all_phenome_protein_figures.py
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"I:\Plasma immune atalas"
GC = os.path.join(ROOT, "06_genetic_causality")
FIG = os.path.join(ROOT, "08_figures", "nature")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

mr_all = pd.read_csv(os.path.join(GC, "cis_MR_ALL_finngen_results.tsv"), sep="\t",
                     usecols=["gene_symbol", "phenotype", "OR", "FDR", "chr"])
hits = mr_all[mr_all.FDR < 0.05]
coloc = pd.read_csv(os.path.join(GC, "coloc_ALL_finngen_results.tsv"), sep="\t")
pmr = pd.read_csv(os.path.join(GC, "pqtl_MR_ALL_finngen_results.tsv"), sep="\t")
pcol = pd.read_csv(os.path.join(GC, "pqtl_coloc_ALL_finngen_results.tsv"), sep="\t")
nov = pd.read_csv(os.path.join(GC, "novelty_engine_ranked_ALL.tsv"), sep="\t")
il = pd.read_csv(os.path.join(GC, "intelligence_layer_final_table_ALL.tsv"), sep="\t")
uk = pd.read_csv(os.path.join(GC, "uk_panphenome_concordance_ALL.tsv"), sep="\t")


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", p, flush=True)


# ---------------------------------------------------------------- A ------- #
fig, ax = plt.subplots(figsize=(7.4, 4.2))
layers = ["cis-MR\n(causal hits)", "Colocalization", "Protein pQTL-MR",
          "Protein pQTL\ncoloc", "Novelty engine", "Intelligence\nlayer",
          "Two-population\nUK replication"]
after = [hits.phenotype.nunique(), coloc.disease.nunique(), pmr.disease.nunique(),
         pcol.disease.nunique(), nov.disease.nunique(), il.Disease.nunique(),
         uk.disease.nunique()]
before = [28, 24, 13, 13, 25, 25, 48]
x = np.arange(len(layers))
ax.bar(x - 0.2, before, 0.4, label="previous release (curated core)", color="#c9c9c9",
       edgecolor="#222")
ax.bar(x + 0.2, after, 0.4, label="this release (phenome-wide)", color="#2b6cb0",
       edgecolor="#222")
for xi, (b, a) in enumerate(zip(before, after)):
    ax.text(xi - 0.2, b + 4, str(b), ha="center", fontsize=7.5)
    ax.text(xi + 0.2, a + 4, str(a), ha="center", fontsize=7.5, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(layers, fontsize=7.6)
ax.set_ylabel("diseases covered")
ax.set_title("Every evidence layer now spans the FinnGen R12 phenome, not 28 curated diseases",
             fontsize=10, fontweight="bold")
ax.legend(fontsize=8, frameon=False)
save(fig, "Figure17_layer_disease_coverage.png")

# ---------------------------------------------------------------- B ------- #
fig, ax = plt.subplots(figsize=(6.4, 4.6))
d = pmr.dropna(subset=["OR", "MR_p"]).copy()
d["lnOR"] = np.log(d.OR.clip(lower=1e-6))
d["nlp"] = -np.log10(d.MR_p.clip(lower=1e-300))
sig = d.FDR < 0.05
ax.scatter(d.lnOR[~sig], d.nlp[~sig], s=9, color="#bbb", label="ns")
ax.scatter(d.lnOR[sig], d.nlp[sig], s=14, color="#c0392b", label="FDR<0.05")
for _, r in d[sig].nlargest(8, "nlp").iterrows():
    ax.annotate(f"{r.gene}→{str(r.disease)[:22]}", (r.lnOR, r.nlp), fontsize=6.4,
                xytext=(3, 3), textcoords="offset points")
ax.axvline(0, color="#444", lw=0.8, ls="--")
ax.set_xlabel("protein-level MR effect, ln(OR) per SD plasma protein")
ax.set_ylabel("-log10 P")
ax.set_xlim(np.percentile(d.lnOR, 0.5), np.percentile(d.lnOR, 99.5))
ax.set_title(f"Pan-phenome protein-level (INTERVAL cis-pQTL) MR\n"
             f"{len(d)} tests · {int(sig.sum())} FDR<0.05 · {d.disease.nunique()} diseases · "
             f"{d.gene.nunique()} proteins", fontsize=9.5, fontweight="bold")
ax.legend(fontsize=8, frameon=False)
save(fig, "Figure18_panphenome_pqtl_volcano.png")

# ---------------------------------------------------------------- C ------- #
m = nov.dropna(subset=["pqtl_OR"]).copy()
fig, ax = plt.subplots(figsize=(5.4, 5.0))
xe = np.log(m.OR.clip(lower=1e-6)); yp = np.log(m.pqtl_OR.clip(lower=1e-6))
conc = ((xe > 0) & (yp > 0)) | ((xe < 0) & (yp < 0))
ax.scatter(xe[~conc], yp[~conc], s=16, color="#e08a8a", label="discordant")
ax.scatter(xe[conc], yp[conc], s=16, color="#1f6f3f", label="concordant")
lim = max(abs(xe).max(), abs(yp).max()) * 1.05
ax.plot([-lim, lim], [-lim, lim], color="#888", lw=0.8, ls="--")
ax.axhline(0, color="#444", lw=0.6); ax.axvline(0, color="#444", lw=0.6)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("transcript-level cis-eQTL MR, ln(OR)")
ax.set_ylabel("protein-level cis-pQTL MR, ln(OR)")
ax.set_title(f"Transcript vs protein causal direction\n{int(conc.sum())}/{len(m)} "
             f"({100 * conc.mean():.0f}%) same direction", fontsize=9.5, fontweight="bold")
ax.legend(fontsize=8, frameon=False)
save(fig, "Figure19_eqtl_pqtl_concordance_all.png")

# ---------------------------------------------------------------- D ------- #
sysmap = nov.set_index(["gene_symbol", "disease"]).disease_system.to_dict()
il2 = il.copy()
il2["System"] = [sysmap.get((g.upper(), d), "Other")
                 for g, d in zip(il2.Gene, il2.Disease)]
ct = (il2.groupby(["System", "Novelty_tier"]).size().unstack(fill_value=0)
      .sort_values(by=list(range(1, 6))[::-1], ascending=False))
fig, ax = plt.subplots(figsize=(7.2, 4.4))
colors = {1: "#8e8e8e", 2: "#d9b382", 3: "#9ec5e8", 4: "#2b6cb0", 5: "#c0392b"}
bottom = np.zeros(len(ct))
for t in sorted(ct.columns):
    ax.bar(ct.index, ct[t], bottom=bottom, color=colors.get(t, "#ccc"),
           edgecolor="#222", linewidth=0.4, label=f"Tier {t}")
    bottom += ct[t].values
ax.set_ylabel("causal gene–disease pairs")
ax.set_xticklabels(ct.index, rotation=35, ha="right", fontsize=7.6)
ax.set_title(f"Evidence tiers across all {len(il)} causal pairs and {il.Disease.nunique()} diseases",
             fontsize=10, fontweight="bold")
ax.legend(fontsize=7.6, frameon=False, ncol=5)
save(fig, "Figure20_tiers_by_system_all.png")

# ---------------------------------------------------------------- E ------- #
top = nov.nlargest(28, "novelty_priority").iloc[::-1]
fig, ax = plt.subplots(figsize=(7.6, 6.4))
cmap = {"NOVEL protein-confirmed": "#c0392b", "NOVEL colocalized": "#2b6cb0",
        "novel nomination": "#9ec5e8", "recovered known axis": "#8e8e8e",
        "MHC-caution": "#d9b382"}
cols = [cmap.get(c, "#ccc") for c in top.category_label]
lbl = [f"{g} → {str(d)[:38]}" for g, d in zip(top.gene_symbol, top.disease)]
ax.barh(range(len(top)), top.novelty_priority, color=cols, edgecolor="#222", linewidth=0.4)
ax.set_yticks(range(len(top))); ax.set_yticklabels(lbl, fontsize=7.2)
ax.set_xlabel("novelty-priority score")
ax.set_title("Top pan-phenome novelty-priority targets", fontsize=10, fontweight="bold")
handles = [plt.Rectangle((0, 0), 1, 1, color=v, ec="#222") for v in cmap.values()]
ax.legend(handles, cmap.keys(), fontsize=7.2, frameon=False, loc="lower right")
save(fig, "Figure21_novelty_priority_all.png")

# ---------------------------------------------------------------- F ------- #
best = (pcol.dropna(subset=["PP_H4"]).sort_values("PP_H4", ascending=False)
        .drop_duplicates(["gene", "disease"]).head(25).iloc[::-1])
fig, ax = plt.subplots(figsize=(7.4, 6.0))
ax.barh(range(len(best)), best.PP_H4, color="#1f6f3f", edgecolor="#222", linewidth=0.4)
ax.set_yticks(range(len(best)))
ax.set_yticklabels([f"{g} → {str(d)[:40]}" for g, d in zip(best.gene, best.disease)],
                   fontsize=7.2)
ax.axvline(0.8, color="#c0392b", ls="--", lw=1)
ax.set_xlim(0, 1.02); ax.set_xlabel("protein-level colocalization PP.H4")
ax.set_title(f"Strongest plasma-protein colocalizations across the phenome\n"
             f"{len(pcol)} loci tested · {int((pcol.PP_H4 >= 0.8).sum())} with PP.H4≥0.8 · "
             f"{pcol.disease.nunique()} diseases", fontsize=9.5, fontweight="bold")
save(fig, "Figure22_pqtl_coloc_all.png")

print(f"\n[51] done — 6 panels in {FIG}")
