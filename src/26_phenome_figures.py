#!/usr/bin/env python
"""
HDDM Layer 54 - Step 26
Pan-phenome FIGURE SUITE. Turns the expanded immune cis-MR + analysis into a
large, publishable figure set. All real-data. Two destinations:
  08_figures/nature/        -> Figure11..Figure14 (main, multi-panel)
  08_figures/phenome/       -> per-disease + per-category supplementary panels
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality")
NAT  = os.path.join(ROOT, "08_figures", "nature")
PHE  = os.path.join(ROOT, "08_figures", "phenome")
os.makedirs(PHE, exist_ok=True)

mr    = pd.read_csv(os.path.join(GEN, "cis_MR_phenome_results.tsv"), sep="\t")
sig   = pd.read_csv(os.path.join(GEN, "phenome_hits.tsv"), sep="\t")
pleio = pd.read_csv(os.path.join(GEN, "phenome_pleiotropy_axes.tsv"), sep="\t")
cell  = pd.read_csv(os.path.join(GEN, "phenome_cellsource_map.tsv"), sep="\t")
dm    = pd.read_csv(os.path.join(GEN, "phenome_direction_map.tsv"), sep="\t")

CATCOL = {"Autoimmune":"#c0392b","Cardiovascular":"#1f4e79","Metabolic":"#e67e22",
          "Renal":"#16a085","Neuro/Aging":"#8e44ad","Other":"#888888"}
n = 0
def sv(fig, path, name):
    global n; n += 1
    fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(os.path.join(path,name),dpi=150,bbox_inches="tight")
    plt.close(fig)

# ================= FIGURE 11: pan-phenome global volcano =================
fig, ax = plt.subplots(figsize=(11,7))
for cat, d in mr.groupby("disease_category"):
    ns = d[d.FDR >= 0.05]
    ax.scatter(np.log(ns.OR), -np.log10(ns.MR_p.clip(lower=1e-300)), s=6,
               c="#dddddd", alpha=0.5, zorder=1)
for cat, d in mr[mr.FDR<0.05].groupby("disease_category"):
    ax.scatter(np.log(d.OR), -np.log10(d.MR_p.clip(lower=1e-300)), s=34,
               c=CATCOL.get(cat,"#888"), edgecolor="k", lw=.3, label=cat, zorder=3)
top = mr[mr.FDR<0.05].sort_values("MR_p").head(22)
for _, r in top.iterrows():
    ax.annotate(f"{r.gene_symbol}·{r.disease[:10]}", (np.log(r.OR), -np.log10(max(r.MR_p,1e-300))),
                fontsize=6.3, xytext=(3,2), textcoords="offset points")
ax.axvline(0, ls="--", c="#888", lw=.7)
ax.set_xlabel("cis-MR log(OR) per SD expression"); ax.set_ylabel("-log10(P)")
ax.legend(fontsize=8, title="disease category")
fig.suptitle("Figure 11  |  Pan-phenome immune cis-MR: causal immune proteins across "
             "autoimmune, cardiovascular, metabolic, renal & aging disease",
             fontsize=12, fontweight="bold", x=0.02, ha="left")
sv(fig, NAT, "Figure11_phenome_volcano.png")

# ================= FIGURE 12: cross-category pleiotropy heatmap =================
if len(pleio):
    genes = pleio.gene_symbol.unique().tolist()
    diss  = sig[sig.gene_symbol.isin(genes)].disease.unique().tolist()
    M = pd.DataFrame(index=genes, columns=diss, dtype=float)
    for _, r in sig[sig.gene_symbol.isin(genes)].iterrows():
        M.loc[r.gene_symbol, r.disease] = np.log(r.OR)
    fig, ax = plt.subplots(figsize=(max(8,0.5*len(diss)+3), max(4,0.4*len(genes)+2)))
    vmax = np.nanmax(np.abs(M.values.astype(float)))
    im = ax.imshow(M.values.astype(float), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(diss))); ax.set_xticklabels(diss, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=8)
    for i in range(len(genes)):
        for j in range(len(diss)):
            v = M.values[i,j]
            if not np.isnan(v): ax.text(j,i,"+" if v>0 else "\u2212",ha="center",va="center",
                                        fontsize=8,color="w",fontweight="bold")
    cb = fig.colorbar(im, ax=ax, shrink=0.7); cb.set_label("cis-MR log(OR)")
    dcat = sig.drop_duplicates("disease").set_index("disease").disease_category
    for j,dd in enumerate(diss):
        ax.get_xticklabels()[j].set_color(CATCOL.get(dcat.get(dd,"Other"),"#000"))
    fig.suptitle("Figure 12  |  Shared causal immune axes across disease categories "
                 "(genes causal in \u22652 categories; +=risk, \u2212=protective)",
                 fontsize=11.5, fontweight="bold", x=0.02, ha="left")
    sv(fig, NAT, "Figure12_pleiotropy_heatmap.png")

# ================= FIGURE 13: cell-source causal map =================
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,6),gridspec_kw={"width_ratios":[1,1]})
c = cell.sort_values("odds_ratio")
y = np.arange(len(c))
ax1.barh(y, -np.log10(c.p), color=["#c0392b" if p<0.05 else "#5b9bd5" for p in c.p])
ax1.set_yticks(y); ax1.set_yticklabels([f"{s} ({int(a)}/{int(nn)})" for s,a,nn in
                 zip(c.source_cell,c.n_causal,c.n_source)], fontsize=8)
ax1.axvline(-np.log10(0.05), ls="--", c="green", lw=.8)
ax1.set_xlabel("-log10(P) Fisher enrichment"); ax1.set_title("a  Causal-target enrichment by blood-cell source", fontsize=10.5, loc="left")
# panel b: causal counts by category x source (stacked)
src_cat = sig.dropna(subset=["immune_source_cells"]).copy()
src_cat = src_cat.assign(src=src_cat.immune_source_cells.str.split(r"[;,]")).explode("src")
src_cat["src"]=src_cat.src.str.strip()
piv = src_cat.pivot_table(index="src", columns="disease_category", values="gene_symbol",
                          aggfunc="nunique", fill_value=0)
piv = piv.loc[piv.sum(1).sort_values().index]
bottom = np.zeros(len(piv))
for cat in piv.columns:
    ax2.barh(piv.index, piv[cat], left=bottom, color=CATCOL.get(cat,"#888"), label=cat)
    bottom += piv[cat].values
ax2.set_xlabel("causal immune targets (unique genes)"); ax2.legend(fontsize=7.5)
ax2.set_title("b  Cell-source of causal targets by disease category", fontsize=10.5, loc="left")
fig.suptitle("Figure 13  |  Cell-origin map of genetically-supported causal immune proteins",
             fontsize=12, fontweight="bold", x=0.02, ha="left")
sv(fig, NAT, "Figure13_cellsource_map.png")

# ================= FIGURE 14: direction-aware therapeutic map =================
fig, ax = plt.subplots(figsize=(11, max(5,0.32*len(dm)+1)))
dm2 = dm.sort_values(["disease_category","OR"]).reset_index(drop=True)
y = np.arange(len(dm2))
col = [CATCOL.get(c,"#888") for c in dm2.disease_category]
ax.scatter(np.log(dm2.OR), y, c=col, s=42, edgecolor="k", lw=.3, zorder=3)
ax.axvline(0, ls="--", c="red", lw=1)
for yi,(_,r) in enumerate(dm2.iterrows()):
    tag = "\u25b2block" if r.OR>1 else "\u25bcagonise"
    star = " *MHC" if str(r.status).startswith("MHC") else (" [known]" if r.status=="recovers known axis" else "")
    ax.text(np.log(r.OR)+(0.03 if r.OR>1 else -0.03), yi,
            f"{r.gene_symbol}\u2192{r.disease[:12]} {tag}{star}",
            va="center", ha="left" if r.OR>1 else "right", fontsize=6.5)
ax.set_yticks([]); ax.set_xlabel("cis-MR log(OR)  (>0 risk\u2192block ; <0 protective\u2192agonise)")
handles=[plt.Line2D([0],[0],marker='o',ls='',mfc=CATCOL[c],mec='k',label=c) for c in CATCOL if c in dm2.disease_category.values]
ax.legend(handles=handles, fontsize=8, loc="lower right")
fig.suptitle("Figure 14  |  Direction-aware therapeutic map across the phenome "
             "(genetic effect direction \u2192 modality)", fontsize=11.5, fontweight="bold", x=0.02, ha="left")
sv(fig, NAT, "Figure14_direction_map.png")

# ================= SUPPLEMENTARY: per-disease volcano + forest =================
for dis, d in mr.groupby("disease"):
    cat = d.disease_category.iloc[0]
    fig, ax = plt.subplots(figsize=(6.8,5))
    s = d.FDR<0.05
    ax.scatter(np.log(d.OR[~s]), -np.log10(d.MR_p[~s].clip(lower=1e-300)), s=9, c="#ccc")
    ax.scatter(np.log(d.OR[s]), -np.log10(d.MR_p[s].clip(lower=1e-300)), s=30,
               c=CATCOL.get(cat,"#888"), edgecolor="k", lw=.3)
    for _, r in d[s].sort_values("MR_p").head(8).iterrows():
        ax.annotate(r.gene_symbol,(np.log(r.OR),-np.log10(max(r.MR_p,1e-300))),
                    fontsize=7, xytext=(3,2), textcoords="offset points")
    ax.axvline(0, ls="--", c="#888", lw=.7)
    ax.set_xlabel("cis-MR log(OR)"); ax.set_ylabel("-log10(P)")
    ax.set_title(f"{dis}  [{cat}]  \u2014 {int(s.sum())} FDR<0.05", fontsize=10)
    safe = dis.replace(" ","_").replace("/","_")[:18]
    sv(fig, PHE, f"PFig_volcano_{safe}.png")

# per-category forest of top hits
for cat, d in sig.groupby("disease_category"):
    top = d.sort_values("MR_p").head(20).sort_values("OR")
    if len(top)<1: continue
    fig, ax = plt.subplots(figsize=(7, max(3,0.33*len(top)+1)))
    y = np.arange(len(top))
    ax.errorbar(top.OR, y, xerr=[top.OR-top.OR_l95, top.OR_u95-top.OR], fmt="o",
                ecolor="#bbb", capsize=2, ls="none", mfc="w", mec="w")
    ax.scatter(top.OR, y, c=CATCOL.get(cat,"#888"), s=32, zorder=3, edgecolor="k", lw=.3)
    ax.axvline(1, ls="--", c="red", lw=1); ax.set_yticks(y)
    ax.set_yticklabels([f"{g}\u2192{dd[:12]}" for g,dd in zip(top.gene_symbol,top.disease)], fontsize=7)
    ax.set_xscale("log"); ax.set_xlabel("OR per SD cis-expr (95% CI)")
    ax.set_title(f"Top causal immune targets \u2014 {cat}", fontsize=10)
    sv(fig, PHE, f"PFig_forest_{cat.replace('/','_')}.png")

# category summary bar
fig, ax = plt.subplots(figsize=(7,4.5))
cc = sig.disease_category.value_counts()
ax.bar(cc.index, cc.values, color=[CATCOL.get(c,"#888") for c in cc.index])
for i,v in enumerate(cc.values): ax.text(i,v+0.3,str(v),ha="center",fontsize=9)
ax.set_ylabel("FDR<0.05 gene-disease hits"); ax.set_title("Causal immune hits per disease category", fontsize=10)
sv(fig, PHE, "PFig_hits_per_category.png")

print(f"pan-phenome figures written: {n}")
print("  main (nature/): Figure11-14 ; supplementary (phenome/): rest")
print("  nature dir total:", len([f for f in os.listdir(NAT) if f.endswith('.png')]))
print("  phenome dir total:", len(os.listdir(PHE)))
