#!/usr/bin/env python
"""
HDDM Layer 54 - Step 29
Figures for the integrated novelty priority engine (src/28).
Writes Figure15 (ranking) + Figure16 (evidence-component & novel-vs-known) to
08_figures/nature/.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality"); NAT = os.path.join(ROOT, "08_figures", "nature")

r = pd.read_csv(os.path.join(GEN, "novelty_engine_ranked.tsv"), sep="\t")
CATCOL = {"Autoimmune":"#c0392b","Cardiovascular":"#1f4e79","Metabolic":"#e67e22",
          "Renal":"#16a085","Neuro/Aging":"#8e44ad","Other":"#888"}
LABCOL = {"NOVEL colocalized":"#c0392b","novel nomination":"#e67e22",
          "recovered known axis":"#888","MHC-caution":"#bdc3c7"}

# ---- Figure 15: top-30 ranked bar (stacked evidence components) ----
top = r.head(30).iloc[::-1]
comp = ["s_causal","s_coloc","s_pleio","s_drug","s_cell"]
compcol = {"s_causal":"#1f4e79","s_coloc":"#16a085","s_pleio":"#e67e22",
           "s_drug":"#8e44ad","s_cell":"#95a5a6"}
fig, ax = plt.subplots(figsize=(11, 9))
y = np.arange(len(top)); left = np.zeros(len(top))
for c in comp:
    ax.barh(y, top[c], left=left, color=compcol[c], label=c.replace("s_",""))
    left += top[c].values
# penalties as negative
for pcol,pc in [("p_known","#7f8c8d"),("p_mhc","#c0392b")]:
    ax.barh(y, -top[pcol], color=pc, alpha=0.5, label=pcol.replace("p_","-"))
ax.set_yticks(y); ax.set_yticklabels([f"{g}\u2192{d[:14]}" for g,d in zip(top.gene_symbol,top.disease)], fontsize=7.5)
ax.axvline(0, c="k", lw=.6)
for yi,(_,rr) in enumerate(top.iterrows()):
    ax.text(rr.novelty_priority+0.05, yi, f"{rr.novelty_priority:.2f}", va="center", fontsize=6.5)
ax.set_xlabel("novelty priority score (stacked evidence components; penalties negative)")
ax.legend(fontsize=7.5, ncol=4, loc="lower right")
fig.suptitle("Figure 15  |  Integrated novelty-priority ranking of causal immune targets across the phenome",
             fontsize=12, fontweight="bold", x=0.02, ha="left")
fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig(os.path.join(NAT,"Figure15_novelty_ranking.png"),dpi=150,bbox_inches="tight")
plt.close(fig)

# ---- Figure 16: 2-panel (evidence scatter + novel-vs-known) ----
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,6.5))
# a: coloc vs causal, size=druggability, color=category
for cat,d in r.groupby("disease_category"):
    ax1.scatter(d.s_causal, d.PP_H4.fillna(0), s=30+40*d.druggability, alpha=0.75,
                c=CATCOL.get(cat,"#888"), edgecolor="k", lw=.3, label=cat)
for _,rr in r[r.category_label=="NOVEL colocalized"].head(14).iterrows():
    ax1.annotate(f"{rr.gene_symbol}\u2192{rr.disease[:8]}",(rr.s_causal,rr.PP_H4),
                 fontsize=6.5, xytext=(3,2), textcoords="offset points")
ax1.axhline(0.8, ls="--", c="green", lw=.8)
ax1.set_xlabel("causal evidence (\u2212log10 FDR, scaled)"); ax1.set_ylabel("transcript coloc PP.H4")
ax1.legend(fontsize=7.5, title="category"); ax1.set_title("a  Evidence landscape (size = druggability)", fontsize=10.5, loc="left")
# b: novel vs known priority distribution
order=["NOVEL colocalized","novel nomination","recovered known axis","MHC-caution"]
data=[r[r.category_label==o].novelty_priority.values for o in order]
bp=ax2.boxplot(data, vert=False, patch_artist=True, widths=0.6)
for patch,o in zip(bp["boxes"],order): patch.set_facecolor(LABCOL[o])
ax2.set_yticklabels([f"{o}\n(n={len(r[r.category_label==o])})" for o in order], fontsize=8)
ax2.set_xlabel("novelty priority score")
ax2.set_title("b  Priority separates novel targets from controls & MHC", fontsize=10.5, loc="left")
fig.suptitle("Figure 16  |  Evidence composition and novel-versus-known target separation",
             fontsize=12, fontweight="bold", x=0.02, ha="left")
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(os.path.join(NAT,"Figure16_novelty_evidence.png"),dpi=150,bbox_inches="tight")
plt.close(fig)
print("wrote Figure15_novelty_ranking.png + Figure16_novelty_evidence.png")
