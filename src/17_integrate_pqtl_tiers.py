#!/usr/bin/env python
"""
HDDM Layer 54 - Step 17
Integrate the protein-level (INTERVAL pQTL) evidence into the claim-strength
ladder, producing the FINAL evidence table + a protein-confirmation figure.

Protein-level rule (applied ON TOP of the eQTL tier):
  tier 5 "protein-level causal target": eQTL tier>=4 AND pQTL-MR concordant
        direction AND pQTL coloc PP.H4>=0.8
  "+ protein-supported": concordant AND pQTL-MR nominal p<0.05
  "transcript/protein DISCORDANT (caution)": pQTL available but opposite sign
  else keep eQTL tier / "transcript-level only"
Outputs:
  06_genetic_causality/FINAL_evidence_tiers.tsv
  08_figures/nature/Figure7_pqtl_confirmation.png
"""
import os
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=r"I:\Plasma immune atalas"
GEN=os.path.join(ROOT,"06_genetic_causality")
FIG=os.path.join(ROOT,"08_figures","nature")

m=pd.read_csv(os.path.join(GEN,"evidence_tiered_with_pqtl.tsv"),sep="\t")

def concordant(r):
    if pd.isna(r.pQTL_OR): return None
    return (r.OR-1)*(r.pQTL_OR-1)>0

def final_label(r):
    c=concordant(r)
    base_tier=int(r.evidence_tier)
    if c is None:
        return base_tier, r.claim_label, "no plasma pQTL in INTERVAL (transcript-level only)"
    if not c:
        return base_tier, r.claim_label+" \u2014 transcript/protein DISCORDANT (caution)", \
               "eQTL and pQTL effects oppose \u2014 regulation may differ (e.g. soluble receptor)"
    # concordant
    if base_tier>=4 and pd.notna(r.pQTL_PPH4) and r.pQTL_PPH4>=0.8:
        return 5, "protein-level causal target (pQTL-MR + pQTL-coloc)", \
               "colocalized at BOTH transcript and protein level, concordant direction"
    if pd.notna(r.pQTL_MR_p) and r.pQTL_MR_p<0.05:
        return base_tier, r.claim_label+" \u2014 protein-supported (concordant pQTL-MR)", \
               "direction concordant and pQTL-MR nominally significant"
    return base_tier, r.claim_label+" \u2014 protein direction concordant", \
           "pQTL direction concordant but not independently significant"

res=m.apply(final_label,axis=1,result_type="expand")
m["final_tier"]=res[0]; m["final_claim"]=res[1]; m["protein_note"]=res[2]
m["pQTL_concordant"]=m.apply(lambda r: "" if concordant(r) is None else ("yes" if concordant(r) else "NO"),axis=1)

cols=["gene_symbol","disease","immune_class","chr","OR","OR_l95","OR_u95","FDR","PP_H4",
      "pQTL_OR","pQTL_MR_p","pQTL_PPH4","pQTL_concordant","final_tier","final_claim","protein_note"]
out=m.sort_values(["final_tier","FDR"],ascending=[False,True])
out[cols].to_csv(os.path.join(GEN,"FINAL_evidence_tiers.tsv"),sep="\t",index=False)

print("=== FINAL evidence tiers ===")
print(out[["gene_symbol","disease","OR","pQTL_OR","pQTL_PPH4","pQTL_concordant","final_tier","final_claim"]].to_string(index=False))
print("\nfinal tier counts:"); print(out.final_tier.value_counts().sort_index().to_string())
t5=out[out.final_tier==5]
print("\nTIER 5 protein-level causal targets:", list(t5.gene_symbol+"\u2192"+t5.disease))

# ---------- Figure 7 ----------
avail=m[m.pQTL_OR.notna()].copy()
avail["lab"]=avail.gene_symbol+"\u2192"+avail.disease.str[:12]
fig=plt.figure(figsize=(13,7.5))
gs=fig.add_gridspec(1,2,width_ratios=[1.15,1],wspace=0.32)

# A: eQTL vs pQTL effect concordance (log-OR scatter)
ax=fig.add_subplot(gs[0,0])
x=np.log(avail.OR); y=np.log(avail.pQTL_OR)
col=["#2e75b6" if (a-0)*(b-0)>0 else "#c0392b" for a,b in zip(x,y)]
ax.axhline(0,c="#999",lw=.8); ax.axvline(0,c="#999",lw=.8)
ax.scatter(x,y,c=col,s=70,edgecolor="k",lw=.5,zorder=3)
for _,r in avail.iterrows():
    ax.annotate(r.gene_symbol,(np.log(r.OR),np.log(r.pQTL_OR)),fontsize=6.5,
                xytext=(3,3),textcoords="offset points")
lim=max(abs(x).max(),abs(y).max())*1.15
ax.plot([-lim,lim],[-lim,lim],ls="--",c="#2e75b6",lw=.8,alpha=.6)
ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.set_xlabel("transcript (cis-eQTL) MR log-OR"); ax.set_ylabel("protein (cis-pQTL) MR log-OR")
ax.set_title("a  Transcript vs protein MR concordance\n(blue=concordant, red=discordant)",fontsize=10,loc="left")

# B: protein-level confirmation bars (pQTL coloc PP.H4 for eQTL hits)
ax=fig.add_subplot(gs[0,1])
b=avail.sort_values("pQTL_PPH4")
colors=["#1f4e79" if p>=0.8 else ("#5b9bd5" if p>=0.5 else "#cccccc") for p in b.pQTL_PPH4.fillna(0)]
ax.barh(b.lab,b.pQTL_PPH4.fillna(0),color=colors)
ax.axvline(0.8,ls="--",c="green",lw=1)
ax.set_xlabel("protein-level (pQTL) coloc PP.H4")
ax.set_title("b  Protein-level colocalization\nof transcript-level hits",fontsize=10,loc="left")
ax.tick_params(axis="y",labelsize=7)
fig.suptitle("Figure 7  |  Protein-level (INTERVAL plasma pQTL) validation of genetic nominations",
             fontsize=12.5,fontweight="bold",x=0.02,ha="left")
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(os.path.join(FIG,"Figure7_pqtl_confirmation.png"),dpi=150,bbox_inches="tight")
print("\nwrote Figure7_pqtl_confirmation.png")
