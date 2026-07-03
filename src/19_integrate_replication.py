#!/usr/bin/env python
"""
HDDM Layer 54 - Step 19
Fold the independent-GWAS replication (OpenGWAS) into the FINAL evidence table
and draw Figure 8 (replication). Replication status per hit:
  replicated     : found, directionally concordant, p<0.05
  concordant     : found, concordant, not significant
  discordant     : found, opposite direction
  not covered    : SNP absent from independent GWAS (or none exists)
Outputs:
  06_genetic_causality/FINAL_evidence_tiers_repl.tsv
  08_figures/nature/Figure8_replication.png
"""
import os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=r"I:\Plasma immune atalas"
GEN=os.path.join(ROOT,"06_genetic_causality"); FIG=os.path.join(ROOT,"08_figures","nature")

fin=pd.read_csv(os.path.join(GEN,"FINAL_evidence_tiers.tsv"),sep="\t")
rep=pd.read_csv(os.path.join(GEN,"opengwas_replication.tsv"),sep="\t")

def status(r):
    if not bool(r.rep_found): return "not covered"
    if r.rep_concordant=="yes" and bool(r.rep_sig): return "replicated"
    if r.rep_concordant=="yes": return "concordant (ns)"
    if r.rep_concordant=="NO": return "discordant"
    return "not covered"
rep["rep_status"]=rep.apply(status,axis=1)

m=fin.merge(rep[["gene","disease","rep_gwas","rep_cite","rep_p","rep_status"]].rename(
    columns={"gene":"gene_symbol"}),on=["gene_symbol","disease"],how="left")
m["rep_status"]=m["rep_status"].fillna("not covered")
m.to_csv(os.path.join(GEN,"FINAL_evidence_tiers_repl.tsv"),sep="\t",index=False)

n_rep=int((m.rep_status=="replicated").sum())
n_conc=int((m.rep_status.isin(["replicated","concordant (ns)"])).sum())
n_cov=int((m.rep_status!="not covered").sum())
print("=== FINAL with replication ===")
print(m[["gene_symbol","disease","final_tier","PP_H4","pQTL_PPH4","rep_status","rep_p"]].to_string(index=False))
print(f"\ncovered by independent GWAS: {n_cov}/{len(m)} | replicated(p<0.05): {n_rep} | concordant: {n_conc}")

# -------- Figure 8 --------
cov=m[m.rep_status!="not covered"].copy()
cov["lab"]=cov.gene_symbol+"\u2192"+cov.disease.str[:12]
cov["mlog"]=-np.log10(cov.rep_p.clip(lower=1e-300))
cov=cov.sort_values("mlog")
cmap={"replicated":"#1f4e79","concordant (ns)":"#5b9bd5","discordant":"#c0392b"}
colors=[cmap.get(s,"#ccc") for s in cov.rep_status]

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,7),gridspec_kw={"width_ratios":[1.3,1]})
ax1.barh(cov.lab,cov.mlog,color=colors)
ax1.axvline(-np.log10(0.05),ls="--",c="green",lw=1,label="p=0.05")
ax1.set_xlabel("independent-GWAS replication  \u2212log10(P)")
ax1.set_title("a  Replication in independent consortium GWAS\n(color = status; direction concordant unless red)",
              fontsize=10,loc="left"); ax1.tick_params(axis="y",labelsize=7.5)
ax1.legend(fontsize=8,loc="lower right")
# annotate cohort at bar end
for y,(_,r) in enumerate(cov.iterrows()):
    ax1.text(r.mlog+0.5,y,r.rep_cite,fontsize=6,va="center")

# panel b: summary donut of replication outcome among covered hits
ax2.axis("equal")
vc=m.rep_status.value_counts()
order=["replicated","concordant (ns)","discordant","not covered"]
vals=[int(vc.get(k,0)) for k in order]
cols=[cmap.get(k,"#dddddd") for k in order]
w,_=ax2.pie(vals,colors=cols,startangle=90,wedgeprops=dict(width=0.42,edgecolor="w"))
ax2.legend([f"{k} ({v})" for k,v in zip(order,vals)],loc="center",fontsize=8,frameon=False)
ax2.set_title("b  Replication outcome across all\nsignificant gene\u2013disease hits",fontsize=10,loc="left")
fig.suptitle("Figure 8  |  Independent replication of genetic nominations (OpenGWAS, non-FinnGen)",
             fontsize=12.5,fontweight="bold",x=0.02,ha="left")
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(os.path.join(FIG,"Figure8_replication.png"),dpi=150,bbox_inches="tight")
print("wrote Figure8_replication.png")
