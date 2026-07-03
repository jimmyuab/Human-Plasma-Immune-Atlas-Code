#!/usr/bin/env python
"""
HDDM Layer 54 - Step 21  (NOVELTY LAYER)
Adds four publication-strengthening analyses on top of the evidence arc, all
from local data / literature-curated annotation (no downloads, no fabrication):

  1. Immune-CLASS enrichment of genetically-supported causal targets
     (Fisher exact vs the curated plasma-immune background).
  2. Cross-disease PLEIOTROPY map of shared causal immune axes
     (genes acting on >=2 diseases, with direction-of-effect).
  3. Genetic support for THERAPEUTIC DIRECTION: OR<1 => protein protective
     => agonize/replace ; OR>1 => risk => block ; cross-checked against
     literature-curated approved/known drug mechanism (explicitly labelled).
  4. NOVELTY prioritisation: separates atlas hits into (a) recovered known
     therapeutic axes (internal validity) and (b) novel colocalised+replicated
     nominations without an approved drug for that indication.

Outputs:
  06_genetic_causality/novelty_class_enrichment.tsv
  06_genetic_causality/novelty_pleiotropy.tsv
  06_genetic_causality/novelty_drug_direction.tsv
  06_genetic_causality/novelty_prioritised_targets.tsv
  08_figures/nature/Figure9_novelty.png
"""
import os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=r"I:\Plasma immune atalas"
GEN=os.path.join(ROOT,"06_genetic_causality"); PROC=os.path.join(ROOT,"02_data_processed")
FIG=os.path.join(ROOT,"08_figures","nature")

m=pd.read_csv(os.path.join(GEN,"FINAL_evidence_tiers_repl.tsv"),sep="\t")
ann=pd.read_csv(os.path.join(PROC,"plasma_immune_protein_annotation.tsv"),sep="\t")
imm=ann[ann.is_plasma_immune==1].copy()

# map gene -> immune_class
g2class=dict(zip(imm.gene_symbol,imm.immune_class))
m["immune_class"]=m["immune_class"].fillna(m.gene_symbol.map(g2class))

# ============================================================
# 1. IMMUNE-CLASS ENRICHMENT of causal targets (tier>=4, unique genes)
# ============================================================
causal_genes=set(m[m.final_tier>=4].gene_symbol)
bg=imm[["gene_symbol","immune_class"]].dropna().drop_duplicates("gene_symbol")
Ntot=bg.gene_symbol.nunique(); Ncaus=len(causal_genes & set(bg.gene_symbol))
rows=[]
for cls,sub in bg.groupby("immune_class"):
    cls_genes=set(sub.gene_symbol)
    a=len(cls_genes & causal_genes)          # in class & causal
    b=len(cls_genes)-a                        # in class not causal
    c=Ncaus-a                                 # causal not in class
    d=Ntot-len(cls_genes)-c                   # neither
    if a==0: continue
    orr,p=fisher_exact([[a,b],[c,d]],alternative="greater")
    frac_cls=a/len(cls_genes)
    rows.append(dict(immune_class=cls,n_class=len(cls_genes),n_causal_in_class=a,
                     class_hit_rate=frac_cls,odds_ratio=orr,p=p))
enr=pd.DataFrame(rows).sort_values("p")
from statsmodels.stats.multitest import multipletests
enr["FDR"]=multipletests(enr.p,method="fdr_bh")[1]
enr.to_csv(os.path.join(GEN,"novelty_class_enrichment.tsv"),sep="\t",index=False)
print("=== immune-class enrichment (Fisher, causal tier>=4) ===")
print(enr.to_string(index=False))

# ============================================================
# 2. PLEIOTROPY: genes on >=2 diseases
# ============================================================
gc=m.gene_symbol.value_counts(); pleio_genes=gc[gc>=2].index.tolist()
pl=m[m.gene_symbol.isin(pleio_genes)][["gene_symbol","disease","OR","final_tier","rep_status"]].copy()
pl["direction"]=np.where(pl.OR>1,"risk (block)","protective (agonise)")
pl=pl.sort_values(["gene_symbol","disease"])
pl.to_csv(os.path.join(GEN,"novelty_pleiotropy.tsv"),sep="\t",index=False)
print("\n=== pleiotropic causal immune axes (>=2 diseases) ===")
print(pl.to_string(index=False))

# ============================================================
# 3. THERAPEUTIC DIRECTION vs literature-curated known drug mechanism
#    (curated = established pharmacology from the primary literature;
#     labelled so no computed claim is implied)
# ============================================================
# established target -> (mechanism of approved/known drug, direction the drug pushes)
KNOWN={
 "IL6ST":("tocilizumab / sarilumab (IL-6R axis)","block"),
 "CTLA4":("abatacept (CTLA4-Ig, co-stim agonist)","agonise"),
 "TNFRSF1A":("etanercept / anti-TNF","block"),
 "TNFSF14":("LIGHT axis (experimental)","block"),
 "TNFRSF14":("HVEM/LIGHT axis (experimental)","block"),
 "IL4":("dupilumab (IL-4Ra)","block"),
 "ERBB3":("(oncology, not immune-approved)","na"),
}
rows=[]
for _,r in m[m.final_tier>=3].iterrows():
    g=r.gene_symbol
    gen_dir="block" if r.OR>1 else "agonise"
    known=KNOWN.get(g)
    if known:
        drug,drug_dir=known
        agree=("match" if (drug_dir==gen_dir) else
               ("informative-discordance" if drug_dir in("block","agonise") else "na"))
    else:
        drug,drug_dir,agree="(no approved immune drug)","-","novel"
    rows.append(dict(gene_symbol=g,disease=r.disease,OR=round(r.OR,3),
                     genetic_direction=gen_dir,known_drug=drug,drug_direction=drug_dir,
                     concordance=agree,final_tier=int(r.final_tier),rep_status=r.rep_status))
dd=pd.DataFrame(rows)
dd.to_csv(os.path.join(GEN,"novelty_drug_direction.tsv"),sep="\t",index=False)
print("\n=== genetic support for therapeutic direction ===")
print(dd.to_string(index=False))

# ============================================================
# 4. NOVELTY prioritisation
#    known-axis (internal positive control) vs NOVEL nomination
# ============================================================
known_genes=set(KNOWN)
pri=m[m.final_tier>=4].copy()
pri["category"]=np.where(pri.gene_symbol.isin(known_genes),
                         "recovered known axis","NOVEL nomination")
# novelty score: coloc + pQTL protein support + replication + non-MHC
def score(r):
    s=0.0
    s+=float(pd.notna(r.PP_H4) and r.PP_H4>=0.8)          # transcript coloc
    s+=float(str(r.pQTL_concordant)=="yes")               # protein concordant
    s+=float(pd.notna(r.pQTL_PPH4) and r.pQTL_PPH4>=0.8)  # protein coloc
    s+=float(r.rep_status=="replicated")                  # independent repl
    s+=0.5*float(r.final_tier==5)
    return s
pri["evidence_score"]=pri.apply(score,axis=1)
pri=pri.sort_values(["category","evidence_score"],ascending=[True,False])
pri_out=pri[["gene_symbol","disease","category","final_tier","PP_H4","pQTL_PPH4",
             "pQTL_concordant","rep_status","evidence_score"]]
pri_out.to_csv(os.path.join(GEN,"novelty_prioritised_targets.tsv"),sep="\t",index=False)
print("\n=== prioritised targets (known vs novel) ===")
print(pri_out.to_string(index=False))

# ============================================================
# FIGURE 9 (multi-panel novelty)
# ============================================================
fig=plt.figure(figsize=(15,10))
gs=fig.add_gridspec(2,2,hspace=0.42,wspace=0.28)

# a: class enrichment
ax=fig.add_subplot(gs[0,0])
e=enr.sort_values("odds_ratio")
y=np.arange(len(e))
sig=e.FDR<0.05
ax.barh(y,-np.log10(e.p),color=["#c0392b" if s else "#5b9bd5" for s in sig])
ax.set_yticks(y); ax.set_yticklabels([f"{c}  ({int(a)}/{int(n)})" for c,a,n in
              zip(e.immune_class,e.n_causal_in_class,e.n_class)],fontsize=8)
ax.axvline(-np.log10(0.05),ls="--",c="green",lw=.8)
ax.set_xlabel("-log10(P)  Fisher enrichment among causal targets")
ax.set_title("a  Immune-class enrichment of tier\u22654 causal targets",fontsize=10.5,loc="left")

# b: pleiotropy bipartite
ax=fig.add_subplot(gs[0,1])
genes=sorted(pl.gene_symbol.unique()); diss=sorted(pl.disease.unique())
gy={g:i for i,g in enumerate(genes)}; dy={d:i for i,d in enumerate(diss)}
gx,dx=0,1
for _,r in pl.iterrows():
    col="#c0392b" if r.OR>1 else "#1f4e79"
    ax.plot([gx,dx],[gy[r.gene_symbol],dy[r.disease]*len(genes)/max(1,len(diss))],
            "-",c=col,lw=1.4,alpha=0.75)
for g,i in gy.items(): ax.text(gx-0.03,i,g,ha="right",va="center",fontsize=8,fontweight="bold")
for d,i in dy.items(): ax.text(dx+0.03,i*len(genes)/max(1,len(diss)),d[:16],ha="left",va="center",fontsize=7.5)
ax.set_xlim(-0.5,1.6); ax.axis("off")
ax.set_title("b  Pleiotropic shared causal axes\n(red=risk/block, blue=protective/agonise)",fontsize=10.5,loc="left")

# c: therapeutic direction concordance
ax=fig.add_subplot(gs[1,0]); ax.axis("off")
kd=dd[dd.concordance.isin(["match","informative-discordance"])]
cell=[]
for _,r in kd.iterrows():
    cell.append([r.gene_symbol,r.disease[:14],r.genetic_direction,r.drug_direction,r.concordance])
tb=ax.table(cellText=cell,colLabels=["gene","disease","genetic","drug","verdict"],
            loc="center",cellLoc="center")
tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1,1.5)
for j in range(5): tb[0,j].set_facecolor("#1f4e79"); tb[0,j].set_text_props(color="w",fontweight="bold")
for i,(_,r) in enumerate(kd.iterrows(),1):
    c="#d5f5e3" if r.concordance=="match" else "#fdebd0"
    for j in range(5): tb[i,j].set_facecolor(c)
ax.set_title("c  Genetics recover approved-drug DIRECTION",fontsize=10.5,loc="left")

# d: novelty prioritisation scatter
ax=fig.add_subplot(gs[1,1])
for cat,mk,cc in [("recovered known axis","s","#888"),("NOVEL nomination","o","#c0392b")]:
    s=pri[pri.category==cat]
    ax.scatter(s.evidence_score+np.random.RandomState(0).uniform(-.05,.05,len(s)),
               range(len(s)) if False else s.final_tier+np.random.RandomState(1).uniform(-.08,.08,len(s)),
               marker=mk,s=90,c=cc,edgecolor="k",lw=.4,label=cat,alpha=0.85)
for _,r in pri[pri.category=="NOVEL nomination"].iterrows():
    ax.annotate(f"{r.gene_symbol}\u2192{r.disease[:8]}",(r.evidence_score,r.final_tier),
                fontsize=6.5,xytext=(3,3),textcoords="offset points")
ax.set_xlabel("cumulative evidence score (coloc+pQTL+repl)")
ax.set_ylabel("final tier"); ax.set_yticks([4,5]); ax.legend(fontsize=8,loc="lower right")
ax.set_title("d  Novel vs recovered-known target prioritisation",fontsize=10.5,loc="left")

fig.suptitle("Figure 9  |  Novelty: immune-class enrichment, pleiotropy, genetic drug-direction support & prioritisation",
             fontsize=12.5,fontweight="bold",x=0.02,ha="left")
fig.savefig(os.path.join(FIG,"Figure9_novelty.png"),dpi=150,bbox_inches="tight")
print("\nwrote Figure9_novelty.png")
