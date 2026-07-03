#!/usr/bin/env python
"""
HDDM Layer 54 - Step 20
Generate the SUPPLEMENTARY / EXTENDED DATA figure suite (~70 figures) that a
high-impact paper carries alongside the 8 main multi-panel figures. Every panel
answers a specific per-disease / per-gene / QC question. All from local data
(no downloads). Output: 08_figures/supplementary/SFig_*.png
"""
import os, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=r"I:\Plasma immune atalas"
GEN=os.path.join(ROOT,"06_genetic_causality")
PROC=os.path.join(ROOT,"02_data_processed")
CIS=os.path.join(ROOT,"01_data_raw","INTERVAL_pQTL","cis")
SUP=os.path.join(ROOT,"08_figures","supplementary")
os.makedirs(SUP,exist_ok=True)

mr=pd.read_csv(os.path.join(GEN,"cis_MR_immune_results.tsv"),sep="\t")
ann=pd.read_csv(os.path.join(PROC,"plasma_immune_protein_annotation.tsv"),sep="\t")
imm=ann[ann.is_plasma_immune==1]
inst=pd.read_csv(os.path.join(GEN,"immune_cis_eqtl_instruments.tsv"),sep="\t")
rep_path=os.path.join(GEN,"opengwas_replication.tsv")
rep=pd.read_csv(rep_path,sep="\t") if os.path.exists(rep_path) else None

n=0
def save(fig,name,title):
    global n; n+=1
    num=f"S{n:02d}"
    fig.suptitle(f"Supplementary Figure {num}  |  {title}",fontsize=10.5,
                 fontweight="bold",x=0.02,ha="left")
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(os.path.join(SUP,f"SFig{num}_{name}.png"),dpi=140,bbox_inches="tight")
    plt.close(fig)

# ================================================================
# GROUP A — QC / overview
# ================================================================
# A1 immune class composition
cls=imm.immune_class.value_counts()
fig,ax=plt.subplots(figsize=(7,4.5))
ax.barh(cls.index[::-1],cls.values[::-1],color="#5b9bd5")
ax.set_xlabel("plasma immune proteins");
for i,v in enumerate(cls.values[::-1]): ax.text(v+2,i,str(v),va="center",fontsize=8)
save(fig,"immune_class_composition","Immune-class composition of the curated plasma immunome")

# A2 source cell lineage
if "immune_source_cells" in imm:
    src=imm.immune_source_cells.dropna().str.split(r"[;,]").explode().str.strip()
    src=src[src!=""].value_counts().head(12)
    fig,ax=plt.subplots(figsize=(7,4.5))
    ax.barh(src.index[::-1],src.values[::-1],color="#2e75b6")
    ax.set_xlabel("proteins mapped to lineage")
    save(fig,"source_cell_lineage","Blood-cell lineage of origin for immune proteins")

# A3 instrument Z-score strength distribution
fig,ax=plt.subplots(figsize=(6.5,4.2))
ax.hist(inst.Zscore.abs(),bins=60,color="#1f4e79")
ax.axvline(5.45,ls="--",c="red",label="|Z|=5.45 (p=5e-8)")
ax.set_xlabel("|cis-eQTL Z| (instrument strength)"); ax.set_ylabel("genes"); ax.legend(fontsize=8)
save(fig,"instrument_strength","cis-eQTL instrument strength across immune genes")

# A4 instrument sample size
fig,ax=plt.subplots(figsize=(6.5,4.2))
ax.hist(inst.NrSamples,bins=40,color="#8e44ad")
ax.set_xlabel("eQTLGen NrSamples per instrument"); ax.set_ylabel("genes")
save(fig,"instrument_samplesize","eQTLGen discovery sample size per instrument")

# A5 tests per disease
tc=mr.disease.value_counts().sort_values()
fig,ax=plt.subplots(figsize=(7,4.5))
ax.barh(tc.index,tc.values,color="#16a085")
ax.set_xlabel("MR tests (genes with instrument)")
save(fig,"tests_per_disease","Number of MR tests per disease")

# A6 significant hits per disease
sc=mr[mr.FDR<0.05].disease.value_counts()
allc=pd.Series(0,index=mr.disease.unique());
for k,v in sc.items(): allc[k]=v
allc=allc.sort_values()
fig,ax=plt.subplots(figsize=(7,4.5))
ax.barh(allc.index,allc.values,color="#c0392b")
ax.set_xlabel("FDR<0.05 gene-disease hits")
for i,v in enumerate(allc.values):
    if v: ax.text(v+0.1,i,str(int(v)),va="center",fontsize=8)
save(fig,"hits_per_disease","Significant (FDR<0.05) hits per disease")

# ================================================================
# GROUP B — per-disease MR volcano
# ================================================================
for dis,d in mr.groupby("disease"):
    fig,ax=plt.subplots(figsize=(6.8,5))
    x=np.log(d.OR); y=-np.log10(d.MR_p.clip(lower=1e-300))
    sig=d.FDR<0.05
    ax.scatter(x[~sig],y[~sig],s=10,c="#cccccc",label="ns")
    ax.scatter(x[sig],y[sig],s=28,c="#c0392b",edgecolor="k",lw=.3,label="FDR<0.05")
    for _,r in d[sig].sort_values("MR_p").head(8).iterrows():
        ax.annotate(r.gene_symbol,(np.log(r.OR),-np.log10(max(r.MR_p,1e-300))),
                    fontsize=7,xytext=(3,2),textcoords="offset points")
    ax.axvline(0,ls="--",c="#888",lw=.7)
    ax.set_xlabel("cis-MR log(OR) per SD expr"); ax.set_ylabel("-log10(P)")
    ax.legend(fontsize=8)
    save(fig,f"volcano_{dis[:14].replace(' ','_')}",f"cis-MR volcano \u2014 {dis}")

# ================================================================
# GROUP C — per-disease QQ plot
# ================================================================
for dis,d in mr.groupby("disease"):
    p=np.sort(d.MR_p.clip(lower=1e-300).values)
    exp=-np.log10((np.arange(1,len(p)+1))/(len(p)+1))
    obs=-np.log10(p)
    # genomic inflation lambda
    from scipy.stats import chi2
    z2=chi2.isf(d.MR_p.clip(lower=1e-300),1)
    lam=np.median(z2)/chi2.ppf(0.5,1)
    fig,ax=plt.subplots(figsize=(5.2,5))
    ax.scatter(exp,obs,s=9,c="#1f4e79")
    m=max(exp.max(),obs.max()); ax.plot([0,m],[0,m],ls="--",c="red",lw=1)
    ax.set_xlabel("expected -log10(P)"); ax.set_ylabel("observed -log10(P)")
    ax.text(0.05,0.92,f"\u03bb = {lam:.2f}",transform=ax.transAxes,fontsize=9)
    save(fig,f"qq_{dis[:14].replace(' ','_')}",f"MR P-value QQ \u2014 {dis}")

# ================================================================
# GROUP D — per-disease forest (top hits)
# ================================================================
for dis,d in mr.groupby("disease"):
    top=d.sort_values("MR_p").head(15).sort_values("OR")
    fig,ax=plt.subplots(figsize=(6.8,max(3,0.32*len(top)+1)))
    y=np.arange(len(top))
    colors=["#c0392b" if f<0.05 else "#888" for f in top.FDR]
    ax.errorbar(top.OR,y,xerr=[top.OR-top.OR_l95,top.OR_u95-top.OR],fmt="o",
                ecolor="#bbb",capsize=2,ls="none",mfc="w",mec="w")
    ax.scatter(top.OR,y,c=colors,s=30,zorder=3,edgecolor="k",lw=.3)
    ax.axvline(1,ls="--",c="red",lw=1); ax.set_yticks(y); ax.set_yticklabels(top.gene_symbol,fontsize=7.5)
    ax.set_xscale("log"); ax.set_xlabel("OR per SD cis-expr (95% CI)")
    save(fig,f"forest_{dis[:14].replace(' ','_')}",f"Top cis-MR effects \u2014 {dis}")

# ================================================================
# GROUP E — per-gene cis-pQTL regional (INTERVAL)
# ================================================================
for fp in sorted(glob.glob(os.path.join(CIS,"*.tsv"))):
    base=os.path.basename(fp)[:-4]; gene,acc=base.split("__")
    try:
        df=pd.read_csv(fp,sep="\t")
        df["p_value"]=pd.to_numeric(df["p_value"],errors="coerce")
        df["base_pair_location"]=pd.to_numeric(df["base_pair_location"],errors="coerce")
        df=df.dropna(subset=["p_value","base_pair_location"])
        if len(df)<10: continue
        fig,ax=plt.subplots(figsize=(7.5,4))
        y=-np.log10(df.p_value.clip(lower=1e-300))
        ax.scatter(df.base_pair_location/1e6,y,s=7,c="#5b9bd5")
        lead=df.loc[df.p_value.idxmin()]
        ax.scatter(lead.base_pair_location/1e6,-np.log10(max(lead.p_value,1e-300)),
                   s=60,c="#c0392b",edgecolor="k",zorder=5,label=f"lead {lead['hm_rsid'] if 'hm_rsid' in df else ''}")
        ax.axhline(-np.log10(5e-8),ls="--",c="green",lw=.8)
        ax.set_xlabel(f"chr{int(lead.chromosome)} position (Mb)"); ax.set_ylabel("cis-pQTL -log10(P)")
        ax.legend(fontsize=7)
        save(fig,f"pqtl_regional_{gene}",f"Plasma cis-pQTL regional association \u2014 {gene} ({acc})")
    except Exception as e:
        print("skip",gene,e)

# ================================================================
# GROUP F — discovery vs replication per covered disease
# ================================================================
if rep is not None:
    fnd=rep[rep.rep_found==True].copy()
    dm=mr[["gene_symbol","disease","MR_p"]].rename(columns={"gene_symbol":"gene"})
    fnd=fnd.merge(dm,on=["gene","disease"],how="left")
    for dis,d in fnd.groupby("disease"):
        d=d.dropna(subset=["MR_p","rep_p"])
        if len(d)<1: continue
        fig,ax=plt.subplots(figsize=(5.6,5))
        dx=-np.log10(d.MR_p.clip(lower=1e-300)); dy=-np.log10(pd.to_numeric(d.rep_p,errors="coerce").clip(lower=1e-300))
        col=["#1f4e79" if s else "#c0392b" for s in (d.rep_concordant=="yes")]
        ax.scatter(dx,dy,c=col,s=45,edgecolor="k",lw=.3)
        for _,r in d.iterrows():
            ax.annotate(r.gene,(-np.log10(max(r.MR_p,1e-300)),-np.log10(max(float(r.rep_p),1e-300))),
                        fontsize=7,xytext=(3,2),textcoords="offset points")
        ax.axhline(-np.log10(0.05),ls="--",c="green",lw=.8)
        ax.set_xlabel("discovery (FinnGen) -log10(P)"); ax.set_ylabel("replication -log10(P)")
        save(fig,f"repl_{dis[:14].replace(' ','_')}",f"Discovery vs independent replication \u2014 {dis}")

print(f"\nGENERATED {n} supplementary figures -> {SUP}")
print("total files:",len(os.listdir(SUP)))
