#!/usr/bin/env python
"""
HDDM Layer 54 - Step 18
Independent replication of significant cis-MR hits in OpenGWAS (MRC-IEU)
disease GWAS that are INDEPENDENT of the FinnGen discovery set.

For each hit: take the discovery eQTL instrument (SNP + effect/other allele +
expression effect b_exp), query the SNP's disease association in the matched
independent GWAS (POST /api/associations {"variant":[...],"id":[gid]}),
harmonise to the eQTL expression-increasing allele, and compute the
replication MR direction. A hit REPLICATES if the independent disease effect
is directionally concordant with discovery (and, stronger, nominally p<0.05).

Token from .secrets/opengwas_token.txt (not committed).
Output: 06_genetic_causality/opengwas_replication.tsv
"""
import os, json, time, urllib.request
import pandas as pd, numpy as np

ROOT=r"I:\Plasma immune atalas"
GEN=os.path.join(ROOT,"06_genetic_causality")
TOK=open(os.path.join(ROOT,".secrets","opengwas_token.txt")).read().strip()
API="https://api.opengwas.io/api"

REP={
 "Multiple sclerosis":       "ieu-b-18",
 "Rheumatoid arthritis":     "ieu-a-833",
 "Ankylosing spondylitis":   "ebi-a-GCST005529",
 "Coeliac disease":          "ieu-a-1058",
 "Psoriasis":                "ebi-a-GCST90019017",
 "Sarcoidosis":              "ebi-a-GCST005538",
 "Crohn's disease":          "ieu-a-12",
 "Autoimmune hyperthyroidism":"ebi-a-GCST90018860",
}
REP_CITE={"ieu-b-18":"Patsopoulos/IMSGC 2019","ieu-a-833":"Okada 2014",
 "ebi-a-GCST005529":"Cortes/IGAS 2013","ieu-a-1058":"Trynka 2011",
 "ebi-a-GCST90019017":"Stuart 2021","ebi-a-GCST005538":"Fischer 2015",
 "ieu-a-12":"Liu/IIBDGC 2015","ebi-a-GCST90018860":"Sakaue 2021"}
COMP={"A":"T","T":"A","C":"G","G":"C"}
def ambiguous(a,b): return COMP.get(a)==b

def post(path,payload,tries=4):
    for i in range(tries):
        try:
            req=urllib.request.Request(API+path,data=json.dumps(payload).encode(),
                headers={"Authorization":f"Bearer {TOK}","Content-Type":"application/json"})
            with urllib.request.urlopen(req,timeout=90) as r:
                return json.load(r)
        except Exception as e:
            print("   retry",i,e); time.sleep(4)
    return None

fin=pd.read_csv(os.path.join(GEN,"FINAL_evidence_tiers.tsv"),sep="\t")
mr=pd.read_csv(os.path.join(GEN,"cis_MR_immune_results.tsv"),sep="\t")
inst=pd.read_csv(os.path.join(GEN,"immune_cis_eqtl_instruments.tsv"),sep="\t")
disc=mr[["gene_symbol","disease","SNP","b_exp"]].drop_duplicates().merge(
     inst[["gene_symbol","SNP","effect_allele","other_allele"]].drop_duplicates(),
     on=["gene_symbol","SNP"],how="left")
hit=fin.merge(disc,on=["gene_symbol","disease"],how="left")

rows=[]
for dis,gid in REP.items():
    sub=hit[hit.disease==dis]
    if len(sub)==0: continue
    rsids=sorted(set(sub.SNP.dropna()))
    assoc=post("/associations",{"variant":rsids,"id":[gid]})
    amap={}
    if assoc:
        for a in assoc: amap[a.get("rsid")]=a
    for _,r in sub.iterrows():
        rec=dict(gene=r.gene_symbol,disease=dis,rep_gwas=gid,rep_cite=REP_CITE[gid],
                 snp=r.SNP,disc_tier=int(r.final_tier),disc_OR=r.OR)
        ar=amap.get(r.SNP)
        if ar is None or pd.isna(r.effect_allele):
            rec.update(rep_found=False,rep_beta=np.nan,rep_p=np.nan,
                       rep_dir_OR=np.nan,rep_concordant="",rep_sig=False)
            rows.append(rec); continue
        ea=str(ar.get("ea","")).upper(); nea=str(ar.get("nea","")).upper()
        try: rb=float(ar.get("beta")); rp=float(ar.get("p"))
        except (TypeError,ValueError):
            rec.update(rep_found=True,rep_beta=np.nan,rep_p=np.nan,rep_dir_OR=np.nan,
                       rep_concordant="",rep_sig=False); rows.append(rec); continue
        e_all=str(r.effect_allele).upper(); o_all=str(r.other_allele).upper()
        # orient replication disease beta to the eQTL effect allele
        if ambiguous(e_all,o_all):
            sign=None
        elif ea==e_all and nea==o_all: sign=1.0
        elif ea==o_all and nea==e_all: sign=-1.0
        else: sign=None
        if sign is None:
            rec.update(rep_found=True,rep_beta=rb,rep_p=rp,rep_dir_OR=np.nan,
                       rep_concordant="allele-mismatch",rep_sig=False); rows.append(rec); continue
        rep_beta_eff=rb*sign                          # disease log-OR per eQTL effect allele
        # discovery direction: OR per +1 SD expr. concordance:
        # expr-increasing direction risk = sign(disc_OR-1). replication same allele MR sign:
        rep_mr_sign=np.sign(rep_beta_eff / r.b_exp) if pd.notna(r.b_exp) and r.b_exp!=0 else np.nan
        disc_sign=np.sign(r.OR-1)
        conc = "yes" if (rep_mr_sign==disc_sign and rep_mr_sign!=0) else "NO"
        rec.update(rep_found=True,rep_beta=rb,rep_p=rp,
                   rep_dir_OR=float(np.exp(rep_beta_eff/abs(r.b_exp))) if r.b_exp else np.nan,
                   rep_concordant=conc,rep_sig=bool(rp<0.05))
        rows.append(rec)
    print(f"  {dis} <- {gid} ({REP_CITE[gid]}): {sum(1 for x in rsids if x in amap)}/{len(rsids)} SNPs")
    time.sleep(0.5)

rep=pd.DataFrame(rows)
rep.to_csv(os.path.join(GEN,"opengwas_replication.tsv"),sep="\t",index=False)

fnd=rep[rep.rep_found==True]
print("\n=== REPLICATION (independent, non-FinnGen GWAS) ===")
print(fnd[["gene","disease","rep_cite","snp","rep_p","rep_concordant","rep_sig","disc_tier"]].to_string(index=False))
print(f"\nSNP found in independent GWAS: {len(fnd)}/{len(rep)}")
print(f"directionally concordant: {(fnd.rep_concordant=='yes').sum()}/{(fnd.rep_concordant.isin(['yes','NO'])).sum()}")
print(f"replicated (concordant AND p<0.05): {((fnd.rep_concordant=='yes')&(fnd.rep_sig)).sum()}")
print("REPLICATED hits:",
      list((fnd[(fnd.rep_concordant=='yes')&(fnd.rep_sig)].gene+'\u2192'+fnd[(fnd.rep_concordant=='yes')&(fnd.rep_sig)].disease)))
print("\nwrote opengwas_replication.tsv")
