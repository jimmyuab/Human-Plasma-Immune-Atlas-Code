#!/usr/bin/env python
"""
HDDM Layer 54 - Step 16
Protein-level cis-pQTL Mendelian randomization + colocalization, using INTERVAL
(Sun 2018) plasma pQTLs as the exposure and FinnGen R12 as the outcome.  This
UPGRADES the transcript-level (eQTL) analysis to protein level - the key
reviewer ask - with real beta+SE on both sides (no Zhu reconstruction needed).

For each aptamer cis file:
  instrument = strongest cis-pQTL (min p) with rsid, |beta|, se, EAF
  For each of 13 FinnGen diseases:
     harmonise instrument -> Wald ratio MR (beta_out/beta_exp), OR, CI
  Coloc (coloc.abf, proper per-SNP variance) over shared cis rsids.
Outputs:
  06_genetic_causality/pqtl_MR_results.tsv
  06_genetic_causality/pqtl_coloc_results.tsv
"""
import os, glob, gzip, math
import pandas as pd, numpy as np
from statsmodels.stats.multitest import multipletests

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality")
CIS  = os.path.join(ROOT, "01_data_raw", "INTERVAL_pQTL", "cis")
SS   = os.path.join(ROOT, "01_data_raw", "FinnGen_GWAS", "sumstats")

DIS = {
 "M13_ANKYLOSPON":"Ankylosing spondylitis","AUTOIMMUNE_HYPERTHYROIDISM":"Autoimmune hyperthyroidism",
 "E4_THYROIDITAUTOIM":"Autoimmune thyroiditis","K11_COELIAC":"Coeliac disease","CHRONNAS":"Crohn's disease",
 "G6_GUILBAR":"Guillain-Barre","G6_MS":"Multiple sclerosis","L12_PSORIASIS":"Psoriasis",
 "M13_RHEUMA":"Rheumatoid arthritis","D3_SARCOIDOSIS":"Sarcoidosis","M13_SJOGREN":"Sjogren syndrome",
 "SLE_FG":"Systemic lupus erythematosus","L12_VITILIGO":"Vitiligo",
}
COMP = {"A":"T","T":"A","C":"G","G":"C"}
def ambiguous(a,b): return COMP.get(a)==b

# ---- load cis pQTL per aptamer ----
def load_cis(fp):
    df = pd.read_csv(fp, sep="\t", dtype=str)
    df = df.rename(columns={"hm_rsid":"rsid"})
    for c in ["beta","standard_error","p_value","effect_allele_frequency","base_pair_location"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["rsid","beta","standard_error","p_value","effect_allele","other_allele"])
    df["effect_allele"] = df["effect_allele"].str.upper()
    df["other_allele"]  = df["other_allele"].str.upper()
    df = df[df.rsid.str.startswith("rs")]
    return df

cis_files = sorted(glob.glob(os.path.join(CIS, "*.tsv")))
apt = {}   # (gene,acc) -> df
for fp in cis_files:
    base = os.path.basename(fp)[:-4]
    gene, acc = base.split("__")
    apt[(gene,acc)] = load_cis(fp)
print("loaded", len(apt), "aptamer cis files")

# collect all rsids we may need from FinnGen
need = set()
for df in apt.values(): need |= set(df.rsid)

# ---- pull FinnGen rows for needed rsids per disease ----
def load_fg(code):
    fn = os.path.join(SS, f"finngen_R12_{code}.gz")
    out = {}
    with gzip.open(fn, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); jx={c:i for i,c in enumerate(h)}
        rsi,refi,alti = jx["rsids"],jx["ref"],jx["alt"]
        bi,sei,afi = jx["beta"],jx["sebeta"],jx["af_alt_controls"]
        for line in f:
            if "rs" not in line: continue
            p=line.rstrip("\n").split("\t"); rs=p[rsi]
            if rs not in need: continue
            try:
                b=float(p[bi]); se=float(p[sei]); af=float(p[afi])
                if se<=0: continue
                out[rs]=(b,se,p[refi].upper(),p[alti].upper(),af)
            except (ValueError,IndexError): continue
    return out

# ---- coloc helpers (per-SNP variance Wakefield ABF) ----
W_pqtl=0.15**2; W_gwas=0.2**2
def lABF(Z,V,W):
    r=W/(W+V); return 0.5*(math.log(1-r)+r*Z*Z)
def lse(a):
    m=np.max(a); return m+math.log(np.sum(np.exp(a-m)))

mr_rows=[]; coloc_rows=[]
for code,dname in DIS.items():
    fg = load_fg(code)
    for (gene,acc),df in apt.items():
        # ---- instrument: strongest cis-pQTL present in FinnGen ----
        cand = df[df.rsid.isin(fg)].sort_values("p_value")
        if len(cand)==0: continue
        inst=None
        for _,r in cand.iterrows():
            if ambiguous(r.effect_allele,r.other_allele): continue
            inst=r; break
        if inst is not None:
            b_o,se_o,ref,alt,af = fg[inst.rsid]
            # orient outcome to pQTL effect allele
            if alt==inst.effect_allele and ref==inst.other_allele: sign=1.0
            elif ref==inst.effect_allele and alt==inst.other_allele: sign=-1.0
            else: sign=None
            if sign is not None:
                b_exp=inst.beta; b_out=b_o*sign
                mr_beta=b_out/b_exp; mr_se=abs(se_o/b_exp)
                z=mr_beta/mr_se; p=2*(1-0.5*(1+math.erf(abs(z)/math.sqrt(2))))
                mr_rows.append(dict(gene=gene,acc=acc,disease=dname,disease_code=code,
                    rsid=inst.rsid,pqtl_p=inst.p_value,b_exp=b_exp,b_out=b_out,se_out=se_o,
                    MR_beta=mr_beta,MR_se=mr_se,MR_p=p,OR=math.exp(mr_beta),
                    OR_l95=math.exp(mr_beta-1.96*mr_se),OR_u95=math.exp(mr_beta+1.96*mr_se)))
        # ---- coloc over shared cis SNPs ----
        labf=[]
        for _,r in df.iterrows():
            if r.rsid not in fg: continue
            if ambiguous(r.effect_allele,r.other_allele): continue
            bo,seo,ref,alt,af = fg[r.rsid]
            if not ((alt==r.effect_allele and ref==r.other_allele) or (ref==r.effect_allele and alt==r.other_allele)):
                continue
            Ze=r.beta/r.standard_error; Ve=r.standard_error**2
            Zg=bo/seo; Vg=seo**2
            labf.append((lABF(Ze,Ve,W_pqtl), lABF(Zg,Vg,W_gwas)))
        if len(labf)>=5:
            le=np.array([x[0] for x in labf]); lg=np.array([x[1] for x in labf])
            p1=p2=1e-4; p12=1e-5
            l1=lse(le); l2=lse(lg); l4=lse(le+lg); l12=l1+l2
            diff=1-math.exp(l4-l12)
            lH=[0.0, math.log(p1)+l1, math.log(p2)+l2,
                (math.log(p1)+math.log(p2)+l12+math.log(diff)) if diff>0 else -np.inf,
                math.log(p12)+l4]
            arr=np.array(lH); pp=np.exp(arr-lse(arr))
            coloc_rows.append(dict(gene=gene,acc=acc,disease=dname,disease_code=code,nsnp=len(labf),
                PP_H0=pp[0],PP_H1=pp[1],PP_H2=pp[2],PP_H3=pp[3],PP_H4=pp[4]))
    print("  done", code)

mrd=pd.DataFrame(mr_rows)
if len(mrd):
    mrd["FDR"]=multipletests(mrd.MR_p,method="fdr_bh")[1]
    mrd=mrd.sort_values("MR_p")
mrd.to_csv(os.path.join(GEN,"pqtl_MR_results.tsv"),sep="\t",index=False)
cld=pd.DataFrame(coloc_rows).sort_values("PP_H4",ascending=False)
cld.to_csv(os.path.join(GEN,"pqtl_coloc_results.tsv"),sep="\t",index=False)

print("\n=== pQTL-MR: FDR<0.05 ===")
print(mrd[mrd.FDR<0.05][["gene","disease","OR","OR_l95","OR_u95","MR_p","FDR"]].to_string(index=False))
print("\n=== pQTL coloc PP.H4>=0.8 ===")
print(cld[cld.PP_H4>=0.8][["gene","disease","nsnp","PP_H4"]].to_string(index=False))
print(f"\npQTL-MR tests={len(mrd)} sig={int((mrd.FDR<0.05).sum())} | coloc loci={len(cld)} strong={int((cld.PP_H4>=0.8).sum())}")
