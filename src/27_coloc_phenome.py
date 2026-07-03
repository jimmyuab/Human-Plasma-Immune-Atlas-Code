#!/usr/bin/env python
"""
HDDM Layer 54 - Step 27
Colocalization for the NEW pan-phenome FDR<0.05 cis-MR hits (src/24), using the
identical corrected per-SNP Wakefield ABF coloc.abf as src/11. Non-MHC only.
Extends the cached eQTLGen cis extract with any newly-needed genes (one full-file
scan, then re-cached).

Output: 06_genetic_causality/coloc_phenome_results.tsv
"""
import os, gzip, math, pickle
import pandas as pd, numpy as np
from collections import defaultdict

ROOT = r"I:\Plasma immune atalas"
RAW  = os.path.join(ROOT, "01_data_raw")
GEN  = os.path.join(ROOT, "06_genetic_causality")
SS   = os.path.join(RAW, "FinnGen_GWAS", "sumstats")
EQTL_FULL = os.path.join(RAW, "eQTLGen", "cis-eQTLs-full.txt.gz")
CACHE = os.path.join(GEN, "coloc_eqtl_cis_cache.pkl")

mr = pd.read_csv(os.path.join(GEN, "cis_MR_phenome_results.tsv"), sep="\t")
sig = mr[(mr.FDR < 0.05) & (mr.chr != 6)].copy()
targets = sig[["gene_symbol", "disease_code", "disease", "chr"]].drop_duplicates()
genes = set(targets.gene_symbol)
print("non-MHC pan-phenome coloc targets:", len(targets), "| genes:", len(genes))

# ---- load/extend eQTL cis cache ----
eq = {}
if os.path.exists(CACHE):
    with open(CACHE, "rb") as fh: eq = pickle.load(fh)
    print(f"loaded cache: {len(eq)} genes")
missing = genes - set(eq.keys())
if missing:
    print(f"scanning eQTLGen full for {len(missing)} new genes ...")
    add = defaultdict(dict)
    with gzip.open(EQTL_FULL, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t"); ix = {h:i for i,h in enumerate(hdr)}
        gi, si, zi = ix["GeneSymbol"], ix["SNP"], ix["Zscore"]
        a1i, a2i, ni = ix["AssessedAllele"], ix["OtherAllele"], ix["NrSamples"]
        nrow = 0
        for line in f:
            nrow += 1
            p = line.rstrip("\n").split("\t")
            if p[gi] not in missing: continue
            add[p[gi]][p[si]] = (float(p[zi]), p[a1i].upper(), p[a2i].upper(), float(p[ni]))
    for g,v in add.items(): eq[g] = v
    print(f"  scanned {nrow:,} rows; cache now {len(eq)} genes")
    with open(CACHE, "wb") as fh: pickle.dump(eq, fh)

# ---- coloc ----
COMP = {"A":"T","T":"A","C":"G","G":"C"}
def ambiguous(a,b): return COMP.get(a)==b
W_eqtl = 0.15**2; W_gwas = 0.2**2
def lABF(Z,V,W): r=W/(W+V); return 0.5*(math.log(1-r)+r*Z*Z)
def logsumexp(a): m=np.max(a); return m+math.log(np.sum(np.exp(a-m)))

results = []
for code in targets.disease_code.unique():
    fn = os.path.join(SS, f"finngen_R12_{code}.gz")
    if not os.path.exists(fn): print("  MISSING", fn); continue
    gset = set(targets.loc[targets.disease_code==code, "gene_symbol"])
    wanted = set()
    for g in gset:
        if g in eq: wanted |= set(eq[g].keys())
    out = {}
    with gzip.open(fn, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); jx = {c:i for i,c in enumerate(h)}
        rsi, refi, alti = jx["rsids"], jx["ref"], jx["alt"]
        bi, sei, afi = jx["beta"], jx["sebeta"], jx["af_alt_controls"]
        for line in f:
            if "rs" not in line: continue
            p = line.rstrip("\n").split("\t"); rs = p[rsi]
            if rs not in wanted: continue
            try:
                b=float(p[bi]); se=float(p[sei]); af=float(p[afi])
                if se<=0 or not (0<af<1): continue
                out[rs]=(b, se, p[refi].upper(), p[alti].upper(), af)
            except (ValueError,IndexError): continue
    for g in gset:
        if g not in eq: continue
        dz = targets.loc[(targets.disease_code==code)&(targets.gene_symbol==g),"disease"].iloc[0]
        labf_snp = []
        for rs,(ze,a1,a2,N) in eq[g].items():
            if rs not in out: continue
            b_o, se_o, ref, alt, af = out[rs]
            if ambiguous(a1,a2): continue
            if alt==a1 and ref==a2: sign=1.0
            elif ref==a1 and alt==a2: sign=-1.0
            else: continue
            Zg=b_o/se_o; Vg=se_o*se_o
            denom=2.0*af*(1.0-af)*(N+ze*ze)
            if denom<=0: continue
            se_e=1.0/math.sqrt(denom); Ve=se_e*se_e; Ze=ze*sign
            labf_snp.append((lABF(Ze,Ve,W_eqtl), lABF(Zg,Vg,W_gwas)))
        nsnp=len(labf_snp)
        if nsnp<5:
            results.append(dict(gene=g,disease=dz,disease_code=code,nsnp=nsnp,
                                PP_H0=np.nan,PP_H1=np.nan,PP_H2=np.nan,PP_H3=np.nan,PP_H4=np.nan)); continue
        le=np.array([x[0] for x in labf_snp]); lg=np.array([x[1] for x in labf_snp])
        p1=p2=1e-4; p12=1e-5
        l1=logsumexp(le); l2=logsumexp(lg); l4=logsumexp(le+lg)
        l12=l1+l2; diff=1.0-math.exp(l4-l12)
        lH3=(-np.inf if diff<=0 else math.log(p1)+math.log(p2)+l12+math.log(diff))
        arr=np.array([0.0, math.log(p1)+l1, math.log(p2)+l2, lH3, math.log(p12)+l4])
        pp=np.exp(arr-logsumexp(arr))
        results.append(dict(gene=g,disease=dz,disease_code=code,nsnp=nsnp,
                            PP_H0=pp[0],PP_H1=pp[1],PP_H2=pp[2],PP_H3=pp[3],PP_H4=pp[4]))
    print(f"  {code}: {len(gset)} genes ({len(out)} outcome SNPs)", flush=True)

cd = pd.DataFrame(results).sort_values("PP_H4", ascending=False)
cd.to_csv(os.path.join(GEN, "coloc_phenome_results.tsv"), sep="\t", index=False)
print("\n=== PAN-PHENOME COLOC ===")
print(cd[["gene","disease","nsnp","PP_H4","PP_H3"]].head(25).to_string(index=False))
print("\nPP.H4>=0.8:", int((cd.PP_H4>=0.8).sum()), "| >=0.5:", int((cd.PP_H4>=0.5).sum()))
print("wrote coloc_phenome_results.tsv")
