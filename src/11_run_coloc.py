#!/usr/bin/env python
"""
HDDM Layer 54 - Step 11 (corrected)
Colocalization (coloc.abf / Wakefield 2009 approximate Bayes factors) for the
FDR-significant NON-MHC cis-MR hits: shared causal variant between
  eQTLGen cis-eQTL (exposure)  and  FinnGen disease GWAS (outcome).

Proper per-SNP-variance Wakefield ABF (fixes earlier V=1 dilution bug):
  For each shared cis SNP:
    - FinnGen (outcome): real beta, sebeta  ->  Z_g = beta/sebeta,  V_g = sebeta^2
    - eQTLGen (exposure): only Z + N given.  Reconstruct SE via Zhu et al. 2016
        se_e = 1 / sqrt( 2*af*(1-af)*(N + z^2) )        (af from FinnGen controls)
      so V_e = se_e^2,  Z_e = eQTL Zscore (given directly).
    - lABF = 0.5*( log(1-r) + r*Z^2 ),  r = W/(W+V)   per SNP
        W_e = 0.15^2 (eQTL prior var),  W_g = 0.2^2 (log-OR prior var)
  Posterior for 5 hypotheses (H0..H4), priors p1=p2=1e-4, p12=1e-5
  (Giambartolomei 2014 combinatorics).

eQTL cis-SNP extract is cached to a pickle so the 127M-row full file is
scanned only once.
Outputs 06_genetic_causality/coloc_results.tsv (PP.H0..H4 per locus).
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

# ---- the non-MHC FDR<0.05 hits to colocalize (gene, disease_code) ----
mr = pd.read_csv(os.path.join(GEN, "cis_MR_immune_results.tsv"), sep="\t")
sig = mr[(mr.FDR < 0.05) & (mr.chr != 6)].copy()
targets = sig[["gene_symbol", "disease_code", "disease", "chr"]].drop_duplicates()
genes = set(targets.gene_symbol)
print("non-MHC coloc targets:", len(targets), "| genes:", len(genes))

# ---- 1. pull ALL cis SNPs for these genes from eQTLGen full (cached) ----
# eq: gene -> {rsid: (Zscore, AssessedAllele, OtherAllele, NrSamples)}
if os.path.exists(CACHE):
    with open(CACHE, "rb") as fh:
        eq = pickle.load(fh)
    # if cache doesn't cover all current target genes, rebuild
    if not genes.issubset(set(eq.keys())):
        eq = None
    else:
        print(f"loaded eQTL cis cache: {len(eq)} genes")
else:
    eq = None

if eq is None:
    print("scanning eQTLGen full for cis SNPs of target genes ...")
    eq = defaultdict(dict)
    with gzip.open(EQTL_FULL, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        gi, si, zi = ix["GeneSymbol"], ix["SNP"], ix["Zscore"]
        a1i, a2i, ni = ix["AssessedAllele"], ix["OtherAllele"], ix["NrSamples"]
        n = 0
        for line in f:
            n += 1
            p = line.rstrip("\n").split("\t")
            g = p[gi]
            if g not in genes:
                continue
            eq[g][p[si]] = (float(p[zi]), p[a1i].upper(), p[a2i].upper(), float(p[ni]))
    eq = dict(eq)
    print(f"  scanned {n:,} rows; captured cis SNPs for {len(eq)} genes",
          "| median cis SNPs:", int(np.median([len(v) for v in eq.values()])) if eq else 0)
    with open(CACHE, "wb") as fh:
        pickle.dump(eq, fh)
    print("  cached ->", CACHE)

# ---- 2. coloc per disease ----
COMP = {"A":"T","T":"A","C":"G","G":"C"}
def ambiguous(a,b): return COMP.get(a)==b

W_eqtl = 0.15**2
W_gwas = 0.2**2

def lABF(Z, V, W):
    r = W / (W + V)
    return 0.5 * (math.log(1 - r) + r * Z * Z)

def logsumexp(a):
    m = np.max(a)
    return m + math.log(np.sum(np.exp(a - m)))

results = []
codes = targets.disease_code.unique()
for code in codes:
    fn = os.path.join(SS, f"finngen_R12_{code}.gz")
    if not os.path.exists(fn):
        print("  MISSING", fn); continue
    gset = set(targets.loc[targets.disease_code==code, "gene_symbol"])
    wanted = set()
    for g in gset: wanted |= set(eq[g].keys())
    # pull outcome beta/se/af for wanted rsids
    out = {}   # rsid -> (beta, sebeta, ref, alt, af_ctrl)
    with gzip.open(fn, "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        jx = {c:i for i,c in enumerate(h)}
        rsi, refi, alti = jx["rsids"], jx["ref"], jx["alt"]
        bi, sei, afi = jx["beta"], jx["sebeta"], jx["af_alt_controls"]
        for line in f:
            if "rs" not in line: continue
            p = line.rstrip("\n").split("\t")
            rs = p[rsi]
            if rs not in wanted: continue
            try:
                b=float(p[bi]); se=float(p[sei]); af=float(p[afi])
                if se<=0 or not (0<af<1): continue
                out[rs]=(b, se, p[refi].upper(), p[alti].upper(), af)
            except (ValueError,IndexError): continue
    for g in gset:
        dz = targets.loc[(targets.disease_code==code)&(targets.gene_symbol==g),"disease"].iloc[0]
        labf_snp = []
        for rs,(ze,a1,a2,N) in eq[g].items():
            if rs not in out: continue
            b_o, se_o, ref, alt, af = out[rs]
            if ambiguous(a1,a2): continue
            # orient eQTL assessed allele a1 to FinnGen alt
            if alt==a1 and ref==a2:   sign = 1.0
            elif ref==a1 and alt==a2: sign = -1.0
            else: continue
            # outcome (FinnGen): real per-SNP variance
            Zg = (b_o / se_o)
            Vg = se_o * se_o
            # exposure (eQTL): reconstruct SE via Zhu using FinnGen control AF + eQTL N
            denom = 2.0 * af * (1.0 - af) * (N + ze*ze)
            if denom <= 0: continue
            se_e = 1.0 / math.sqrt(denom)
            Ve = se_e * se_e
            Ze = ze * sign   # oriented to FinnGen alt allele (sign only affects Z^2 -> no effect, kept for clarity)
            le = lABF(Ze, Ve, W_eqtl)
            lg = lABF(Zg, Vg, W_gwas)
            labf_snp.append((le, lg))
        nsnp = len(labf_snp)
        if nsnp < 5:
            results.append(dict(gene=g, disease=dz, disease_code=code, nsnp=nsnp,
                                PP_H0=np.nan,PP_H1=np.nan,PP_H2=np.nan,PP_H3=np.nan,PP_H4=np.nan))
            continue
        le_arr = np.array([x[0] for x in labf_snp]); lg_arr = np.array([x[1] for x in labf_snp])
        p1=1e-4; p2=1e-4; p12=1e-5
        l1 = logsumexp(le_arr)
        l2 = logsumexp(lg_arr)
        l4 = logsumexp(le_arr+lg_arr)
        lH0 = 0.0
        lH1 = math.log(p1) + l1
        lH2 = math.log(p2) + l2
        # H3: independent causal variants = p1*p2*[ (sum_i ABFe_i)(sum_j ABFg_j) - sum_k ABFe_k*ABFg_k ]
        # log( exp(l1+l2) - exp(l4) ) computed stably (l4 <= l1+l2)
        l12 = l1 + l2
        diff = 1.0 - math.exp(l4 - l12)
        if diff <= 0:
            lH3 = -np.inf
        else:
            lH3 = math.log(p1)+math.log(p2) + l12 + math.log(diff)
        lH4 = math.log(p12) + l4
        arr = np.array([lH0,lH1,lH2,lH3,lH4])
        denom = logsumexp(arr)
        pp = np.exp(arr-denom)
        results.append(dict(gene=g, disease=dz, disease_code=code, nsnp=nsnp,
                            PP_H0=pp[0],PP_H1=pp[1],PP_H2=pp[2],PP_H3=pp[3],PP_H4=pp[4]))
    print(f"  {code}: coloc done for {len(gset)} genes ({len(out)} outcome SNPs matched)")

cd = pd.DataFrame(results).sort_values("PP_H4", ascending=False)
cd.to_csv(os.path.join(GEN, "coloc_results.tsv"), sep="\t", index=False)
print("\n=== COLOC RESULTS (top) ===")
print(cd[["gene","disease","nsnp","PP_H4","PP_H3"]].head(20).to_string(index=False))
print("\nPP.H4>=0.8 (strong coloc):", int((cd.PP_H4>=0.8).sum()),
      "| >=0.5:", int((cd.PP_H4>=0.5).sum()))
print("wrote", os.path.join(GEN, "coloc_results.tsv"))
