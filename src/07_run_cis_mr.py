#!/usr/bin/env python
"""
HDDM Layer 54 - Step 7
Real cis-MR: immune-gene cis-eQTL (eQTLGen) -> immune-disease risk (FinnGen R12).
Single strongest cis instrument per gene => Wald-ratio MR.

eQTLGen gives Z + N (no beta/AF). Convert with Zhu et al. (2016):
    b   = z / sqrt(2*p*(1-p)*(n + z^2))
    se  = 1 / sqrt(2*p*(1-p)*(n + z^2))
using effect-allele freq p taken from FinnGen controls for the matched SNP.

Wald ratio:  theta = b_out / b_exp ;  se_theta = se_out / |b_exp|
Output: 06_genetic_causality/cis_MR_immune_results.tsv
"""
import os, gzip, math
import pandas as pd
import numpy as np

ROOT = r"I:\Plasma immune atalas"
RAW  = os.path.join(ROOT, "01_data_raw")
GEN  = os.path.join(ROOT, "06_genetic_causality")
SS   = os.path.join(RAW, "FinnGen_GWAS", "sumstats")

inst = pd.read_csv(os.path.join(GEN, "immune_cis_eqtl_instruments.tsv"), sep="\t")
inst["SNP"] = inst["SNP"].astype(str)
snp2genes = {}   # one SNP can instrument several nearby immune genes
for _, r in inst.iterrows():
    snp2genes.setdefault(r["SNP"], []).append(r.to_dict())
target_snps = set(inst["SNP"])
print("gene instruments:", len(inst), "| unique SNPs:", len(target_snps))

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def ambiguous(a, b):
    return COMP.get(a) == b  # palindromic A/T or C/G

DISEASE_NAME = {
 "AUTOIMMUNE_HYPERTHYROIDISM": "Autoimmune hyperthyroidism",
 "CHRONNAS": "Crohn's disease", "D3_SARCOIDOSIS": "Sarcoidosis",
 "E4_THYROIDITAUTOIM": "Autoimmune thyroiditis", "G6_GUILBAR": "Guillain-Barre",
 "G6_MS": "Multiple sclerosis", "K11_COELIAC": "Coeliac disease",
 "L12_PSORIASIS": "Psoriasis", "L12_VITILIGO": "Vitiligo",
 "M13_ANKYLOSPON": "Ankylosing spondylitis", "M13_RHEUMA": "Rheumatoid arthritis",
 "M13_SJOGREN": "Sjogren syndrome", "SLE_FG": "Systemic lupus erythematosus",
}

def zhu_beta_se(z, n, p):
    denom = 2 * p * (1 - p) * (n + z * z)
    if denom <= 0:
        return None, None
    s = math.sqrt(denom)
    return z / s, 1.0 / s

rows = []
files = sorted([f for f in os.listdir(SS) if f.endswith(".gz")])
for fi, fn in enumerate(files, 1):
    code = fn.replace("finngen_R12_", "").replace(".gz", "")
    dname = DISEASE_NAME.get(code, code)
    path = os.path.join(SS, fn)
    found = 0
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        rsi = ix["rsids"]; refi = ix["ref"]; alti = ix["alt"]
        betai = ix["beta"]; sei = ix["sebeta"]; pvi = ix["pval"]
        afci = ix.get("af_alt_controls", ix.get("af_alt"))
        for line in f:
            # cheap pre-filter: rsid column appears; only split if 'rs' present
            if "rs" not in line:
                continue
            p = line.rstrip("\n").split("\t")
            rs = p[rsi]
            if rs not in target_snps:
                continue
            ref = p[refi].upper(); alt = p[alti].upper()  # FinnGen beta is per-alt
            try:
                b_out = float(p[betai]); se_out = float(p[sei]); pv_out = float(p[pvi])
                af = float(p[afci])
            except (ValueError, IndexError):
                continue
            for inst_row in snp2genes[rs]:
                ea = inst_row["effect_allele"].upper()   # eQTLGen assessed allele
                oa = inst_row["other_allele"].upper()
                if ambiguous(ea, oa):
                    continue
                # orient FinnGen beta to eQTLGen effect allele
                if alt == ea and ref == oa:
                    bo = b_out; pe = af                 # af is freq of alt = ea
                elif ref == ea and alt == oa:
                    bo = -b_out; pe = 1 - af             # flip
                else:
                    continue  # allele mismatch (multiallelic/strand)
                if not (0 < pe < 1):
                    continue
                z = inst_row["Zscore"]; n = inst_row["NrSamples"]
                b_exp, se_exp = zhu_beta_se(z, n, pe)
                if b_exp is None or b_exp == 0:
                    continue
                theta = bo / b_exp
                se_theta = se_out / abs(b_exp)
                zt = theta / se_theta
                pval = math.erfc(abs(zt) / math.sqrt(2))
                rows.append({
                    "gene_symbol": inst_row["gene_symbol"], "disease_code": code,
                    "disease": dname, "SNP": rs, "chr": inst_row["chr"],
                    "eqtl_Z": z, "eqtl_N": n, "eaf": round(pe, 4),
                    "b_exp": b_exp, "b_out_oriented": bo, "se_out": se_out,
                    "MR_beta": theta, "MR_se": se_theta, "MR_z": zt, "MR_p": pval,
                    "OR": math.exp(theta), "OR_l95": math.exp(theta - 1.96 * se_theta),
                    "OR_u95": math.exp(theta + 1.96 * se_theta),
                })
                found += 1
    print(f"[{fi}/{len(files)}] {dname:30s} instruments matched: {found}")

res = pd.DataFrame(rows)
# Benjamini-Hochberg FDR across all gene x disease tests
res = res.sort_values("MR_p").reset_index(drop=True)
m = len(res)
res["FDR"] = (res["MR_p"] * m / (res.index + 1)).clip(upper=1.0)
res["FDR"] = res["FDR"][::-1].cummin()[::-1]
res = res.merge(
    pd.read_csv(os.path.join(ROOT, "02_data_processed", "plasma_immune_protein_annotation.tsv"),
                sep="\t")[["gene_symbol", "immune_class", "immune_source_cells"]],
    on="gene_symbol", how="left")
out = os.path.join(GEN, "cis_MR_immune_results.tsv")
res.to_csv(out, sep="\t", index=False)

print("\n=== cis-MR RESULTS ===")
print("total gene x disease tests:", m)
print("FDR<0.05:", int((res.FDR < 0.05).sum()), "| nominal p<0.05:", int((res.MR_p < 0.05).sum()))
print("\nTop 15 by FDR:")
cols = ["gene_symbol", "disease", "immune_class", "OR", "OR_l95", "OR_u95", "MR_p", "FDR"]
print(res[cols].head(15).to_string(index=False))
print("\nwrote", out)
