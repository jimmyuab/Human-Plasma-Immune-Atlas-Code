#!/usr/bin/env python
"""
HDDM Layer 54 - Step 24
PAN-PHENOME immune cis-MR. Extends the validated cis-eQTL (eQTLGen) -> disease
(FinnGen R12) Wald-ratio MR from the 13 autoimmune diseases to the full expanded
phenome (cardiovascular / metabolic / renal / neuro / aging). Identical, already
validated Zhu (2016) Z->beta/se conversion and allele harmonisation as src/07 —
only the disease set (and a disease_category tag) is broadened.

Output: 06_genetic_causality/cis_MR_phenome_results.tsv
"""
import os, gzip, math
import pandas as pd

ROOT = r"I:\Plasma immune atalas"
RAW  = os.path.join(ROOT, "01_data_raw")
GEN  = os.path.join(ROOT, "06_genetic_causality")
SS   = os.path.join(RAW, "FinnGen_GWAS", "sumstats")

inst = pd.read_csv(os.path.join(GEN, "immune_cis_eqtl_instruments.tsv"), sep="\t")
inst["SNP"] = inst["SNP"].astype(str)
snp2genes = {}
for _, r in inst.iterrows():
    snp2genes.setdefault(r["SNP"], []).append(r.to_dict())
target_snps = set(inst["SNP"])
print("gene instruments:", len(inst), "| unique SNPs:", len(target_snps))

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def ambiguous(a, b): return COMP.get(a) == b

# code -> (readable name, category). Autoimmune block kept identical to src/07.
DISEASE = {
 "AUTOIMMUNE_HYPERTHYROIDISM": ("Autoimmune hyperthyroidism","Autoimmune"),
 "CHRONNAS": ("Crohn's disease","Autoimmune"), "D3_SARCOIDOSIS": ("Sarcoidosis","Autoimmune"),
 "E4_THYROIDITAUTOIM": ("Autoimmune thyroiditis","Autoimmune"),
 "G6_GUILBAR": ("Guillain-Barre","Autoimmune"), "G6_MS": ("Multiple sclerosis","Autoimmune"),
 "K11_COELIAC": ("Coeliac disease","Autoimmune"), "L12_PSORIASIS": ("Psoriasis","Autoimmune"),
 "L12_VITILIGO": ("Vitiligo","Autoimmune"), "M13_ANKYLOSPON": ("Ankylosing spondylitis","Autoimmune"),
 "M13_RHEUMA": ("Rheumatoid arthritis","Autoimmune"), "M13_SJOGREN": ("Sjogren syndrome","Autoimmune"),
 "SLE_FG": ("Systemic lupus erythematosus","Autoimmune"),
 # expanded non-immune phenome
 "I9_HYPTENS": ("Hypertension","Cardiovascular"), "I9_AF": ("Atrial fibrillation","Cardiovascular"),
 "I9_CHD": ("Coronary heart disease","Cardiovascular"), "I9_HEARTFAIL": ("Heart failure","Cardiovascular"),
 "I9_STR": ("Stroke","Cardiovascular"), "I9_VTE": ("Venous thromboembolism","Cardiovascular"),
 "T2D": ("Type 2 diabetes","Metabolic"), "T1D": ("Type 1 diabetes","Metabolic"),
 "E4_OBESITY": ("Obesity","Metabolic"),
 "N14_CHRONKIDNEYDIS": ("Chronic kidney disease","Renal"),
 "F5_DEMENTIA": ("Dementia","Neuro/Aging"), "G6_AD_WIDE": ("Alzheimer's disease","Neuro/Aging"),
 "G6_EPLEPSY": ("Epilepsy","Neuro/Aging"), "M13_OSTEOPOROSIS": ("Osteoporosis","Neuro/Aging"),
 "H7_GLAUCOMA": ("Glaucoma","Neuro/Aging"),
}

def zhu_beta_se(z, n, p):
    denom = 2 * p * (1 - p) * (n + z * z)
    if denom <= 0: return None, None
    s = math.sqrt(denom); return z / s, 1.0 / s

rows = []
files = sorted([f for f in os.listdir(SS) if f.endswith(".gz")])
for fi, fn in enumerate(files, 1):
    code = fn.replace("finngen_R12_", "").replace(".gz", "")
    dname, dcat = DISEASE.get(code, (code, "Other"))
    path = os.path.join(SS, fn); found = 0
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        rsi = ix["rsids"]; refi = ix["ref"]; alti = ix["alt"]
        betai = ix["beta"]; sei = ix["sebeta"]; pvi = ix["pval"]
        afci = ix.get("af_alt_controls", ix.get("af_alt"))
        for line in f:
            if "rs" not in line: continue
            p = line.rstrip("\n").split("\t")
            rs = p[rsi]
            if rs not in target_snps: continue
            ref = p[refi].upper(); alt = p[alti].upper()
            try:
                b_out = float(p[betai]); se_out = float(p[sei]); pv_out = float(p[pvi]); af = float(p[afci])
            except (ValueError, IndexError): continue
            for inst_row in snp2genes[rs]:
                ea = inst_row["effect_allele"].upper(); oa = inst_row["other_allele"].upper()
                if ambiguous(ea, oa): continue
                if alt == ea and ref == oa: bo = b_out; pe = af
                elif ref == ea and alt == oa: bo = -b_out; pe = 1 - af
                else: continue
                if not (0 < pe < 1): continue
                z = inst_row["Zscore"]; n = inst_row["NrSamples"]
                b_exp, se_exp = zhu_beta_se(z, n, pe)
                if b_exp is None or b_exp == 0: continue
                theta = bo / b_exp; se_theta = se_out / abs(b_exp)
                zt = theta / se_theta; pval = math.erfc(abs(zt) / math.sqrt(2))
                rows.append({
                    "gene_symbol": inst_row["gene_symbol"], "disease_code": code,
                    "disease": dname, "disease_category": dcat, "SNP": rs, "chr": inst_row["chr"],
                    "eqtl_Z": z, "eqtl_N": n, "eaf": round(pe, 4),
                    "b_exp": b_exp, "b_out_oriented": bo, "se_out": se_out,
                    "MR_beta": theta, "MR_se": se_theta, "MR_z": zt, "MR_p": pval,
                    "OR": math.exp(theta), "OR_l95": math.exp(theta - 1.96 * se_theta),
                    "OR_u95": math.exp(theta + 1.96 * se_theta),
                })
                found += 1
    print(f"[{fi}/{len(files)}] {dname:28s} [{dcat:13s}] matched: {found}", flush=True)

res = pd.DataFrame(rows).sort_values("MR_p").reset_index(drop=True)
m = len(res)
res["FDR"] = (res["MR_p"] * m / (res.index + 1)).clip(upper=1.0)
res["FDR"] = res["FDR"][::-1].cummin()[::-1]
ann = pd.read_csv(os.path.join(ROOT, "02_data_processed", "plasma_immune_protein_annotation.tsv"), sep="\t")
res = res.merge(ann[["gene_symbol", "immune_class", "immune_source_cells"]], on="gene_symbol", how="left")
out = os.path.join(GEN, "cis_MR_phenome_results.tsv")
res.to_csv(out, sep="\t", index=False)

print("\n=== PAN-PHENOME cis-MR ===")
print("diseases:", res.disease.nunique(), "| tests:", m,
      "| FDR<0.05:", int((res.FDR < 0.05).sum()), "| nominal p<0.05:", int((res.MR_p < 0.05).sum()))
print("\nFDR<0.05 by category:")
print(res[res.FDR < 0.05].disease_category.value_counts().to_string())
print("\nTop 20 by FDR:")
cols = ["gene_symbol", "disease", "disease_category", "OR", "MR_p", "FDR"]
print(res[res.FDR < 0.05][cols].head(20).to_string(index=False))
print("\nwrote", out)
