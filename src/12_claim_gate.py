#!/usr/bin/env python
"""
HDDM Layer 54 - Step 12
Claim-strength gate: assign every MR hit an evidence tier + the exact language
that may be used for it, following the rubric:
  cis-eQTL only            -> "transcript-level genetic proxy"
  cis-MR FDR<0.05          -> "genetically-supported nomination"
  + colocalization PP.H4>=0.8 -> "prioritized causal target (transcript-level)"
  MHC-region               -> forced down to "nomination (MHC LD caution)"
Output: 06_genetic_causality/evidence_tiered_targets.tsv
"""
import os
import pandas as pd, numpy as np

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality")

mr = pd.read_csv(os.path.join(GEN, "cis_MR_immune_results.tsv"), sep="\t")
cd_path = os.path.join(GEN, "coloc_results.tsv")
cd = pd.read_csv(cd_path, sep="\t") if os.path.exists(cd_path) else pd.DataFrame(
    columns=["gene","disease","PP_H4"])
cd = cd.rename(columns={"gene": "gene_symbol"})

m = mr.merge(cd[["gene_symbol","disease","PP_H4","nsnp"]], on=["gene_symbol","disease"], how="left")
m["MHC"] = (m["chr"] == 6)

def tier(r):
    if r.FDR >= 0.05:
        return 1, "transcript-level genetic proxy", "associated with (nominal)"
    # FDR-significant
    if r.MHC:
        return 2, "genetically-supported nomination (MHC LD caution)", \
               "genetically supported, MHC \u2014 requires colocalization"
    if pd.notna(r.PP_H4) and r.PP_H4 >= 0.8:
        return 4, "prioritized causal target (transcript-level)", \
               "colocalized cis-eQTL/GWAS signal (PP.H4\u22650.8)"
    if pd.notna(r.PP_H4) and r.PP_H4 >= 0.5:
        return 3, "genetically-supported nomination (suggestive coloc)", \
               "genetically supported; suggestive colocalization"
    return 3, "genetically-supported nomination", \
           "genetically supported; colocalization pending/absent"

tiers = m.apply(tier, axis=1, result_type="expand")
m["evidence_tier"] = tiers[0]
m["claim_label"] = tiers[1]
m["allowed_language"] = tiers[2]

sig = m[m.FDR < 0.05].sort_values(["evidence_tier","FDR"], ascending=[False, True])
cols = ["gene_symbol","disease","immune_class","chr","OR","OR_l95","OR_u95",
        "MR_p","FDR","PP_H4","nsnp","evidence_tier","claim_label","allowed_language"]
sig[cols].to_csv(os.path.join(GEN, "evidence_tiered_targets.tsv"), sep="\t", index=False)

print("=== EVIDENCE-TIERED SIGNIFICANT TARGETS ===")
print(sig[["gene_symbol","disease","OR","FDR","PP_H4","evidence_tier","claim_label"]].to_string(index=False))
print("\nTier counts (FDR<0.05):")
print(sig["evidence_tier"].value_counts().sort_index().to_string())
print("\nPrioritized causal (tier 4):",
      list(sig.loc[sig.evidence_tier==4, "gene_symbol"] + "\u2192" + sig.loc[sig.evidence_tier==4, "disease"]))
