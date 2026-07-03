#!/usr/bin/env python
"""
HDDM Layer 54 - Step 28
Integrated NOVELTY PRIORITY ENGINE across the full phenome. Combines, per
gene-disease pair, every real-data evidence layer into one auditable score:

  NoveltyPriority = causal_evidence + coloc_evidence + pleiotropy
                    + druggability + cellsource_specificity
                    - known_drug_penalty - MHC_LD_penalty

All components are computed from local files (MR, coloc, HPA annotation). Known
positive controls are DOWN-weighted (they validate the pipeline but are not
novel); MHC/LD-confounded loci are penalised. Output ranks novel, high-evidence
immune targets across autoimmune + cardiometabolic + aging disease.

Outputs:
  06_genetic_causality/novelty_engine_ranked.tsv
  09_tables/T5_novelty_priority_targets.tsv
"""
import os
import numpy as np, pandas as pd

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality"); PROC = os.path.join(ROOT, "02_data_processed")
TAB  = os.path.join(ROOT, "09_tables"); os.makedirs(TAB, exist_ok=True)

hits = pd.read_csv(os.path.join(GEN, "phenome_hits.tsv"), sep="\t")
ann  = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")

# coloc (combine autoimmune + pan-phenome)
def load_coloc(fp):
    if not os.path.exists(fp): return pd.DataFrame(columns=["gene","disease","PP_H4"])
    return pd.read_csv(fp, sep="\t")[["gene","disease","PP_H4"]]
col = pd.concat([load_coloc(os.path.join(GEN,"coloc_results.tsv")),
                 load_coloc(os.path.join(GEN,"coloc_phenome_results.tsv"))],
                ignore_index=True).dropna(subset=["PP_H4"])
col = col.sort_values("PP_H4",ascending=False).drop_duplicates(["gene","disease"])
col = col.rename(columns={"gene":"gene_symbol"})

df = hits.merge(col, on=["gene_symbol","disease"], how="left")

# --- druggability from HPA protein class ---
acls = ann.set_index("gene_symbol").hpa_protein_class.fillna("").astype(str)
def drug_flag(g, key): return int(key.lower() in acls.get(g,"").lower())
df["fda_target"]   = df.gene_symbol.map(lambda g: drug_flag(g,"FDA approved drug target"))
df["is_secreted"]  = df.gene_symbol.map(lambda g: drug_flag(g,"secreted"))
df["is_membrane"]  = df.gene_symbol.map(lambda g: drug_flag(g,"membrane"))
df["druggability"] = (2*df.fda_target + df.is_secreted + df.is_membrane).clip(upper=3)

# --- cross-category pleiotropy bonus ---
gcat = hits.groupby("gene_symbol").disease_category.nunique()
df["n_categories"] = df.gene_symbol.map(gcat).fillna(1).astype(int)

# --- cell-source specificity (has a defined blood-cell source) ---
src = ann.set_index("gene_symbol").immune_source_cells.fillna("")
df["has_cellsource"] = df.gene_symbol.map(lambda g: int(str(src.get(g,"")).strip()!=""))

# --- known-drug / MHC penalties ---
KNOWN = {"IL6ST","CTLA4","TNFRSF1A","IL4","TNFSF14","TNFRSF14","TNFRSF1B","IL6R","TYK2","JAK2"}
df["known_axis"]  = df.gene_symbol.isin(KNOWN).astype(int)
df["mhc_penalty"] = df.get("mhc", pd.Series(False,index=df.index)).astype(int)

# ============================================================
# scoring (each term bounded, transparent)
# ============================================================
df["s_causal"]  = (-np.log10(df.FDR.clip(lower=1e-300))).clip(upper=10) / 5.0   # 0..2
df["s_coloc"]   = df.PP_H4.fillna(0)                                            # 0..1
df["s_pleio"]   = (df.n_categories - 1).clip(upper=3) * 0.5                     # 0..1.5
df["s_drug"]    = df.druggability * 0.5                                         # 0..1.5
df["s_cell"]    = df.has_cellsource * 0.5                                       # 0..0.5
df["p_known"]   = df.known_axis * 1.5
df["p_mhc"]     = df.mhc_penalty * 2.0
df["novelty_priority"] = (df.s_causal + df.s_coloc + df.s_pleio + df.s_drug
                          + df.s_cell - df.p_known - df.p_mhc)

df["category_label"] = np.where(df.known_axis==1, "recovered known axis",
                        np.where(df.mhc_penalty==1, "MHC-caution",
                        np.where(df.PP_H4.fillna(0)>=0.8, "NOVEL colocalized", "novel nomination")))

rank = df.sort_values("novelty_priority", ascending=False)
keep = ["gene_symbol","disease","disease_category","immune_class","OR","FDR","PP_H4",
        "n_categories","druggability","fda_target","category_label","novelty_priority",
        "s_causal","s_coloc","s_pleio","s_drug","s_cell","p_known","p_mhc"]
rank[keep].to_csv(os.path.join(GEN,"novelty_engine_ranked.tsv"), sep="\t", index=False)

# publication table: top novel targets
pub = rank[rank.category_label.str.startswith("NOVEL") | (rank.category_label=="novel nomination")]
pubcols = ["gene_symbol","disease","disease_category","immune_class","OR","FDR","PP_H4",
           "n_categories","druggability","category_label","novelty_priority"]
pub[pubcols].head(40).to_csv(os.path.join(TAB,"T5_novelty_priority_targets.tsv"), sep="\t", index=False)

print("=== NOVELTY PRIORITY ENGINE ===")
print(f"scored pairs: {len(df)} | novel colocalized: {(df.category_label=='NOVEL colocalized').sum()} "
      f"| known-axis: {df.known_axis.sum()} | MHC-caution: {df.mhc_penalty.sum()}")
print("\nTop 25 by novelty priority:")
show = ["gene_symbol","disease","disease_category","OR","FDR","PP_H4","n_categories",
        "druggability","category_label","novelty_priority"]
with pd.option_context("display.width",220,"display.max_columns",30):
    print(rank[show].head(25).to_string(index=False))
print("\nwrote novelty_engine_ranked.tsv + T5_novelty_priority_targets.tsv")
