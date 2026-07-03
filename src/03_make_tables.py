#!/usr/bin/env python
"""
HDDM Layer 54 - Step 3
Build supplementary tables for the report:
  T1 plasma immune universe (full annotation, immune only)
  T2 immune class summary
  T3 immune communication-axis shortlist (cytokines/chemokines/checkpoints/complement)
  T4 candidate druggable immune targets (drug-target flagged, with source cell + disease)
"""
import os, json
import pandas as pd

ROOT = r"I:\Plasma immune atalas"
PROC = os.path.join(ROOT, "02_data_processed")
TAB  = os.path.join(ROOT, "09_tables")
os.makedirs(TAB, exist_ok=True)

ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
imm = ann[ann.is_plasma_immune == 1].copy()

# T1 full immune universe
keep_cols = ["olink_id","gene_symbol","protein_name","immune_class","immune_score",
             "cytokine","chemokine","interferon_axis","TNF_axis","complement","coagulation",
             "checkpoint","CD_marker","acute_phase","HLA_antigen","immunoglobulin",
             "soluble_receptor","immune_cell_enriched","immune_source_cells",
             "secreted_to_blood","msigdb_c7_immune_gene","rna_blood_cell_specificity",
             "hpa_disease","hpa_protein_class"]
keep_cols = [c for c in keep_cols if c in imm.columns]
imm[keep_cols].sort_values(["immune_class","gene_symbol"]).to_csv(
    os.path.join(TAB, "T1_plasma_immune_universe.tsv"), sep="\t", index=False)

# T2 class summary
t2 = imm.groupby("immune_class").agg(
    n_proteins=("gene_symbol","size"),
    pct_drug_target=("hpa_protein_class",
        lambda s: round(100*s.fillna("").str.contains("drug target",case=False).mean(),1)),
    pct_secreted=("secreted_to_blood", lambda s: round(100*s.mean(),1)),
    pct_immune_cell_enriched=("immune_cell_enriched", lambda s: round(100*s.mean(),1)),
).sort_values("n_proteins", ascending=False)
t2.to_csv(os.path.join(TAB, "T2_immune_class_summary.tsv"), sep="\t")

# T3 communication-axis shortlist
axis = imm[(imm.cytokine|imm.chemokine|imm.checkpoint|imm.complement|imm.TNF_axis|imm.interferon_axis)>0]
axis[["gene_symbol","protein_name","immune_class","immune_source_cells","hpa_disease"]]\
    .sort_values("immune_class").to_csv(
    os.path.join(TAB, "T3_communication_axis_shortlist.tsv"), sep="\t", index=False)

# T4 druggable target shortlist
imm["is_drug_target"] = imm["hpa_protein_class"].fillna("").str.contains("drug target", case=False)
t4 = imm[imm.is_drug_target & (imm.immune_class!="Immune-cell-enriched")].copy()
t4 = t4[["gene_symbol","protein_name","immune_class","immune_source_cells",
         "secreted_to_blood","hpa_disease"]].sort_values(["immune_class","gene_symbol"])
t4.to_csv(os.path.join(TAB, "T4_druggable_immune_targets.tsv"), sep="\t", index=False)

print("T1 universe rows:", len(imm))
print("T2 classes:", len(t2))
print("T3 axis shortlist:", len(axis))
print("T4 druggable targets:", len(t4))
print(t2.to_string())
