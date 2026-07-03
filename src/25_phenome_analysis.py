#!/usr/bin/env python
"""
HDDM Layer 54 - Step 25
Pan-phenome analysis on the expanded immune cis-MR (src/24 output). Produces the
real-data novelty tables that turn the atlas into a multi-disease causal engine:

  1. cross-category PLEIOTROPY: immune genes causal (FDR<0.05) in >=2 disease
     CATEGORIES (e.g. autoimmune + cardiovascular) -> shared immune axes.
  2. CELL-SOURCE causal map: enrichment of causal targets by blood-cell source.
  3. DIRECTION-AWARE therapeutic map across the whole phenome.
  4. per-(gene,disease) table with FDR, OR, direction, category.

Outputs (06_genetic_causality/):
  phenome_pleiotropy_axes.tsv
  phenome_cellsource_map.tsv
  phenome_direction_map.tsv
  phenome_hits.tsv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality")
PROC = os.path.join(ROOT, "02_data_processed")

mr  = pd.read_csv(os.path.join(GEN, "cis_MR_phenome_results.tsv"), sep="\t")
ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
imm = ann[ann.is_plasma_immune == 1].copy()

sig = mr[mr.FDR < 0.05].copy()
sig["direction"] = np.where(sig.OR > 1, "risk (block)", "protective (agonise)")
sig["mhc"] = sig.chr.astype(str).eq("6") & sig.gene_symbol.isin(  # crude MHC flag
    imm[imm.HLA_antigen == 1].gene_symbol) | sig.gene_symbol.isin(["C2","C4A","C4B"])
sig.to_csv(os.path.join(GEN, "phenome_hits.tsv"), sep="\t", index=False)
print(f"pan-phenome FDR<0.05 hits: {len(sig)} across {sig.disease.nunique()} diseases, "
      f"{sig.gene_symbol.nunique()} genes")
print(sig.disease_category.value_counts().to_string())

# ============================================================
# 1. cross-category pleiotropy
# ============================================================
g_cat = sig.groupby("gene_symbol").disease_category.nunique()
pleio_genes = g_cat[g_cat >= 2].index.tolist()
rows = []
for g in pleio_genes:
    s = sig[sig.gene_symbol == g]
    cats = sorted(s.disease_category.unique())
    for _, r in s.iterrows():
        rows.append(dict(gene_symbol=g, n_categories=len(cats), categories=";".join(cats),
                         disease=r.disease, disease_category=r.disease_category,
                         OR=round(r.OR, 3), direction=r.direction, FDR=r.FDR))
pleio = pd.DataFrame(rows).sort_values(["n_categories", "gene_symbol"], ascending=[False, True])
pleio.to_csv(os.path.join(GEN, "phenome_pleiotropy_axes.tsv"), sep="\t", index=False)
print(f"\ncross-category pleiotropic genes (>=2 categories): {len(pleio_genes)}")
if len(pleio): print(pleio.to_string(index=False))

# ============================================================
# 2. cell-source causal enrichment (Fisher)
# ============================================================
def explode_src(df):
    s = df[["gene_symbol", "immune_source_cells"]].dropna().drop_duplicates("gene_symbol")
    s = s.assign(src=s.immune_source_cells.str.split(r"[;,]")).explode("src")
    s["src"] = s.src.str.strip()
    return s[s.src != ""]

bg = explode_src(imm)                       # background: all immune proteins w/ a source
causal_genes = set(sig.gene_symbol)
Ntot = bg.gene_symbol.nunique(); Ncaus = len(causal_genes & set(bg.gene_symbol))
rows = []
for src, sub in bg.groupby("src"):
    cls = set(sub.gene_symbol)
    a = len(cls & causal_genes); b = len(cls) - a; c = Ncaus - a; d = Ntot - len(cls) - c
    if a == 0: continue
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    rows.append(dict(source_cell=src, n_source=len(cls), n_causal=a,
                     hit_rate=a/len(cls), odds_ratio=orr, p=p))
cell = pd.DataFrame(rows).sort_values("p")
cell.to_csv(os.path.join(GEN, "phenome_cellsource_map.tsv"), sep="\t", index=False)
print("\n=== cell-source causal enrichment ===")
print(cell.to_string(index=False))

# ============================================================
# 3. direction-aware therapeutic map
# ============================================================
KNOWN = {  # literature-curated approved/known immune drug direction
 "IL6ST":"block","CTLA4":"agonise","TNFRSF1A":"block","IL4":"block",
 "TNFSF14":"block","TNFRSF14":"block",
}
dm = sig[["gene_symbol","disease","disease_category","OR","direction","FDR","mhc"]].copy()
dm["genetic_strategy"] = np.where(dm.OR > 1, "block/neutralise", "agonise/replace")
dm["known_drug_dir"] = dm.gene_symbol.map(KNOWN).fillna("-")
dm["status"] = np.where(dm.mhc, "MHC-caution (LD)",
                np.where(dm.known_drug_dir != "-", "recovers known axis", "novel nomination"))
dm = dm.sort_values(["disease_category","FDR"])
dm.to_csv(os.path.join(GEN, "phenome_direction_map.tsv"), sep="\t", index=False)
print(f"\ndirection map rows: {len(dm)} | novel: {(dm.status=='novel nomination').sum()} "
      f"| known-axis: {(dm.status=='recovers known axis').sum()} | MHC: {(dm.status=='MHC-caution (LD)').sum()}")
print("\nwrote phenome_pleiotropy_axes / cellsource_map / direction_map / hits")
