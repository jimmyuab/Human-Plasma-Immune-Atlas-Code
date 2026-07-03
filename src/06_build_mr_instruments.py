#!/usr/bin/env python
"""
HDDM Layer 54 - Step 6
Build cis-eQTL instruments for the plasma immune genes from eQTLGen.
For each immune gene keep the single strongest cis-eQTL SNP (max |Z|) -> Wald-ratio cis-MR.
Output: 06_genetic_causality/immune_cis_eqtl_instruments.tsv
"""
import os, gzip
import pandas as pd

ROOT = r"I:\Plasma immune atalas"
RAW  = os.path.join(ROOT, "01_data_raw")
PROC = os.path.join(ROOT, "02_data_processed")
OUT  = os.path.join(ROOT, "06_genetic_causality")
os.makedirs(OUT, exist_ok=True)

ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
immune_genes = set(ann.loc[ann.is_plasma_immune == 1, "gene_symbol"].dropna())
print("immune genes:", len(immune_genes))

eqtl = os.path.join(RAW, "eQTLGen", "cis-eQTLs-FDR0.05.txt.gz")
best = {}   # gene -> best row dict (max |Z|)
n = 0
with gzip.open(eqtl, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(header)}
    for line in f:
        n += 1
        p = line.rstrip("\n").split("\t")
        g = p[idx["GeneSymbol"]]
        if g not in immune_genes:
            continue
        z = abs(float(p[idx["Zscore"]]))
        cur = best.get(g)
        if cur is None or z > cur["absZ"]:
            best[g] = {
                "gene_symbol": g,
                "SNP": p[idx["SNP"]],
                "chr": p[idx["SNPChr"]],
                "pos": p[idx["SNPPos"]],
                "effect_allele": p[idx["AssessedAllele"]],
                "other_allele": p[idx["OtherAllele"]],
                "Zscore": float(p[idx["Zscore"]]),
                "absZ": z,
                "NrSamples": int(p[idx["NrSamples"]]),
                "Pvalue": p[idx["Pvalue"]],
            }
print("rows scanned:", f"{n:,}", "| genes with cis instrument:", len(best))

inst = pd.DataFrame(best.values()).drop(columns=["absZ"])
inst.to_csv(os.path.join(OUT, "immune_cis_eqtl_instruments.tsv"), sep="\t", index=False)
print("wrote", len(inst), "instruments ->", os.path.join(OUT, "immune_cis_eqtl_instruments.tsv"))
print(inst.head())
