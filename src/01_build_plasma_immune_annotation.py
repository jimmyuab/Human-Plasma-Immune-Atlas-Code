#!/usr/bin/env python
"""
HDDM Layer 54 - Step 1
Build the plasma immune protein universe annotation table by integrating:
  - UKB coding 143  : full Olink Explore protein universe (gene;description)
  - HPA proteinatlas.tsv : immune-cell source, blood specificity, secretome, protein class, disease
  - MSigDB C7       : immune signature gene-set membership
Outputs:
  02_data_processed/plasma_immune_protein_annotation.tsv
  02_data_processed/immune_universe_summary.json
"""
import os, re, json, zipfile, io, csv
import pandas as pd
import numpy as np

ROOT = r"I:\Plasma immune atalas"
RAW  = os.path.join(ROOT, "01_data_raw")
OUT  = os.path.join(ROOT, "02_data_processed")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Olink universe (UKB coding 143)
# ---------------------------------------------------------------------------
olink = pd.read_csv(os.path.join(RAW, "UKB_Olink_NPX", "ukb_coding143_olink.tsv"),
                    sep="\t")
olink[["gene_symbol", "protein_name"]] = olink["meaning"].str.split(";", n=1, expand=True)
olink = olink.rename(columns={"coding": "olink_id"})[["olink_id", "gene_symbol", "protein_name"]]
olink["gene_symbol"] = olink["gene_symbol"].str.strip()
print(f"[Olink] universe proteins: {len(olink)}")

# ---------------------------------------------------------------------------
# 2. HPA master annotation
# ---------------------------------------------------------------------------
zf = zipfile.ZipFile(os.path.join(RAW, "HPA_blood_protein", "proteinatlas.tsv.zip"))
with zf.open("proteinatlas.tsv") as fh:
    hpa = pd.read_csv(io.TextIOWrapper(fh, "utf-8"), sep="\t", low_memory=False)
print(f"[HPA] genes: {len(hpa)}  cols: {hpa.shape[1]}")

hpa_cols = {
    "Gene": "gene_symbol",
    "Protein class": "hpa_protein_class",
    "Biological process": "hpa_bio_process",
    "Molecular function": "hpa_mol_function",
    "Disease involvement": "hpa_disease",
    "RNA blood cell specificity": "rna_blood_cell_specificity",
    "RNA blood cell specific nTPM": "rna_blood_cell_specific_ntpm",
    "RNA blood lineage specificity": "rna_blood_lineage_specificity",
    "Secretome location": "secretome_location",
    "RNA tissue specificity": "rna_tissue_specificity",
}
hpa_cols = {k: v for k, v in hpa_cols.items() if k in hpa.columns}
hpa_small = hpa[list(hpa_cols)].rename(columns=hpa_cols)

ann = olink.merge(hpa_small, on="gene_symbol", how="left")
print(f"[merge] HPA-matched: {ann['hpa_protein_class'].notna().sum()} / {len(ann)}")

# ---------------------------------------------------------------------------
# 3. MSigDB C7 immune-gene membership
# ---------------------------------------------------------------------------
c7_genes = set()
c7_path = os.path.join(RAW, "MSigDB_C7_C8", "c7.all.v2024.1.Hs.symbols.gmt")
n_sets = 0
with open(c7_path) as fh:
    for line in fh:
        n_sets += 1
        parts = line.rstrip("\n").split("\t")
        for g in parts[2:]:
            c7_genes.add(g.strip())
print(f"[MSigDB C7] sets={n_sets} unique_genes={len(c7_genes)}")
ann["msigdb_c7_immune_gene"] = ann["gene_symbol"].isin(c7_genes).astype(int)

# ---------------------------------------------------------------------------
# 4. Rule-based immune classification
# ---------------------------------------------------------------------------
def s(x):
    return "" if pd.isna(x) else str(x)

# Cytokine / chemokine / interleukin / interferon / TNF families by symbol
RE_CYTOKINE   = re.compile(r"^(IL\d|IL\dR|CSF|TGFB|TNF|TNFSF|TNFRSF|IFN|LIF|OSM|EBI3|TSLP|LTA|LTB)", re.I)
RE_INTERLEUKIN= re.compile(r"^IL\d", re.I)
RE_INTERFERON = re.compile(r"^(IFN|IRF|ISG|MX\d|OAS\d|STAT1)", re.I)
RE_TNF        = re.compile(r"(TNF|TNFSF|TNFRSF)", re.I)
RE_CHEMOKINE  = re.compile(r"^(CXCL|CCL|CX3CL|XCL|CXCR|CCR|CX3CR|XCR|ACKR|PF4|PPBP)", re.I)
RE_COMPLEMENT = re.compile(r"^(C[1-9]$|C[1-9][A-Z]|CFB|CFH|CFD|CFI|CFP|CR\d|CD55|CD46|CD59|MASP|MBL|SERPING1|C1Q|C4BP)", re.I)
RE_COAG       = re.compile(r"^(F\d{1,2}$|F\d{1,2}[A-Z]|PLG|PLAU|PLAT|SERPINC1|SERPINE1|VWF|PROC|PROS1|THBD|FGA|FGB|FGG|TFPI)", re.I)
RE_CHECKPOINT = re.compile(r"(PDCD1|CD274|CTLA4|LAG3|HAVCR2|TIGIT|ICOS|CD80|CD86|CD28|BTLA|VSIR|CD276|TNFRSF9|TNFRSF4|CD40|CD40LG)", re.I)
RE_CD         = re.compile(r"^CD\d", re.I)
RE_ACUTEPHASE = re.compile(r"^(CRP|SAA\d|SAA|LBP|HP|HPX|SERPINA3|ORM\d|FGA|FGB|FGG|FTL|FTH1|PCT|CALCA)", re.I)
RE_HLA        = re.compile(r"^(HLA|B2M|TAP\d|TAPBP|CD1)", re.I)
RE_IG         = re.compile(r"^(IGH|IGK|IGL|FCGR|FCER|FCRL|JCHAIN|POLR|AICDA|CD19|MS4A1|CD79)", re.I)

# Soluble-receptor heuristic: gene name contains 'receptor' AND on Olink (soluble form)
def is_soluble_receptor(name):
    return int(bool(re.search(r"receptor", s(name), re.I)))

infl_terms = re.compile(r"inflamm|cytokine|chemokine|interleukin|interferon|immun|complement|leukocyte|lymphocyte|T cell|B cell|innate|adaptive|antigen", re.I)

def immune_cell_enriched(spec):
    # HPA categorical field: 'Immune cell enriched' / 'enhanced' / 'Group enriched'
    sp = s(spec).lower()
    return int(("immune cell enriched" in sp) or ("immune cell enhanced" in sp)
               or ("group enriched" in sp))

# Map detailed HPA immune cell -> coarse source lineage
SOURCE_MAP = [
    (re.compile(r"monocyte|myeloid DC|plasmacytoid|macrophage|Kupffer", re.I), "Myeloid/Monocyte"),
    (re.compile(r"neutrophil|basophil|eosinophil|granulocyte", re.I),        "Granulocyte"),
    (re.compile(r"\bT-reg|T-cell|MAIT|gdT|Th1|Th2|Th17|CD4|CD8", re.I),       "T-cell"),
    (re.compile(r"B-cell|plasma cell|plasmablast", re.I),                     "B-cell"),
    (re.compile(r"NK-cell|natural killer", re.I),                            "NK-cell"),
    (re.compile(r"\bDC\b|dendritic", re.I),                                  "Dendritic"),
]
def source_cells(ntpm):
    txt = s(ntpm)
    found = []
    for rgx, lab in SOURCE_MAP:
        if rgx.search(txt) and lab not in found:
            found.append(lab)
    return ";".join(found)

ann["cytokine"]          = ann["gene_symbol"].apply(lambda g: int(bool(RE_CYTOKINE.search(s(g)))))
ann["interleukin"]       = ann["gene_symbol"].apply(lambda g: int(bool(RE_INTERLEUKIN.search(s(g)))))
ann["interferon_axis"]   = ann["gene_symbol"].apply(lambda g: int(bool(RE_INTERFERON.search(s(g)))))
ann["TNF_axis"]          = ann["gene_symbol"].apply(lambda g: int(bool(RE_TNF.search(s(g)))))
ann["chemokine"]         = ann["gene_symbol"].apply(lambda g: int(bool(RE_CHEMOKINE.search(s(g)))))
ann["complement"]        = ann["gene_symbol"].apply(lambda g: int(bool(RE_COMPLEMENT.search(s(g)))))
ann["coagulation"]       = ann["gene_symbol"].apply(lambda g: int(bool(RE_COAG.search(s(g)))))
ann["checkpoint"]        = ann["gene_symbol"].apply(lambda g: int(bool(RE_CHECKPOINT.search(s(g)))))
ann["CD_marker"]         = ann["gene_symbol"].apply(lambda g: int(bool(RE_CD.search(s(g)))))
ann["acute_phase"]       = ann["gene_symbol"].apply(lambda g: int(bool(RE_ACUTEPHASE.search(s(g)))))
ann["HLA_antigen"]       = ann["gene_symbol"].apply(lambda g: int(bool(RE_HLA.search(s(g)))))
ann["immunoglobulin"]    = ann["gene_symbol"].apply(lambda g: int(bool(RE_IG.search(s(g)))))
ann["soluble_receptor"]  = ann["protein_name"].apply(is_soluble_receptor)
ann["immune_cell_enriched"] = ann["rna_blood_cell_specificity"].apply(immune_cell_enriched)
ann["immune_source_cells"]  = ann["rna_blood_cell_specific_ntpm"].apply(source_cells)

# secreted / blood detectable
ann["secreted_to_blood"] = (
    ann["hpa_protein_class"].fillna("").str.contains("Secreted|Blood", case=False) |
    ann["secretome_location"].fillna("").str.contains("secreted", case=False)
).astype(int)

# Olink Inflammation panel proxy: HPA process / class flagged as inflammation/immune
ann["inflammation_immune_annotated"] = (
    ann["hpa_protein_class"].apply(lambda x: int(bool(infl_terms.search(s(x))))) |
    ann["hpa_bio_process"].apply(lambda x: int(bool(infl_terms.search(s(x)))))
)

# ---------------------------------------------------------------------------
# 5. Composite immune score & keep rule (per blueprint)
# ---------------------------------------------------------------------------
# NOTE: MSigDB C7 (ImmuneSigDB) covers 2,872/2,923 Olink genes -> too permissive to
# drive inclusion; kept as annotation only. Same for secreted_to_blood. Inclusion is
# driven by hard immune-biology flags + HPA categorical immune-cell enrichment.
ann["immune_score"] = (
      2 * ((ann["cytokine"] | ann["chemokine"] | ann["interferon_axis"] | ann["TNF_axis"]) > 0).astype(int)
    + 2 * ((ann["complement"] | ann["checkpoint"]) > 0).astype(int)
    + 2 * ann["immune_cell_enriched"]
    + 1 * ((ann["CD_marker"] | ann["HLA_antigen"] | ann["immunoglobulin"] |
            ann["acute_phase"] | ann["coagulation"]) > 0).astype(int)
    + 1 * ann["inflammation_immune_annotated"]
)
ann["is_plasma_immune"] = (ann["immune_score"] >= 2).astype(int)

# primary immune class label (priority order)
def immune_class(r):
    if r["complement"]:     return "Complement/Coagulation"
    if r["checkpoint"]:     return "Immune checkpoint"
    if r["chemokine"]:      return "Chemokine axis"
    if r["interferon_axis"]:return "Interferon axis"
    if r["TNF_axis"]:       return "TNF superfamily"
    if r["interleukin"] or r["cytokine"]: return "Cytokine/Interleukin"
    if r["acute_phase"]:    return "Acute-phase"
    if r["HLA_antigen"]:    return "HLA/Antigen presentation"
    if r["immunoglobulin"]: return "Ig/B-cell/Fc"
    if r["CD_marker"]:      return "CD/Leukocyte surface"
    if r["immune_cell_enriched"]: return "Immune-cell-enriched"
    if r["soluble_receptor"]:return "Soluble receptor"
    if r["inflammation_immune_annotated"]: return "Other immune-annotated"
    return "Non-immune"
ann["immune_class"] = ann.apply(immune_class, axis=1)

out_path = os.path.join(OUT, "plasma_immune_protein_annotation.tsv")
ann.to_csv(out_path, sep="\t", index=False)

immune = ann[ann["is_plasma_immune"] == 1]
summary = {
    "olink_universe": int(len(ann)),
    "hpa_matched": int(ann["hpa_protein_class"].notna().sum()),
    "msigdb_c7_immune_genes_in_universe": int(ann["msigdb_c7_immune_gene"].sum()),
    "plasma_immune_proteins": int(len(immune)),
    "immune_class_counts": immune["immune_class"].value_counts().to_dict(),
    "source_cell_counts": immune.loc[immune["immune_source_cells"]!="", "immune_source_cells"]
                                .str.split(";").explode().value_counts().to_dict(),
    "flag_counts": {c: int(ann[c].sum()) for c in
        ["cytokine","chemokine","interferon_axis","TNF_axis","complement","coagulation",
         "checkpoint","CD_marker","acute_phase","HLA_antigen","immunoglobulin",
         "soluble_receptor","immune_cell_enriched","secreted_to_blood","msigdb_c7_immune_gene"]},
}
with open(os.path.join(OUT, "immune_universe_summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)

print("\n=== SUMMARY ===")
print(json.dumps(summary, indent=2))
print(f"\nWrote {out_path}")
