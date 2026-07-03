#!/usr/bin/env python
"""
HDDM Layer 54 - Step 2
Generate the descriptive figure set from downloaded PUBLIC data:
  plasma immune universe + HPA + MSigDB + FinnGen manifest + GWAS Catalog index.
All figures -> 08_figures/main_figures/*.png (300 dpi) and a combined panel.
"""
import os, json, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import seaborn as sns

sns.set_theme(style="whitegrid", context="talk")
PAL = sns.color_palette("Set2")
ROOT = r"I:\Plasma immune atalas"
RAW  = os.path.join(ROOT, "01_data_raw")
PROC = os.path.join(ROOT, "02_data_processed")
FIG  = os.path.join(ROOT, "08_figures", "main_figures")
os.makedirs(FIG, exist_ok=True)

ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
imm = ann[ann.is_plasma_immune == 1].copy()
summary = json.load(open(os.path.join(PROC, "immune_universe_summary.json")))

def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  wrote", name)

# ---------------------------------------------------------------- Fig 1 funnel
def fig1_funnel():
    stages = ["Olink Explore\nuniverse", "HPA-annotated", "MSigDB C7\nimmune genes",
              "Plasma immune\nproteins (kept)"]
    vals = [summary["olink_universe"], summary["hpa_matched"],
            summary["msigdb_c7_immune_genes_in_universe"], summary["plasma_immune_proteins"]]
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(range(len(stages))[::-1], vals, color=sns.color_palette("crest", len(stages)))
    for i, (s, v) in enumerate(zip(stages, vals)):
        y = len(stages) - 1 - i
        ax.text(v + 30, y, f"{v:,}", va="center", fontsize=14, fontweight="bold")
    ax.set_yticks(range(len(stages))[::-1]); ax.set_yticklabels(stages)
    ax.set_xlabel("Number of proteins")
    ax.set_title("Fig.1  Definition of the plasma immune protein universe", fontweight="bold")
    ax.set_xlim(0, max(vals) * 1.15)
    save(fig, "Fig01_immune_universe_funnel.png")

# ---------------------------------------------------------------- Fig 2 classes
def fig2_classes():
    cc = pd.Series(summary["immune_class_counts"]).sort_values()
    fig, ax = plt.subplots(figsize=(10, 7))
    cc.plot.barh(ax=ax, color=sns.color_palette("flare", len(cc)))
    for i, v in enumerate(cc.values):
        ax.text(v + 2, i, str(v), va="center", fontsize=12, fontweight="bold")
    ax.set_xlabel("Number of plasma immune proteins")
    ax.set_title("Fig.2  Plasma immune protein classification", fontweight="bold")
    save(fig, "Fig02_immune_class_composition.png")

# ---------------------------------------------------------------- Fig 3 sources
def fig3_sources():
    sc = pd.Series(summary["source_cell_counts"]).sort_values()
    fig, ax = plt.subplots(figsize=(9, 6))
    sc.plot.barh(ax=ax, color=sns.color_palette("mako", len(sc)))
    for i, v in enumerate(sc.values):
        ax.text(v + 3, i, str(v), va="center", fontsize=12, fontweight="bold")
    ax.set_xlabel("Number of immune proteins enriched in cell type")
    ax.set_title("Fig.3  Immune-cell source map (HPA blood-cell enrichment)", fontweight="bold")
    save(fig, "Fig03_immune_cell_source_map.png")

# ---------------------------------------------------------------- Fig 4 heatmap class x source
def fig4_class_source_heatmap():
    rows = []
    for _, r in imm.iterrows():
        if isinstance(r.immune_source_cells, str) and r.immune_source_cells:
            for sccell in r.immune_source_cells.split(";"):
                rows.append((r.immune_class, sccell))
    d = pd.DataFrame(rows, columns=["immune_class", "source"])
    ct = pd.crosstab(d.immune_class, d.source)
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(ct, annot=True, fmt="d", cmap="rocket_r", ax=ax, cbar_kws={"label": "proteins"})
    ax.set_title("Fig.4  Immune class \u00d7 source-cell matrix", fontweight="bold")
    ax.set_xlabel("Source cell lineage"); ax.set_ylabel("Immune class")
    save(fig, "Fig04_class_by_source_heatmap.png")

# ---------------------------------------------------------------- Fig 5 flag co-occurrence
def fig5_flag_corr():
    flags = ["cytokine","chemokine","interferon_axis","TNF_axis","complement","coagulation",
             "checkpoint","CD_marker","acute_phase","HLA_antigen","immunoglobulin",
             "soluble_receptor","immune_cell_enriched","secreted_to_blood"]
    sub = ann[flags]
    corr = sub.corr()
    fig, ax = plt.subplots(figsize=(11, 9))
    sns.heatmap(corr, cmap="vlag", center=0, square=True, ax=ax,
                cbar_kws={"label": "Pearson r"}, linewidths=.5)
    ax.set_title("Fig.5  Co-occurrence structure of immune annotation flags", fontweight="bold")
    save(fig, "Fig05_flag_cooccurrence.png")

# ---------------------------------------------------------------- Fig 6 FinnGen landscape
def fig6_finngen():
    man = pd.read_csv(os.path.join(RAW, "FinnGen_GWAS", "finngen_R12_manifest.tsv"), sep="\t")
    man["chapter"] = man["category"].astype(str).str.replace(r"\(.*\)", "", regex=True).str.strip()
    top = man["chapter"].value_counts().head(15).sort_values()
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    top.plot.barh(ax=axes[0], color=sns.color_palette("viridis", len(top)))
    axes[0].set_title("FinnGen R12: endpoints per disease chapter (top 15)", fontsize=14, fontweight="bold")
    axes[0].set_xlabel("number of disease endpoints")
    cases = man["num_cases"].clip(lower=1)
    axes[1].hist(np.log10(cases), bins=40, color=PAL[2], edgecolor="white")
    axes[1].set_title("FinnGen R12: case-count distribution", fontsize=14, fontweight="bold")
    axes[1].set_xlabel("log10(num_cases)"); axes[1].set_ylabel("number of endpoints")
    fig.suptitle(f"Fig.6  FinnGen R12 disease-endpoint landscape (N={len(man):,} endpoints) for MR/coloc",
                 fontweight="bold", y=1.02)
    save(fig, "Fig06_finngen_landscape.png")

# ---------------------------------------------------------------- Fig 7 GWAS immune traits
def fig7_gwas_immune():
    g = pd.read_csv(os.path.join(RAW, "GWAS_Catalog", "gwas_catalog_studies.tsv"),
                    sep="\t", low_memory=False)
    col = "DISEASE/TRAIT" if "DISEASE/TRAIT" in g.columns else g.columns[7]
    imm_kw = re.compile(r"immun|inflamm|cytokine|interleukin|lymph|leukocyte|neutrophil|"
                        r"monocyte|eosinophil|C-reactive|autoimmun|rheumat|lupus|"
                        r"psoriasis|colitis|crohn|asthma|allerg|sepsis|infection", re.I)
    gi = g[g[col].astype(str).str.contains(imm_kw)]
    top = gi[col].value_counts().head(18).sort_values()
    fig, ax = plt.subplots(figsize=(11, 8))
    top.plot.barh(ax=ax, color=sns.color_palette("rocket", len(top)))
    ax.set_xlabel("number of GWAS studies")
    ax.set_title(f"Fig.7  Immune-related GWAS in GWAS Catalog\n"
                 f"({len(gi):,} immune studies of {len(g):,} total)", fontweight="bold")
    save(fig, "Fig07_gwas_immune_traits.png")

# ---------------------------------------------------------------- Fig 8 drug targets
def fig8_drug_targets():
    imm["is_drug_target"] = imm["hpa_protein_class"].fillna("").str.contains("drug target", case=False).astype(int)
    by = imm.groupby("immune_class")["is_drug_target"].agg(["sum", "count"])
    by["pct"] = 100 * by["sum"] / by["count"]
    by = by.sort_values("pct")
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(by.index, by["pct"], color=sns.color_palette("crest", len(by)))
    for i, (s, c) in enumerate(zip(by["sum"], by["count"])):
        ax.text(by["pct"].iloc[i] + 0.5, i, f"{int(s)}/{int(c)}", va="center", fontsize=11)
    ax.set_xlabel("% proteins flagged as potential drug target (HPA)")
    ax.set_title("Fig.8  Druggability of plasma immune classes", fontweight="bold")
    save(fig, "Fig08_druggability_by_class.png")

# ---------------------------------------------------------------- Fig 9 disease involvement
def fig9_disease():
    terms = []
    for d in imm["hpa_disease"].dropna():
        for t in str(d).split(","):
            t = t.strip()
            if t and t.lower() not in ("disease related genes",):
                terms.append(t)
    top = pd.Series(terms).value_counts().head(18).sort_values()
    fig, ax = plt.subplots(figsize=(11, 8))
    top.plot.barh(ax=ax, color=sns.color_palette("flare", len(top)))
    ax.set_xlabel("number of plasma immune proteins")
    ax.set_title("Fig.9  Disease involvement of plasma immune proteins (HPA)", fontweight="bold")
    save(fig, "Fig09_disease_involvement.png")

# ---------------------------------------------------------------- Fig 10 MSigDB coverage
def fig10_msigdb():
    by = imm.groupby("immune_class")["msigdb_c7_immune_gene"].mean().mul(100).sort_values()
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(by.index, by.values, color=sns.color_palette("mako", len(by)))
    ax.set_xlabel("% in MSigDB C7 ImmuneSigDB")
    ax.set_title("Fig.10  Immune-signature (MSigDB C7) coverage by class", fontweight="bold")
    ax.set_xlim(0, 105)
    save(fig, "Fig10_msigdb_coverage.png")

# ---------------------------------------------------------------- Fig 11 ligand-receptor families
def fig11_axes():
    axes_def = {
        "IL-6 / gp130 axis": r"^(IL6|IL6R|IL6ST|OSM|LIF|CNTF|IL11|IL27)",
        "TNF superfamily": r"^(TNF|TNFSF|TNFRSF|LTA|LTB)",
        "CXC chemokines": r"^CXCL",
        "CC chemokines": r"^CCL",
        "Interleukin-1 family": r"^(IL1|IL18|IL33|IL36|IL37|IL38)",
        "Interferon axis": r"^(IFN|IFNG|IFNA|IFNB)",
        "Complement": r"^(C[1-9]|CFB|CFH|CFD|CFI|MASP)",
        "Checkpoints": r"(PDCD1|CD274|CTLA4|LAG3|HAVCR2|TIGIT|ICOS|CD40)",
    }
    counts = {}
    for lab, pat in axes_def.items():
        counts[lab] = imm["gene_symbol"].str.contains(pat, case=False, regex=True, na=False).sum()
    s = pd.Series(counts).sort_values()
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(s.index, s.values, color=sns.color_palette("Set2", len(s)))
    for i, v in enumerate(s.values):
        ax.text(v + 0.3, i, str(int(v)), va="center", fontweight="bold")
    ax.set_xlabel("number of measured plasma immune proteins")
    ax.set_title("Fig.11  Key immune communication axes covered by Olink", fontweight="bold")
    save(fig, "Fig11_immune_axes.png")

# ---------------------------------------------------------------- Fig 12 schematic
def fig12_schematic():
    fig, ax = plt.subplots(figsize=(14, 8)); ax.axis("off")
    steps = ["Plasma immune\nprotein", "Prevalent &\nincident disease", "Pre-disease\ntrajectory (15y)",
             "Immune ageing\nwaves", "PIRS\nrisk score", "pQTL / MR /\ncoloc",
             "Immune-cell\nsource", "Receptor /\nreceiver cell", "scATAC disease\nregulation",
             "Therapeutic\ntarget"]
    cols = sns.color_palette("Spectral", len(steps))
    x = 0.04
    for i, (st, c) in enumerate(zip(steps, cols)):
        box = FancyBboxPatch((x, 0.42), 0.082, 0.16, boxstyle="round,pad=0.01",
                             fc=c, ec="black", lw=1.2, alpha=.9)
        ax.add_patch(box)
        ax.text(x + 0.041, 0.50, st, ha="center", va="center", fontsize=9, fontweight="bold")
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + 0.082, 0.50), (x + 0.096, 0.50),
                         arrowstyle="-|>", mutation_scale=14, color="grey"))
        x += 0.096
    ax.text(0.5, 0.85, "Fig.12  HDDM Layer 54 \u2014 Plasma Immunome\u2013Phenome Atlas model",
            ha="center", fontsize=16, fontweight="bold")
    ax.text(0.5, 0.20, "Public-data prototype built here: protein universe \u2192 class \u2192 source cell "
            "\u2192 communication axis \u2192 disease/GWAS \u2192 druggability\n"
            "Greyed steps (trajectory, ageing, PIRS, MR) require controlled UK Biobank Olink access.",
            ha="center", fontsize=10, style="italic", color="#444")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    save(fig, "Fig12_atlas_schematic.png")

if __name__ == "__main__":
    print("Generating figures ->", FIG)
    fig1_funnel(); fig2_classes(); fig3_sources(); fig4_class_source_heatmap()
    fig5_flag_corr(); fig6_finngen(); fig7_gwas_immune(); fig8_drug_targets()
    fig9_disease(); fig10_msigdb(); fig11_axes(); fig12_schematic()
    print("DONE:", len(os.listdir(FIG)), "files")
