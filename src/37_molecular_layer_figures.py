#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
37 - Immune molecular-LAYER figure gallery (antibody / RNA / secretome-EV)
==========================================================================
The plasma-immunome atlas is not only a "protein" resource: each curated
plasma-immune protein already carries, from Human Protein Atlas + our curation,
several orthogonal *molecular layers* that are themselves immune biology:

  * the ANTIBODY layer  -> immunoglobulin / B-cell / Fc genes (Ig axis)
  * functional immune AXES -> cytokine, interleukin, interferon, TNF, chemokine,
    complement, coagulation, checkpoint, CD/leukocyte-surface, acute-phase, HLA,
    soluble-receptor
  * the RNA layer -> HPA blood-cell & lineage RNA specificity of each gene
  * the SECRETOME / EV route -> how each protein reaches plasma
    (secreted-to-blood vs membrane / intracellular = vesicle / leakage route)

This module surfaces those layers as publication-level figures. It is purely
additive: it reads existing tables only, writes new PNGs into
08_figures/paper_style/, and modifies nothing. No new molecular data are
fabricated -- every value is already present in the curated annotation or the
causal-MR result table.

Run:
    python src/37_molecular_layer_figures.py
"""

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GC   = os.path.join(ROOT, "06_genetic_causality")
PROC = os.path.join(ROOT, "02_data_processed")
OUT  = os.path.join(ROOT, "08_figures", "paper_style")
os.makedirs(OUT, exist_ok=True)

FDR_SIG = 0.05

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333",
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "legend.frameon": False,
})

CAT_COLORS = {
    "Autoimmune":     "#D1495B",
    "Cardiovascular": "#2E6F95",
    "Metabolic":      "#E9963A",
    "Neuro/Aging":    "#6A4C93",
    "Renal":          "#2A9D8F",
}


def load():
    ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
    mr  = pd.read_csv(os.path.join(GC, "cis_MR_phenome_results.tsv"), sep="\t")
    return ann, mr


def _flag(series):
    """Coerce a mixed 0/1/'True' annotation column to a boolean Series."""
    s = series.copy()
    if s.dtype == bool:
        return s.fillna(False)
    return s.astype(str).str.lower().isin(["1", "true", "yes", "y"])


# --------------------------------------------------------------------------- #
def panel_ml1_axes(ann, mr):
    """ML1 - the atlas already spans every immune molecular axis (antibody axis
    included), and each axis is genetically instrumented + causally tested."""
    im = ann[ann["is_plasma_immune"] == True].copy()
    inst_genes = set(mr["gene_symbol"].unique())
    hit_genes  = set(mr.loc[mr["FDR"] < FDR_SIG, "gene_symbol"].unique())

    axes = [
        ("immunoglobulin",  "Antibody (Ig / B-cell / Fc)"),
        ("cytokine",        "Cytokine"),
        ("interleukin",     "Interleukin"),
        ("interferon_axis", "Interferon axis"),
        ("TNF_axis",        "TNF superfamily"),
        ("chemokine",       "Chemokine"),
        ("complement",      "Complement"),
        ("coagulation",     "Coagulation"),
        ("checkpoint",      "Immune checkpoint"),
        ("CD_marker",       "CD / leukocyte surface"),
        ("acute_phase",     "Acute-phase"),
        ("HLA_antigen",     "HLA / antigen presentation"),
        ("soluble_receptor","Soluble receptor"),
    ]
    rows = []
    for col, label in axes:
        if col not in im.columns:
            continue
        genes = set(im.loc[_flag(im[col]), "gene_symbol"])
        rows.append((label,
                     len(genes),
                     len(genes & inst_genes),
                     len(genes & hit_genes)))
    df = pd.DataFrame(rows, columns=["axis", "curated", "instrumented", "causal"])
    df = df.sort_values("curated", ascending=True)

    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    y = np.arange(len(df))
    ax.barh(y, df["curated"],      color="#DDE3EA", edgecolor="#B7C1CC",
            linewidth=0.5, label="curated in atlas")
    ax.barh(y, df["instrumented"], color="#7FA8C9", edgecolor="white",
            linewidth=0.5, label="cis-eQTL instrumented")
    ax.barh(y, df["causal"],       color="#C44E52", edgecolor="white",
            linewidth=0.5, label="causal hit (FDR<0.05)")
    for i, r in enumerate(df.itertuples()):
        ax.text(r.curated + 1.5, i, f"{r.curated}", va="center",
                fontsize=7.6, color="#555555")
    # highlight the antibody axis (the layer the user asked to add)
    for i, lab in enumerate(df["axis"]):
        if lab.startswith("Antibody"):
            ax.get_yticklabels()  # ensure ticks built
    ax.set_yticks(y)
    ax.set_yticklabels(df["axis"], fontsize=8.4)
    ax.set_xlabel("plasma-immune proteins")
    ax.set_title("Immune molecular axes in the atlas \u2014 antibody axis included,\n"
                 "each axis genetically instrumented and causally tested",
                 loc="left", fontsize=10.5)
    ax.legend(fontsize=8, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "ML1_immune_axes.png"))
    plt.close(fig)
    print("  ML1_immune_axes.png")


# --------------------------------------------------------------------------- #
def panel_ml2_rna(ann, mr):
    """ML2 - the RNA layer: HPA blood-cell + lineage RNA specificity of the
    genetically-instrumented immune proteins."""
    im = ann[ann["is_plasma_immune"] == True].copy()
    inst = im[im["gene_symbol"].isin(set(mr["gene_symbol"].unique()))].copy()

    cell_order = ["Immune cell enriched", "Immune cell enhanced", "Group enriched",
                  "Low immune cell specificity", "Not detected in immune cells"]
    lin_order = ["Lineage enriched", "Lineage enhanced", "Group enriched",
                 "Low lineage specificity", "Not detected"]
    cell_col = "#2E6F95"
    palette = ["#1B4965", "#2E6F95", "#62929E", "#BFD7EA", "#E6ECF2"]

    def counts(col, order):
        c = inst[col].fillna("Not detected").value_counts()
        return [int(c.get(k, 0)) for k in order]

    cell_c = counts("rna_blood_cell_specificity", cell_order)
    lin_c  = counts("rna_blood_lineage_specificity", lin_order)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    for ax, order, vals, ttl in [
        (axes[0], cell_order, cell_c, "Blood-cell RNA specificity"),
        (axes[1], lin_order,  lin_c,  "Blood-lineage RNA specificity"),
    ]:
        y = np.arange(len(order))
        ax.barh(y, vals, color=palette, edgecolor="white", linewidth=0.6)
        for i, v in enumerate(vals):
            if v:
                ax.text(v + max(vals) * 0.01, i, str(v), va="center", fontsize=8)
        ax.set_yticks(y); ax.set_yticklabels(order, fontsize=8.2)
        ax.invert_yaxis()
        ax.set_xlabel("instrumented immune proteins")
        ax.set_title(ttl, loc="left")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("RNA layer \u2014 blood-cell transcriptomic origin of the atlas immune proteins "
                 "(Human Protein Atlas)", fontsize=11, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(OUT, "ML2_rna_bloodcell_layer.png"))
    plt.close(fig)
    print("  ML2_rna_bloodcell_layer.png")


# --------------------------------------------------------------------------- #
def panel_ml3_secretome(ann, mr):
    """ML3 - the secretome / EV route: how each instrumented immune protein
    reaches the plasma (classical secretion vs membrane / intracellular =
    extracellular-vesicle / leakage route)."""
    im = ann[ann["is_plasma_immune"] == True].copy()
    inst = im[im["gene_symbol"].isin(set(mr["gene_symbol"].unique()))].copy()

    def bucket(v):
        s = str(v)
        if s == "nan" or s.strip() == "":
            return "Unannotated route"
        if s.startswith("Secreted to blood"):
            return "Secreted to blood (classical)"
        if s.startswith("Intracellular"):
            return "Membrane / intracellular\n(vesicle / EV / leakage route)"
        if s.startswith("Secreted to extracellular"):
            return "Secreted to extracellular matrix"
        if s.startswith("Immunoglobulin"):
            return "Immunoglobulin genes"
        if s.startswith("Secreted"):
            return "Secreted in other tissues"
        return "Unannotated route"

    b = inst["secretome_location"].map(bucket)
    order = ["Secreted to blood (classical)",
             "Membrane / intracellular\n(vesicle / EV / leakage route)",
             "Secreted to extracellular matrix",
             "Secreted in other tissues",
             "Immunoglobulin genes",
             "Unannotated route"]
    vals = [int((b == k).sum()) for k in order]
    order = [o for o, v in zip(order, vals) if v > 0]
    vals  = [v for v in vals if v > 0]
    cols = ["#2A9D8F", "#E9963A", "#8172B3", "#B7C1CC", "#8C8C8C", "#DDE3EA"][:len(order)]

    fig, ax = plt.subplots(figsize=(7.8, 6.4))
    wedges, _ = ax.pie(vals, colors=cols, startangle=90,
                       wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.4))
    ax.set(aspect="equal")
    tot = sum(vals)
    ax.text(0, 0, f"{tot}\nimmune\nproteins", ha="center", va="center",
            fontsize=12, fontweight="bold")
    leg = [f"{o.replace(chr(10),' ')}  \u2014 {v} ({v/tot*100:.0f}%)"
           for o, v in zip(order, vals)]
    ax.legend(wedges, leg, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              fontsize=8.2, ncol=1)
    ax.set_title("Secretome / EV route \u2014 how atlas immune proteins reach the plasma",
                 loc="center", pad=14)
    fig.savefig(os.path.join(OUT, "ML3_secretome_ev_route.png"))
    plt.close(fig)
    print("  ML3_secretome_ev_route.png")


# --------------------------------------------------------------------------- #
def panel_ml4_antibody(ann, mr):
    """ML4 - antibody-layer spotlight: the immunoglobulin / B-cell / Fc genes,
    whether each is genetically instrumented, and its strongest causal signal."""
    im = ann[ann["is_plasma_immune"] == True].copy()
    ig = im[_flag(im["immunoglobulin"]) | (im.get("immune_class") == "Ig/B-cell/Fc")].copy()
    ig_genes = sorted(set(ig["gene_symbol"]))
    inst_genes = set(mr["gene_symbol"].unique())

    rows = []
    for g in ig_genes:
        sub = mr[mr["gene_symbol"] == g]
        if len(sub):
            best = sub.loc[sub["MR_p"].idxmin()]
            rows.append((g, True, -np.log10(max(best["MR_p"], 1e-300)),
                         best["disease"], best["disease_category"],
                         bool(best["FDR"] < FDR_SIG)))
        else:
            rows.append((g, False, 0.0, "", "", False))
    df = pd.DataFrame(rows, columns=["gene", "inst", "logp", "disease", "cat", "hit"])
    df = df.sort_values(["inst", "logp"], ascending=[True, True])

    fig, ax = plt.subplots(figsize=(8.6, max(4.2, 0.34 * len(df) + 1.4)))
    y = np.arange(len(df))
    for i, r in enumerate(df.itertuples()):
        if not r.inst:
            ax.barh(i, 0.2, color="#E6ECF2", edgecolor="#C4CDD6", linewidth=0.5)
            ax.text(0.28, i, "no cis-eQTL instrument", va="center",
                    fontsize=7, color="#9AA5B1", style="italic")
        else:
            col = CAT_COLORS.get(r.cat, "#7FA8C9")
            ax.barh(i, r.logp, color=col, edgecolor="white", linewidth=0.5,
                    alpha=0.95 if r.hit else 0.5)
            tag = f"{r.disease}" + ("  \u2605" if r.hit else "")
            ax.text(r.logp + 0.08, i, tag, va="center", fontsize=7,
                    color="#222222" if r.hit else "#666666")
    ax.set_yticks(y); ax.set_yticklabels(df["gene"], fontsize=8)
    ax.set_xlabel("strongest causal signal  \u2212log$_{10}$($P_{MR}$)")
    ax.set_title("Antibody layer \u2014 immunoglobulin / B-cell / Fc genes in the causal atlas\n"
                 "(\u2605 = FDR<0.05; bar colour = disease category)",
                 loc="left", fontsize=10.3)
    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                      markerfacecolor=CAT_COLORS[c], markeredgecolor="white", label=c)
               for c in CAT_COLORS if c in set(df.loc[df["inst"], "cat"])]
    if handles:
        ax.legend(handles=handles, fontsize=7.4, loc="lower right",
                  title="disease category", title_fontsize=7.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(OUT, "ML4_antibody_layer.png"))
    plt.close(fig)
    print("  ML4_antibody_layer.png")


# --------------------------------------------------------------------------- #
def panel_ml5_class_heatmap(mr):
    """ML5 - immune molecular-class x disease-category causal heatmap: which
    molecular layer drives which disease chapter."""
    hits = mr[mr["FDR"] < FDR_SIG].copy()
    cats = ["Autoimmune", "Cardiovascular", "Metabolic", "Neuro/Aging", "Renal"]
    cats = [c for c in cats if c in set(hits["disease_category"])]
    classes = (hits.groupby("immune_class").size().sort_values(ascending=False).index.tolist())
    tab = (hits.groupby(["immune_class", "disease_category"]).size()
                .unstack(fill_value=0).reindex(index=classes, columns=cats, fill_value=0))

    cmap = LinearSegmentedColormap.from_list("imm", ["#F7FBFC", "#2E6F95", "#12303F"])
    fig, ax = plt.subplots(figsize=(7.6, 0.5 * len(classes) + 2.2))
    im = ax.imshow(tab.values, cmap=cmap, aspect="auto")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, rotation=30, ha="right", fontsize=8.6)
    for tick, c in zip(ax.get_xticklabels(), cats):
        tick.set_color(CAT_COLORS.get(c, "#333"))
        tick.set_fontweight("bold")
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes, fontsize=8.4)
    mx = tab.values.max()
    for i in range(len(classes)):
        for j in range(len(cats)):
            v = tab.values[i, j]
            if v:
                ax.text(j, i, int(v), ha="center", va="center", fontsize=8,
                        color="white" if v > mx * 0.5 else "#12303F")
    ax.set_title("Which immune molecular layer drives which disease category\n"
                 "(causal gene\u2013disease hits, FDR<0.05)", loc="left", fontsize=10.3)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("causal hits")
    fig.savefig(os.path.join(OUT, "ML5_class_by_category_heatmap.png"))
    plt.close(fig)
    print("  ML5_class_by_category_heatmap.png")


def main():
    print("[37] Immune molecular-layer figure gallery ...")
    ann, mr = load()
    print(f"    annotation: {len(ann)} rows, {int((ann['is_plasma_immune']==True).sum())} plasma-immune")
    print(f"    MR: {len(mr):,} tests, {mr['gene_symbol'].nunique()} instrumented genes")
    panel_ml1_axes(ann, mr)
    panel_ml2_rna(ann, mr)
    panel_ml3_secretome(ann, mr)
    panel_ml4_antibody(ann, mr)
    panel_ml5_class_heatmap(mr)
    print(f"[37] Done -> {OUT}")


if __name__ == "__main__":
    main()
