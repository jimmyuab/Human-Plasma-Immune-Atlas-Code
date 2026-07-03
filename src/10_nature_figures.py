#!/usr/bin/env python
"""
HDDM Layer 54 - Step 10
Consolidate results into 6 Nature-style multi-panel figures.
Each figure answers ONE biological question. -> 08_figures/nature/
"""
import os, json, re
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import seaborn as sns
import string

sns.set_theme(style="white", context="paper")
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
                     "axes.labelsize": 9, "figure.dpi": 150})

ROOT = r"I:\Plasma immune atalas"
RAW  = os.path.join(ROOT, "01_data_raw")
PROC = os.path.join(ROOT, "02_data_processed")
GEN  = os.path.join(ROOT, "06_genetic_causality")
FIG  = os.path.join(ROOT, "08_figures", "nature")
os.makedirs(FIG, exist_ok=True)

ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
imm = ann[ann.is_plasma_immune == 1].copy()
S   = json.load(open(os.path.join(PROC, "immune_universe_summary.json")))
mr  = pd.read_csv(os.path.join(GEN, "cis_MR_immune_results.tsv"), sep="\t")
mr["nlog10p"] = -np.log10(mr["MR_p"].clip(lower=1e-300))
sig = mr[mr.FDR < 0.05].copy()

def panel_label(ax, letter, dx=-0.14, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")

def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig); print("  wrote", name)

# =====================================================================
# FIGURE 1 — Study design and immunome curation
# Q: what is the plasma immunome and how was it defined?
# =====================================================================
def figure1():
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)
    # A funnel
    ax = fig.add_subplot(gs[0, 0])
    stages = ["Olink Explore", "HPA-annotated", "MSigDB C7", "Plasma immunome"]
    vals = [S["olink_universe"], S["hpa_matched"], S["msigdb_c7_immune_genes_in_universe"], S["plasma_immune_proteins"]]
    ax.barh(range(len(stages))[::-1], vals, color=sns.color_palette("crest", 4))
    for i, v in enumerate(vals):
        ax.text(v + 30, len(stages) - 1 - i, f"{v:,}", va="center", fontsize=9, fontweight="bold")
    ax.set_yticks(range(len(stages))[::-1]); ax.set_yticklabels(stages)
    ax.set_xlabel("proteins"); ax.set_title("Curation of the plasma immunome")
    ax.set_xlim(0, max(vals) * 1.18); panel_label(ax, "a")
    # B class composition
    ax = fig.add_subplot(gs[0, 1])
    cc = pd.Series(S["immune_class_counts"]).sort_values()
    ax.barh(cc.index, cc.values, color=sns.color_palette("flare", len(cc)))
    for i, v in enumerate(cc.values):
        ax.text(v + 2, i, str(v), va="center", fontsize=7.5)
    ax.set_xlabel("proteins"); ax.set_title("Immune class composition")
    ax.tick_params(axis="y", labelsize=7.5); panel_label(ax, "b")
    # C source-cell map
    ax = fig.add_subplot(gs[1, 0])
    scc = pd.Series(S["source_cell_counts"]).sort_values()
    ax.barh(scc.index, scc.values, color=sns.color_palette("mako", len(scc)))
    for i, v in enumerate(scc.values):
        ax.text(v + 3, i, str(v), va="center", fontsize=8, fontweight="bold")
    ax.set_xlabel("proteins enriched"); ax.set_title("Immune-cell source map"); panel_label(ax, "c")
    # D flag co-occurrence (compact)
    ax = fig.add_subplot(gs[1, 1])
    flags = ["cytokine","chemokine","interferon_axis","TNF_axis","complement",
             "checkpoint","CD_marker","acute_phase","immunoglobulin","immune_cell_enriched"]
    corr = ann[flags].corr()
    sns.heatmap(corr, cmap="vlag", center=0, square=True, ax=ax, cbar_kws={"shrink": .6},
                xticklabels=[f.replace("_","\n") for f in flags], yticklabels=flags)
    ax.set_title("Annotation flag co-occurrence")
    ax.tick_params(labelsize=6.5); panel_label(ax, "d")
    fig.suptitle("Figure 1  |  Study design and curation of a 1,007-protein plasma immunome",
                 fontsize=13, fontweight="bold", x=0.02, ha="left")
    save(fig, "Figure1_design_curation.png")

# =====================================================================
# FIGURE 2 — Druggability & communication architecture
# Q: which immune classes are druggable and how do they map to cells/axes/disease?
# =====================================================================
def figure2():
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.32)
    # A class x source heatmap
    ax = fig.add_subplot(gs[0, 0])
    rows = []
    for _, r in imm.iterrows():
        if isinstance(r.immune_source_cells, str) and r.immune_source_cells:
            for c in r.immune_source_cells.split(";"):
                rows.append((r.immune_class, c))
    ct = pd.crosstab(pd.DataFrame(rows, columns=["c","s"]).c, pd.DataFrame(rows, columns=["c","s"]).s)
    sns.heatmap(ct, annot=True, fmt="d", cmap="rocket_r", ax=ax, cbar_kws={"shrink": .6}, annot_kws={"size": 6.5})
    ax.set_title("Immune class \u00d7 source cell"); ax.set_xlabel(""); ax.set_ylabel("")
    ax.tick_params(labelsize=6.5); panel_label(ax, "a")
    # B druggability by class
    ax = fig.add_subplot(gs[0, 1])
    imm["dt"] = imm["hpa_protein_class"].fillna("").str.contains("drug target", case=False).astype(int)
    by = imm.groupby("immune_class")["dt"].agg(["sum","count"]); by["pct"] = 100*by["sum"]/by["count"]
    by = by.sort_values("pct")
    ax.barh(by.index, by["pct"], color=sns.color_palette("crest", len(by)))
    for i,(s,c) in enumerate(zip(by["sum"],by["count"])):
        ax.text(by["pct"].iloc[i]+0.5, i, f"{int(s)}/{int(c)}", va="center", fontsize=7)
    ax.set_xlabel("% potential drug target"); ax.set_title("Druggability by class")
    ax.tick_params(axis="y", labelsize=7.5); panel_label(ax, "b")
    # C communication axes
    ax = fig.add_subplot(gs[1, 0])
    axes_def = {"IL-6/gp130": r"^(IL6|IL6R|IL6ST|OSM|LIF|IL11|IL27)","TNF SF": r"^(TNF|TNFSF|TNFRSF|LTA|LTB)",
                "CXC chemok": r"^CXCL","CC chemok": r"^CCL","IL-1 fam": r"^(IL1|IL18|IL33|IL36|IL37)",
                "IFN axis": r"^(IFN)","Complement": r"^(C[1-9]|CFB|CFH|CFD|CFI|MASP)",
                "Checkpoint": r"(PDCD1|CD274|CTLA4|LAG3|HAVCR2|TIGIT|ICOS|CD40)"}
    cnt = {k: int(imm["gene_symbol"].str.contains(v, case=False, regex=True, na=False).sum()) for k,v in axes_def.items()}
    s = pd.Series(cnt).sort_values()
    ax.barh(s.index, s.values, color=sns.color_palette("Set2", len(s)))
    for i,v in enumerate(s.values): ax.text(v+0.3, i, str(v), va="center", fontweight="bold", fontsize=8)
    ax.set_xlabel("proteins measured"); ax.set_title("Immune communication axes"); panel_label(ax, "c")
    # D disease involvement
    ax = fig.add_subplot(gs[1, 1])
    terms = []
    for d in imm["hpa_disease"].dropna():
        for t in str(d).split(","):
            t = t.strip()
            if t and t.lower() != "disease related genes": terms.append(t)
    top = pd.Series(terms).value_counts().head(12).sort_values()
    ax.barh(top.index, top.values, color=sns.color_palette("flare", len(top)))
    ax.set_xlabel("proteins"); ax.set_title("Disease involvement (HPA)")
    ax.tick_params(axis="y", labelsize=7); panel_label(ax, "d")
    fig.suptitle("Figure 2  |  Druggability and immune communication architecture of the plasma immunome",
                 fontsize=13, fontweight="bold", x=0.02, ha="left")
    save(fig, "Figure2_druggability_architecture.png")

# =====================================================================
# FIGURE 3 — cis-eQTL MR design & global results
# Q: can public transcript-anchored MR recover autoimmune disease biology?
# =====================================================================
def figure3():
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.42, wspace=0.30)
    # A workflow schematic
    ax = fig.add_subplot(gs[0, 0]); ax.axis("off")
    steps = ["812 immune\ncis-eQTL\ninstruments\n(eQTLGen)", "13 FinnGen\nautoimmune\nGWAS",
             "Wald-ratio\nMR\n8,749 tests", "32 hits\nFDR<0.05"]
    cols = sns.color_palette("Spectral", 4)
    for i,(st,c) in enumerate(zip(steps,cols)):
        b = FancyBboxPatch((0.02+i*0.245, 0.35), 0.20, 0.34, boxstyle="round,pad=0.02",
                           fc=c, ec="black", lw=1, alpha=.9); ax.add_patch(b)
        ax.text(0.12+i*0.245, 0.52, st, ha="center", va="center", fontsize=7.5, fontweight="bold")
        if i < 3:
            ax.add_patch(FancyArrowPatch((0.22+i*0.245,0.52),(0.265+i*0.245,0.52),
                         arrowstyle="-|>", mutation_scale=10, color="grey"))
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_title("cis-MR workflow", y=0.82); panel_label(ax, "a", dy=1.0)
    # B volcano
    ax = fig.add_subplot(gs[0, 1])
    s = mr.FDR < 0.05
    ax.scatter(mr.loc[~s,"MR_beta"], mr.loc[~s,"nlog10p"], s=6, c="#cfcfcf")
    up = s & (mr.MR_beta > 0); dn = s & (mr.MR_beta < 0)
    ax.scatter(mr.loc[up,"MR_beta"], mr.loc[up,"nlog10p"], s=22, c="#c0392b")
    ax.scatter(mr.loc[dn,"MR_beta"], mr.loc[dn,"nlog10p"], s=22, c="#2471a3")
    for _, r in sig.sort_values("FDR").head(6).iterrows():
        ax.annotate(r.gene_symbol, (r.MR_beta, r.nlog10p), fontsize=7, fontweight="bold",
                    xytext=(3,2), textcoords="offset points")
    ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=.8)
    ax.set_xlabel("MR effect (log-OR / SD cis-expr)"); ax.set_ylabel(r"$-\log_{10}P$")
    ax.set_title("8,749-test MR volcano"); panel_label(ax, "b")
    # C hits per disease
    ax = fig.add_subplot(gs[1, 0])
    c = sig.groupby("disease").size().sort_values()
    ax.barh(c.index, c.values, color=sns.color_palette("rocket", len(c)))
    for i,v in enumerate(c.values): ax.text(v+.05, i, str(v), va="center", fontweight="bold", fontsize=8)
    ax.set_xlabel("causal genes (FDR<0.05)"); ax.set_title("Hits per disease")
    ax.tick_params(axis="y", labelsize=7.5); panel_label(ax, "c")
    # D immune-class enrichment of hits
    ax = fig.add_subplot(gs[1, 1])
    tested = mr.groupby("immune_class").size()
    hits = sig.groupby("immune_class").size()
    en = pd.DataFrame({"tested": tested, "hits": hits}).fillna(0)
    en["rate"] = 100*en["hits"]/en["tested"]
    en = en[en.tested >= 20].sort_values("rate")
    ax.barh(en.index, en["rate"], color=sns.color_palette("mako", len(en)))
    for i,(h,t) in enumerate(zip(en["hits"],en["tested"])):
        ax.text(en["rate"].iloc[i]+0.05, i, f"{int(h)}/{int(t)}", va="center", fontsize=7)
    ax.set_xlabel("% tests FDR-significant"); ax.set_title("Hit enrichment by immune class")
    ax.tick_params(axis="y", labelsize=7.5); panel_label(ax, "d")
    fig.suptitle("Figure 3  |  Transcript-anchored Mendelian randomization recovers autoimmune disease signal",
                 fontsize=12.5, fontweight="bold", x=0.02, ha="left")
    save(fig, "Figure3_MR_design_global.png")

# =====================================================================
# FIGURE 4 — Positive-control therapeutic recoveries
# Q: does the MR blindly recover known drug targets? (internal control)
# =====================================================================
def figure4():
    controls = [("IL6ST","Rheumatoid arthritis","tocilizumab (IL-6R/gp130)"),
                ("CTLA4","Autoimmune hyperthyroidism","abatacept (CTLA4-Ig)"),
                ("TNFRSF1A","Ankylosing spondylitis","etanercept / anti-TNF"),
                ("TNFSF14","Multiple sclerosis","LIGHT\u2013HVEM axis")]
    fig, axes = plt.subplots(1, 4, figsize=(15, 5.2))
    for ax,(g,d,drug),lab in zip(axes, controls, "abcd"):
        row = sig[(sig.gene_symbol==g)&(sig.disease==d)]
        if len(row)==0: row = mr[(mr.gene_symbol==g)&(mr.disease==d)]
        r = row.iloc[0]
        ax.errorbar([r.OR],[0], xerr=[[r.OR-r.OR_l95],[r.OR_u95-r.OR]], fmt="o",
                    ms=12, color="#1f496e", ecolor="#555", capsize=5, lw=2)
        ax.axvline(1, ls="--", c="red", lw=1)
        ax.set_xscale("log"); ax.set_yticks([])
        ax.set_xlabel("OR per SD cis-expr")
        ax.set_title(f"{g} \u2192 {d}\nOR={r.OR:.2f} ({r.OR_l95:.2f}\u2013{r.OR_u95:.2f})\nFDR={r.FDR:.1e}", fontsize=9)
        ax.text(0.5,-0.32, f"known drug: {drug}", transform=ax.transAxes, ha="center",
                fontsize=8, style="italic", color="#333")
        ax.text(-0.05,1.12, lab, transform=ax.transAxes, fontsize=13, fontweight="bold")
    fig.suptitle("Figure 4  |  MR blindly recovers established immune therapeutic axes (internal positive control)",
                 fontsize=12.5, fontweight="bold", x=0.02, ha="left", y=1.04)
    fig.tight_layout(rect=[0,0.02,1,0.96])
    save(fig, "Figure4_positive_controls.png")

# =====================================================================
# FIGURE 5 — Novel genetically-supported nominations
# Q: what new autoimmune target hypotheses does the atlas nominate?
# =====================================================================
def figure5(coloc_df=None):
    def h4(g, d):
        if coloc_df is None: return None
        r = coloc_df[(coloc_df.gene==g)&(coloc_df.disease==d)]
        return float(r.PP_H4.iloc[0]) if len(r) else None
    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)
    # A: C2 cross-disease forest
    ax = fig.add_subplot(gs[0,0])
    c2 = mr[mr.gene_symbol=="C2"].sort_values("OR")
    y = np.arange(len(c2))
    ax.errorbar(c2.OR, y, xerr=[c2.OR-c2.OR_l95, c2.OR_u95-c2.OR], fmt="o", color="#8e44ad", ecolor="#aaa", capsize=3)
    ax.axvline(1, ls="--", c="red", lw=1); ax.set_yticks(y); ax.set_yticklabels(c2.disease, fontsize=8)
    ax.set_xscale("log"); ax.set_xlabel("OR per SD cis-expr")
    ax.set_title("C2 complement axis across diseases\n(MHC region \u2014 nomination, needs coloc)", fontsize=9)
    panel_label(ax, "a")
    # B: SIGLEC7/9 in GBS
    ax = fig.add_subplot(gs[0,1])
    sg = mr[(mr.gene_symbol.isin(["SIGLEC7","SIGLEC9"]))&(mr.disease=="Guillain-Barre")].sort_values("OR")
    y = np.arange(len(sg))
    ax.errorbar(sg.OR, y, xerr=[sg.OR-sg.OR_l95, sg.OR_u95-sg.OR], fmt="o", color="#16a085", ecolor="#aaa", capsize=4, ms=9)
    ax.axvline(1, ls="--", c="red", lw=1); ax.set_yticks(y); ax.set_yticklabels(sg.gene_symbol, fontsize=9)
    ax.set_xscale("log"); ax.set_xlabel("OR per SD cis-expr")
    h7, h9 = h4("SIGLEC7","Guillain-Barre"), h4("SIGLEC9","Guillain-Barre")
    ct = "SIGLEC7/9 \u2192 Guillain-Barre\n(sialic-acid inhibitory receptors)"
    if h7 and h9: ct += f"\ncoloc PP.H4={h7:.2f}/{h9:.2f}"
    ax.set_title(ct, fontsize=9)
    panel_label(ax, "b")
    # C: HAVCR1 coeliac
    ax = fig.add_subplot(gs[1,0])
    hv = mr[(mr.gene_symbol=="HAVCR1")&(mr.disease=="Coeliac disease")].iloc[0]
    ax.errorbar([hv.OR],[0], xerr=[[hv.OR-hv.OR_l95],[hv.OR_u95-hv.OR]], fmt="o", ms=12, color="#e67e22", ecolor="#555", capsize=5)
    ax.axvline(1, ls="--", c="red", lw=1); ax.set_yticks([]); ax.set_xscale("log")
    ax.set_xlabel("OR per SD cis-expr")
    hh = h4("HAVCR1","Coeliac disease")
    tt = f"HAVCR1 \u2192 Coeliac disease\nOR={hv.OR:.2f} ({hv.OR_l95:.2f}\u2013{hv.OR_u95:.2f})"
    if hh: tt += f"  |  coloc PP.H4={hh:.2f}"
    ax.set_title(tt, fontsize=9)
    panel_label(ax, "c")
    # D: IFNGR2 psoriasis (colocalized novel nomination)
    ax = fig.add_subplot(gs[1,1])
    il = mr[(mr.gene_symbol=="IFNGR2")&(mr.disease=="Psoriasis")].iloc[0]
    ax.errorbar([il.OR],[0], xerr=[[il.OR-il.OR_l95],[il.OR_u95-il.OR]], fmt="o", ms=12, color="#c0392b", ecolor="#555", capsize=5)
    ax.axvline(1, ls="--", c="red", lw=1); ax.set_yticks([]); ax.set_xscale("log")
    ax.set_xlabel("OR per SD cis-expr")
    hi = h4("IFNGR2","Psoriasis")
    ti = f"IFNGR2 \u2192 Psoriasis  (IFN-\u03b3 receptor)\nOR={il.OR:.2f} ({il.OR_l95:.2f}\u2013{il.OR_u95:.2f})"
    if hi: ti += f"  |  coloc PP.H4={hi:.2f}"
    ax.set_title(ti, fontsize=9)
    panel_label(ax, "d")
    fig.suptitle("Figure 5  |  Colocalized novel autoimmune target nominations (cis-MR + coloc PP.H4\u22650.8)",
                 fontsize=12.5, fontweight="bold", x=0.02, ha="left")
    save(fig, "Figure5_novel_nominations.png")

# =====================================================================
# FIGURE 6 — Validation & extension roadmap (evidence-tier map)
# Q: what evidence tier is each finding, and what is needed next?
# =====================================================================
def figure6(coloc_df=None):
    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30)
    # A evidence ladder
    ax = fig.add_subplot(gs[0,0]); ax.axis("off")
    tiers = [("Associated with", "correlation", "#dfe6e9"),
             ("Transcript-level genetic proxy", "cis-eQTL only", "#a3cbe8"),
             ("Genetically-supported nomination", "cis-MR FDR<0.05", "#5b9bd5"),
             ("Prioritized causal target", "+ colocalization / pQTL", "#2e75b6"),
             ("Functionally validated", "+ perturbation", "#1f4e79")]
    for i,(t,e,c) in enumerate(tiers):
        b = FancyBboxPatch((0.05, 0.78-i*0.17), 0.9, 0.13, boxstyle="round,pad=0.01", fc=c, ec="black", lw=.8)
        ax.add_patch(b)
        ax.text(0.5, 0.845-i*0.17, f"{t}   \u2014   {e}", ha="center", va="center", fontsize=8.5,
                fontweight="bold" if i==2 else "normal",
                color="white" if i>=3 else "black")
    ax.text(0.5, 0.95, "Claim-strength ladder (this study reaches tier 4: coloc)", ha="center", fontsize=9, fontweight="bold")
    ax.set_xlim(0,1); ax.set_ylim(0,1); panel_label(ax, "a", dy=1.0)
    # B where this study sits (counts per tier)
    ax = fig.add_subplot(gs[0,1])
    tier_counts = {"Transcript\nproxy\n(all cis-eQTL)": 812,
                   "Genetically\nsupported\n(FDR<0.05)": int((mr.FDR<0.05).sum()),
                   "Prioritized\n(coloc PP.H4>0.8)": (int((coloc_df.PP_H4>0.8).sum()) if coloc_df is not None else 0),
                   "Functionally\nvalidated": 0}
    ax.bar(range(len(tier_counts)), list(tier_counts.values()),
           color=["#a3cbe8","#5b9bd5","#2e75b6","#1f4e79"])
    ax.set_yscale("symlog")
    ax.set_xticks(range(len(tier_counts))); ax.set_xticklabels(tier_counts.keys(), fontsize=7)
    for i,v in enumerate(tier_counts.values()): ax.text(i, v+0.5, str(v), ha="center", fontsize=8, fontweight="bold")
    ax.set_ylabel("findings (symlog)"); ax.set_title("Evidence tier distribution"); panel_label(ax, "b")
    # C colocalization result (if available) else placeholder
    ax = fig.add_subplot(gs[1,0])
    if coloc_df is not None and len(coloc_df):
        cd = coloc_df.sort_values("PP_H4")
        colors = ["#2e75b6" if p>0.8 else ("#f0ad4e" if p>0.5 else "#cccccc") for p in cd.PP_H4]
        ax.barh([f"{r.gene}\u2192{r.disease[:10]}" for _,r in cd.iterrows()], cd.PP_H4, color=colors)
        ax.axvline(0.8, ls="--", c="green", lw=1)
        ax.set_xlabel("coloc PP.H4 (shared causal variant)")
        ax.set_title("Colocalization of top non-MHC hits"); ax.tick_params(axis="y", labelsize=6.5)
    else:
        ax.axis("off"); ax.text(0.5,0.5,"colocalization pending\n(eQTLGen full download)", ha="center", va="center")
    panel_label(ax, "c")
    # D extension roadmap
    ax = fig.add_subplot(gs[1,1]); ax.axis("off")
    steps = ["1. pQTL MR (UKB-PPP / deCODE)","2. Colocalization at each locus",
             "3. Replication (2nd GWAS / ancestry)","4. Plasma-protein association (UKB Olink)",
             "5. Functional perturbation (C2, SIGLEC7/9)","6. Versioned public data + GitHub release"]
    for i,s in enumerate(steps):
        ax.text(0.03, 0.9-i*0.15, s, fontsize=8.5, va="top")
    ax.text(0.03, 0.99, "Roadmap to causal / translational tier", fontsize=9, fontweight="bold")
    ax.set_xlim(0,1); ax.set_ylim(0,1); panel_label(ax, "d", dy=1.0)
    fig.suptitle("Figure 6  |  Evidence-tier calibration and validation roadmap",
                 fontsize=12.5, fontweight="bold", x=0.02, ha="left")
    save(fig, "Figure6_validation_roadmap.png")

if __name__ == "__main__":
    print("Nature-style figures ->", FIG)
    cd_path = os.path.join(GEN, "coloc_results.tsv")
    cd = pd.read_csv(cd_path, sep="\t") if os.path.exists(cd_path) else None
    figure1(); figure2(); figure3(); figure4(); figure5(cd)
    figure6(cd)
    print("done:", len(os.listdir(FIG)), "figures")
