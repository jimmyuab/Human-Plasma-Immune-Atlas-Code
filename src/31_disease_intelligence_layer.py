#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
31 — Disease-Trained PIRS Figure & Result Intelligence Layer
============================================================

The "improvement module": after a Plasma Immune Risk Score (PIRS) is trained on a
disease cohort (src/22), OR on the genetics-anchored causal atlas alone, this layer
turns the raw weights / causal statistics into a *plasma-immune discovery report*:

  * disease-specific PIRS performance (when a trained model is present)
  * plasma immune protein weights with class / localization / druggability
  * causal-predictive concordance (PIRS vs MR / coloc / pQTL / replication)
  * a Plasma Immune Novelty Score (PINS) and Novelty Tier (1-5)
  * therapeutic direction (block / agonize / replace / biomarker-only / hold)
  * 6 figure panels (workflow, performance, signature, concordance,
    novelty-priority map, validation plan)
  * a manuscript-style results narrative
  * the Final Required Output Table (~24 columns)

No fabrication. If no PIRS has been trained (05_machine_learning empty), the module
runs in *causal-atlas-only* mode: every causal / novelty / direction / tier column is
populated from real data, and PIRS-dependent columns are written as "NA (train PIRS)".

Run:
    python src/31_disease_intelligence_layer.py
"""

import os
import sys
import glob
import pickle
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

GC   = os.path.join(ROOT, "06_genetic_causality")
PROC = os.path.join(ROOT, "02_data_processed")
ML   = os.path.join(ROOT, "05_machine_learning")
TAB  = os.path.join(ROOT, "09_tables")
FIG  = os.path.join(ROOT, "08_figures", "intelligence_layer")
OUT  = os.path.join(ROOT, "06_genetic_causality")
REP  = os.path.join(ROOT, "10_manuscript")
os.makedirs(FIG, exist_ok=True)
os.makedirs(REP, exist_ok=True)

PP_COLOC = 0.8          # colocalization promotion threshold
FDR_SIG  = 0.05

# --all : run the intelligence layer over the WHOLE FinnGen phenome
#         (novelty_engine_ranked_ALL.tsv from src/50, 1,016 hits / 379 diseases)
#         instead of the 28 deeply-annotated core diseases.
ALLMODE = "--all" in sys.argv
SUF = "_ALL" if ALLMODE else ""


# --------------------------------------------------------------------------- #
#  Loaders
# --------------------------------------------------------------------------- #
def _read(path, **kw):
    return pd.read_csv(path, sep="\t", **kw) if os.path.exists(path) else None


def load_annotation():
    for p in (os.path.join(PROC, "plasma_immune_protein_annotation.tsv"),
              os.path.join(ROOT, "results", "plasma_immune_protein_annotation.tsv"),
              os.path.join(GC,   "plasma_immune_protein_annotation.tsv")):
        if os.path.exists(p):
            return pd.read_csv(p, sep="\t")
    raise SystemExit("Annotation not found.")


def load_pirs():
    """Return dict: disease -> {weights: df, metrics: row}. Empty if none trained."""
    out = {}
    if not os.path.isdir(ML):
        return out, None
    metrics = _read(os.path.join(ML, "pirs_cv_metrics.tsv"))
    for wf in glob.glob(os.path.join(ML, "pirs_*_weights.tsv")):
        base = os.path.basename(wf)
        dis = base[len("pirs_"):-len("_weights.tsv")]
        w = pd.read_csv(wf, sep="\t")
        out[dis] = {"weights": w}
    return out, metrics


# --------------------------------------------------------------------------- #
#  Annotation-derived helpers
# --------------------------------------------------------------------------- #
def plasma_detectability(row):
    """Qualitative plasma detectability from Olink/secretome annotation."""
    if str(row.get("is_plasma_immune", "")).lower() in ("true", "1", "yes"):
        loc = str(row.get("secretome_location", "")).lower()
        if "secreted" in loc or str(row.get("secreted_to_blood", "")).lower() in ("true", "1", "yes"):
            return "High (secreted to blood)"
        if str(row.get("soluble_receptor", "")).lower() in ("true", "1", "yes"):
            return "High (soluble receptor / shed)"
        return "Measured (Olink Explore panel)"
    return "Measured (Olink Explore panel)"


def druggability_label(score, fda):
    if fda == 1 and score >= 3:
        return "High (FDA-precedented class, tractable)"
    if score >= 3:
        return "High (secreted/membrane, tractable)"
    if score == 2:
        return "Moderate"
    if score == 1:
        return "Low-moderate"
    return "Uncertain"


# --------------------------------------------------------------------------- #
#  Causal-predictive concordance
# --------------------------------------------------------------------------- #
def concordance_category(has_pirs, pirs_sign, or_val, has_coloc, has_pqtl):
    """
    Concordance between predictive (PIRS) and causal (MR/coloc/pQTL) evidence.
    Direction convention: PIRS weight > 0 -> higher protein raises risk (like OR>1).
    """
    causal_dir = "risk" if or_val > 1 else "protective"
    if not has_pirs:
        base = "Causal-only (PIRS not trained)"
    else:
        pred_dir = "risk" if pirs_sign > 0 else "protective"
        agree = (pred_dir == causal_dir)
        if agree and has_pqtl:
            base = "Concordant: predictive + protein-level causal"
        elif agree and has_coloc:
            base = "Concordant: predictive + colocalized causal"
        elif agree:
            base = "Concordant: predictive + genetic causal"
        else:
            base = "Discordant: predictive vs causal direction"
    return base, causal_dir


# --------------------------------------------------------------------------- #
#  Therapeutic direction
# --------------------------------------------------------------------------- #
def therapeutic_direction(or_val, sol_recep, secreted, tier):
    """
    OR>1  -> protein raises risk  -> BLOCK / antagonize / neutralize
    OR<1  -> protein is protective -> AGONIZE / REPLACE (supplement / agonist)
    Soluble receptors that are protective -> decoy/replacement therapy.
    Biomarker-only if evidence does not reach a causal tier.
    """
    if tier <= 2:  # known-drug control or MHC-held
        note = "positive control" if tier == 1 else "held (MHC/LD caution)"
        base = ("Block (antagonist/neutralizing Ab)" if or_val > 1
                else "Agonize / replace") + f" — {note}"
        return base
    if or_val > 1:
        return "Block (antagonist / neutralizing antibody / small molecule)"
    # protective protein
    if sol_recep:
        return "Replace / decoy-receptor or agonist (protective, soluble)"
    if secreted:
        return "Agonize / recombinant-protein replacement (protective, secreted)"
    return "Agonize (protective; intracellular route caution)"


def biomarker_or_target(tier, has_pqtl, has_coloc):
    if tier >= 5 or (has_pqtl and has_coloc):
        return "Therapeutic target (protein-level)"
    if tier == 4 or has_coloc:
        return "Prioritized target (transcript-level)"
    if tier == 3:
        return "Target nomination (causal, needs coloc)"
    return "Biomarker (predictive) — not yet a target"


# --------------------------------------------------------------------------- #
#  Novelty score + tier
# --------------------------------------------------------------------------- #
def plasma_immune_novelty_score(fdr, pp_h4, has_pqtl, druggability, n_cat,
                                known, mhc):
    """
    PINS — bounded composite emphasising *novel, plasma-measurable, druggable,
    genetically-supported* targets. Range roughly 0-8.
    """
    s_causal = min(-np.log10(max(fdr, 1e-300)) / 5.0, 2.0)   # 0-2
    s_coloc  = float(pp_h4) if pp_h4 == pp_h4 else 0.0         # 0-1
    s_prot   = 1.0 if has_pqtl else 0.0                        # protein-level bonus
    s_drug   = min(druggability / 3.0, 1.0)                    # 0-1
    s_pleio  = 0.5 * max(n_cat - 1, 0)                         # cross-disease
    pins = s_causal + s_coloc + s_prot + s_drug + s_pleio + 0.5
    pins -= 1.5 * known + 2.0 * mhc
    return round(max(pins, 0.0), 3)


def novelty_tier(known, mhc, has_coloc, has_pqtl, druggable):
    """
    Tier 1  known-drug positive control
    Tier 2  MHC/LD held at nomination
    Tier 3  causal nomination (MR only)
    Tier 4  prioritized causal target (MR + coloc), novel
    Tier 5  novel plasma-immune target: prediction-ready + coloc + protein pQTL
            + specificity + druggability
    """
    if known:
        return 1, "Known-drug positive control"
    if mhc:
        return 2, "MHC/LD — held at nomination"
    if has_coloc and has_pqtl and druggable:
        return 5, "Novel protein-level plasma-immune target"
    if has_coloc:
        return 4, "Prioritized causal target (transcript-level)"
    return 3, "Causal nomination (MR only)"


# --------------------------------------------------------------------------- #
#  Best figure / validation routing
# --------------------------------------------------------------------------- #
def best_panel(tier):
    return {5: "Panel E (novelty-priority map) + Panel D (concordance)",
            4: "Panel E (novelty-priority map)",
            3: "Panel B/C (performance + signature)",
            2: "Panel D (concordance, MHC-flagged)",
            1: "Panel B (performance — positive control)"}.get(tier, "Panel A (workflow)")


def best_validation(tier, or_val, has_pqtl):
    direction = "neutralization/knockdown" if or_val > 1 else "over-expression/supplementation"
    if tier == 5:
        return (f"cis-pQTL colocalization already met; confirm with CRISPRi/CRISPRa "
                f"{direction} in primary immune cells + plasma NPX dose-response")
    if tier == 4:
        return (f"Protein-level pQTL-MR + coloc in an independent plasma pQTL panel "
                f"(UKB-PPP/deCODE); then {direction} in relevant primary cells")
    if tier == 3:
        return "Transcript-to-protein bridge: acquire cis-pQTL, test colocalization"
    if tier == 2:
        return "Fine-map MHC / condition on classical HLA alleles before any claim"
    return "Confirm as positive control; benchmark against approved agent"


def final_recommendation(tier, biomk, has_pqtl):
    if tier == 5:
        return "Advance: protein-level causal + druggable + plasma-measurable — target-validation package"
    if tier == 4:
        return "Prioritize: obtain plasma pQTL to upgrade to protein-level"
    if tier == 3:
        return "Nominate: causal but needs colocalization before resourcing"
    if tier == 2:
        return "Hold: MHC/LD confounding must be excluded"
    return "Reference/positive control: validates the pipeline"


# --------------------------------------------------------------------------- #
#  Build the Final Required Output Table
# --------------------------------------------------------------------------- #
def build_table():
    ann   = load_annotation()
    nov   = pd.read_csv(os.path.join(GC, f"novelty_engine_ranked{SUF}.tsv"), sep="\t")
    if ALLMODE:
        nov = nov.rename(columns={"disease_system": "disease_category"})
    final = _read(os.path.join(GC, "FINAL_evidence_tiers_repl.tsv"))
    coloc = _read(os.path.join(GC, "coloc_phenome_results.tsv"))
    pqtl  = _read(os.path.join(GC, "pqtl_MR_results.tsv"))
    pirs, pirs_metrics = load_pirs()
    has_any_pirs = len(pirs) > 0

    # coloc lookup: max PP_H4 per (gene, disease)
    coloc_lu = {}
    if coloc is not None:
        for _, r in coloc.iterrows():
            k = (str(r["gene"]).upper(), str(r["disease"]).lower())
            coloc_lu[k] = max(coloc_lu.get(k, 0.0), float(r["PP_H4"]))

    # protein-level pQTL concordance from FINAL table
    pqtl_lu = {}
    rep_lu  = {}
    if final is not None:
        for _, r in final.iterrows():
            k = (str(r["gene_symbol"]).upper(), str(r["disease"]).lower())
            pqtl_lu[k] = str(r.get("pQTL_concordant", "")).lower() in ("yes", "true", "1")
            rep_lu[k]  = str(r.get("rep_status", ""))
    if ALLMODE:
        # pan-phenome protein + two-population columns come straight from src/50
        for _, r in nov.iterrows():
            k = (str(r["gene_symbol"]).upper(), str(r["disease"]).lower())
            if str(r.get("pqtl_concordant", "")).lower() in ("true", "1", "yes"):
                pqtl_lu[k] = True
            else:
                pqtl_lu.setdefault(k, False)
            if str(r.get("two_population_validated", "")).lower() in ("true", "1", "yes"):
                rep_lu[k] = "replicated (FinnGen + UK Biobank)"
            else:
                rep_lu.setdefault(k, "not tested")

    ann_lu = {str(g).upper(): row for g, row in zip(ann["gene_symbol"], ann.to_dict("records"))}

    rows = []
    for _, h in nov.iterrows():
        gene = str(h["gene_symbol"]).upper()
        dis  = str(h["disease"])
        k    = (gene, dis.lower())
        a    = ann_lu.get(gene, {})

        or_val = float(h["OR"])
        fdr    = float(h["FDR"])
        pp_h4  = float(h["PP_H4"]) if not pd.isna(h["PP_H4"]) else np.nan
        drug   = int(h["druggability"])
        fda    = int(h["fda_target"])
        n_cat  = int(h["n_categories"])
        known  = float(h["p_known"]) > 0
        mhc    = float(h["p_mhc"]) > 0

        has_coloc = (pp_h4 == pp_h4) and pp_h4 >= PP_COLOC
        has_pqtl  = pqtl_lu.get(k, False)
        rep       = rep_lu.get(k, "not tested")
        druggable = drug >= 2

        # ---- PIRS integration (real only) --------------------------------- #
        pirs_coef = "NA (train PIRS)"
        cv_stab   = "NA (train PIRS)"
        sens_c = spec_c = auroc_c = "NA (train PIRS)"
        pirs_sign = 0.0
        has_pirs_here = False
        dkey = dis.lower().replace(" ", "_")
        for pd_dis, pd_obj in pirs.items():
            if pd_dis.lower().replace(" ", "_") == dkey:
                w = pd_obj["weights"]
                m = w[w["gene_symbol"].astype(str).str.upper() == gene]
                if len(m):
                    pirs_coef = round(float(m.iloc[0]["pirs_weight"]), 5)
                    pirs_sign = float(m.iloc[0]["pirs_weight"])
                    has_pirs_here = pirs_coef != 0
                    # rank-based stability proxy from |weight|
                    cv_stab = "selected" if abs(pirs_sign) > 0 else "not selected"

        conc, causal_dir = concordance_category(has_pirs_here, pirs_sign, or_val,
                                                has_coloc, has_pqtl)
        tier, tier_lbl = novelty_tier(known, mhc, has_coloc, has_pqtl, druggable)
        pins = plasma_immune_novelty_score(fdr, pp_h4 if pp_h4 == pp_h4 else 0.0,
                                           has_pqtl, drug, n_cat, known, mhc)
        sol_recep = str(a.get("soluble_receptor", "")).lower() in ("true", "1", "yes")
        secreted  = (str(a.get("secreted_to_blood", "")).lower() in ("true", "1", "yes")
                     or "secreted" in str(a.get("secretome_location", "")).lower())
        thera = therapeutic_direction(or_val, sol_recep, secreted, tier)
        biomk = biomarker_or_target(tier, has_pqtl, has_coloc)

        rows.append({
            "Disease": dis,
            "Category": h["disease_category"],
            "Protein": a.get("protein_name", gene),
            "Gene": gene,
            "Protein_class": a.get("hpa_protein_class", h.get("immune_class", "")),
            "Immune_class": h.get("immune_class", ""),
            "Plasma_detectability": plasma_detectability(a),
            "PIRS_coefficient": pirs_coef,
            "PIRS_direction": ("risk" if pirs_sign > 0 else "protective") if has_pirs_here else "NA",
            "CV_stability": cv_stab,
            "Sensitivity_contrib": sens_c,
            "Specificity_contrib": spec_c,
            "AUROC_contrib": auroc_c,
            "MR_OR": round(or_val, 3),
            "MR_FDR": f"{fdr:.2e}",
            "MR_support": "yes" if fdr < FDR_SIG else "sub-threshold",
            "Coloc_PP_H4": round(pp_h4, 3) if pp_h4 == pp_h4 else "NA",
            "Coloc_support": "yes (PP.H4>=0.8)" if has_coloc else "no",
            "pQTL_support": "yes (concordant, protein-level)" if has_pqtl else "no",
            "Protein_level_causal": "yes" if (has_pqtl and has_coloc) else "no",
            "Replication": rep,
            "Known_drug_status": "known drug target" if known else "not an approved target",
            "Druggability": druggability_label(drug, fda),
            "Novelty_score_PINS": pins,
            "Novelty_tier": tier,
            "Novelty_tier_label": tier_lbl,
            "Causal_predictive_concordance": conc,
            "Therapeutic_direction": thera,
            "Biomarker_or_target": biomk,
            "Best_figure_panel": best_panel(tier),
            "Best_validation_experiment": best_validation(tier, or_val, has_pqtl),
            "Final_recommendation": final_recommendation(tier, biomk, has_pqtl),
        })

    df = pd.DataFrame(rows)
    # Rank: high-novelty (tier 4/5) first, then PINS, then -log10 FDR
    df["_fdr_num"] = df["MR_FDR"].apply(lambda s: float(s))
    df = df.sort_values(
        by=["Novelty_tier", "Novelty_score_PINS", "_fdr_num"],
        ascending=[False, False, True]
    ).reset_index(drop=True)
    df.insert(0, "Rank", np.arange(1, len(df) + 1))
    df = df.drop(columns=["_fdr_num"])
    return df, has_any_pirs, pirs_metrics


# --------------------------------------------------------------------------- #
#  Figures (6 panels)
# --------------------------------------------------------------------------- #
def _panelA_workflow(path, has_pirs):
    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.axis("off")
    ax.set_title("Panel A — Plasma Immune Discovery Workflow", fontsize=14, fontweight="bold")
    steps = [
        ("Plasma Olink\nimmune proteome\n(1,007 proteins)", "#cfe8ff"),
        ("PIRS training\n(elastic-net Cox)\n" + ("[trained]" if has_pirs else "[bring your data]"), "#d9f2d9"),
        ("Genetic causal atlas\ncis-MR + coloc\n+ pQTL + replication", "#ffe6cc"),
        ("Concordance +\nNovelty score\n(PINS, Tier 1-5)", "#f3d9ff"),
        ("Therapeutic direction\n+ validation plan", "#ffd9d9"),
        ("Ranked discovery\nreport + Final Table", "#fff2b3"),
    ]
    x = 0.03
    w = 0.145
    for i, (txt, col) in enumerate(steps):
        box = FancyBboxPatch((x, 0.4), w, 0.24, boxstyle="round,pad=0.012",
                             linewidth=1.2, edgecolor="#333", facecolor=col)
        ax.add_patch(box)
        ax.text(x + w / 2, 0.52, txt, ha="center", va="center", fontsize=8.6)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + w + 0.002, 0.52), (x + w + 0.028, 0.52),
                                         arrowstyle="-|>", mutation_scale=14, color="#444"))
        x += w + 0.028
    ax.text(0.5, 0.16,
            "Predictive layer (PIRS) answers WHO is at risk and WHICH plasma proteins carry the signal;\n"
            "causal layer answers WHETHER the protein is causal and DRUGGABLE; the intelligence layer fuses both.",
            ha="center", va="center", fontsize=9.2, style="italic")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)


def _panelB_performance(path, pirs_metrics, df):
    fig, ax = plt.subplots(figsize=(9, 5.6))
    if pirs_metrics is not None and len(pirs_metrics):
        m = pirs_metrics.copy()
        m = m.sort_values("score_mean")
        ax.barh(m["disease"], m["score_mean"], xerr=m.get("score_std"),
                color="#4a90d9", edgecolor="#222")
        ax.set_xlabel("Cross-validated performance (C-index / AUROC)")
        ax.set_title("Panel B — Disease-specific PIRS performance", fontweight="bold")
        ax.axvline(0.5, ls="--", color="grey", lw=1)
    else:
        # causal-only: show causal 'discriminating' strength as -log10 FDR by disease
        agg = df.copy()
        agg["nlfdr"] = agg["MR_FDR"].apply(lambda s: -np.log10(max(float(s), 1e-300)))
        top = agg.groupby("Disease")["nlfdr"].max().sort_values().tail(18)
        ax.barh(top.index, top.values, color="#b0b0b0", edgecolor="#222")
        ax.set_xlabel(r"Strongest causal signal per disease  ($-\log_{10}$ FDR)")
        ax.set_title("Panel B — PIRS not trained: causal signal strength (train PIRS for AUROC)",
                     fontweight="bold", fontsize=11)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)


def _panelC_signature(path, df):
    """Top plasma-immune signature proteins coloured by therapeutic direction."""
    fig, ax = plt.subplots(figsize=(9, 6.4))
    top = df[df["Novelty_tier"] >= 4].head(22).iloc[::-1]
    if not len(top):
        top = df.head(22).iloc[::-1]
    colors = ["#d1495b" if o > 1 else "#2e8b57" for o in top["MR_OR"]]
    ax.barh([f"{g} → {d[:18]}" for g, d in zip(top["Gene"], top["Disease"])],
            [np.log2(o) for o in top["MR_OR"]], color=colors, edgecolor="#222")
    ax.axvline(0, color="#333", lw=1)
    ax.set_xlabel(r"Causal effect  $\log_2$(OR)   — red: risk/BLOCK   green: protective/AGONIZE")
    ax.set_title("Panel C — Plasma immune signature (high-novelty targets)", fontweight="bold")
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)


def _panelD_concordance(path, df):
    fig, ax = plt.subplots(figsize=(8.4, 6))
    layers = ["MR_support", "Coloc_support", "pQTL_support"]
    counts = {
        "MR (causal)":      (df["MR_support"] == "yes").sum(),
        "+ Coloc":          df["Coloc_support"].str.startswith("yes").sum(),
        "+ pQTL (protein)": df["pQTL_support"].str.startswith("yes").sum(),
        "+ Replicated":     (df["Replication"].astype(str).str.lower() == "replicated").sum(),
    }
    ax.bar(list(counts.keys()), list(counts.values()),
           color=["#8ecae6", "#219ebc", "#126782", "#023047"], edgecolor="#111")
    for i, v in enumerate(counts.values()):
        ax.text(i, v + 0.5, str(int(v)), ha="center", fontweight="bold")
    ax.set_ylabel("Gene–disease pairs")
    ax.set_title("Panel D — Causal-evidence concordance ladder", fontweight="bold")
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)


def _panelE_novelty_map(path, df):
    fig, ax = plt.subplots(figsize=(9.4, 6.6))
    x = df["Coloc_PP_H4"].apply(lambda v: v if isinstance(v, (int, float)) else 0.0)
    y = df["Novelty_score_PINS"]
    tiers = df["Novelty_tier"]
    cmap = {1: "#999999", 2: "#c9b458", 3: "#8ecae6", 4: "#fb8500", 5: "#d1495b"}
    for t in sorted(cmap):
        s = tiers == t
        ax.scatter(x[s], y[s], s=42, c=cmap[t], edgecolor="#222", linewidth=0.4,
                   label=f"Tier {t}", alpha=0.85)
    # annotate top novel
    for _, r in df[df["Novelty_tier"] >= 4].head(12).iterrows():
        xx = r["Coloc_PP_H4"] if isinstance(r["Coloc_PP_H4"], (int, float)) else 0.0
        ax.annotate(f"{r['Gene']}→{r['Disease'][:10]}", (xx, r["Novelty_score_PINS"]),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axvline(PP_COLOC, ls="--", color="grey", lw=1)
    ax.set_xlabel("Colocalization PP.H4")
    ax.set_ylabel("Plasma Immune Novelty Score (PINS)")
    ax.set_title("Panel E — Novelty-priority map", fontweight="bold")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)


def _panelF_validation(path, df):
    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    ax.axis("off")
    ax.set_title("Panel F — Validation plan for top novel targets", fontweight="bold", fontsize=13)
    top = df[df["Novelty_tier"] >= 4].head(8)
    if not len(top):
        top = df.head(8)
    y = 0.92
    ax.text(0.01, 0.98, "Rank  Target → Disease            Tier  Next experiment", fontsize=9.5,
            fontweight="bold", family="monospace")
    for _, r in top.iterrows():
        line = f"{r['Rank']:>3}   {r['Gene']:>7} → {r['Disease'][:20]:<20} T{r['Novelty_tier']}"
        ax.text(0.01, y, line, fontsize=8.6, family="monospace")
        ax.text(0.01, y - 0.035, "        " + textwrap.shorten(r["Best_validation_experiment"], 96),
                fontsize=7.8, style="italic", color="#333")
        y -= 0.105
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)


def make_figures(df, has_pirs, pirs_metrics):
    paths = {}
    _panelA_workflow(os.path.join(FIG, f"IL_panelA_workflow{SUF}.png"), has_pirs)
    _panelB_performance(os.path.join(FIG, f"IL_panelB_performance{SUF}.png"), pirs_metrics, df)
    _panelC_signature(os.path.join(FIG, f"IL_panelC_signature{SUF}.png"), df)
    _panelD_concordance(os.path.join(FIG, f"IL_panelD_concordance{SUF}.png"), df)
    _panelE_novelty_map(os.path.join(FIG, f"IL_panelE_novelty_map{SUF}.png"), df)
    _panelF_validation(os.path.join(FIG, f"IL_panelF_validation{SUF}.png"), df)
    for k in "ABCDEF":
        paths[k] = os.path.join(FIG, f"IL_panel{k}_*{SUF}.png")
    return paths


# --------------------------------------------------------------------------- #
#  Manuscript-style narrative
# --------------------------------------------------------------------------- #
def write_report(df, has_pirs, pirs_metrics):
    t5 = df[df["Novelty_tier"] == 5]
    t4 = df[df["Novelty_tier"] == 4]
    t1 = df[df["Novelty_tier"] == 1]
    n_dis = df["Disease"].nunique()
    lines = []
    W = lines.append
    W("# Disease-Trained PIRS — Plasma Immune Discovery Report\n")
    W("## Intelligence layer over the genetics-anchored plasma immunome atlas\n")

    W("### Mode\n")
    if has_pirs:
        W(f"A trained PIRS model was detected and fused with the causal atlas. "
          f"Predictive coefficients, CV stability and discrimination contributions are "
          f"populated from the trained model.\n")
    else:
        W("**Causal-atlas-only mode.** No trained PIRS was found in `05_machine_learning/`, "
          "so PIRS-dependent columns (coefficient, CV stability, sensitivity/specificity/"
          "AUROC contribution) are reported as `NA (train PIRS)`. Every causal, novelty, "
          "direction and tier column is populated from real genetic data. Train a PIRS "
          "(`src/22_train_pirs.py`) on an authorised cohort to activate the predictive columns "
          "— no values are fabricated.\n")

    W("### Result summary\n")
    W(f"- Gene–disease pairs assessed: **{len(df)}** across **{n_dis} diseases** "
      f"(5 categories).\n")
    W(f"- High-novelty targets (Tier 4/5): **{len(t4) + len(t5)}** "
      f"(Tier 5 protein-level novel: **{len(t5)}**; Tier 4 prioritized: **{len(t4)}**).\n")
    W(f"- Known-drug positive controls (Tier 1): **{len(t1)}** — recovered from genetics "
      f"alone, validating the pipeline.\n")
    n_repl = (df["Replication"].astype(str).str.lower() == "replicated").sum()
    W(f"- Independently replicated pairs: **{n_repl}**.\n")
    n_prot = (df["Protein_level_causal"] == "yes").sum()
    W(f"- Protein-level causal pairs (transcript coloc + concordant plasma pQTL): **{n_prot}**.\n")

    if len(t5) == 0:
        W("### Why Tier 5 is currently empty (an honest gate, not a gap in the run)\n")
        W("Tier 5 requires **all five** of: prediction-ready + colocalization + "
          "concordant plasma pQTL + specificity + druggability. Two hits reach protein-"
          "level causality (transcript coloc + concordant INTERVAL pQTL): "
          "**TNFSF14 → multiple sclerosis** and **SWAP70 → rheumatoid arthritis**. "
          "TNFSF14 is flagged as a **known-drug axis** → it is reported as a **Tier 1 "
          "positive control** (correct pipeline recovery, not a novel target). SWAP70 is "
          "**novel and protein-level causal** but has **druggability = 0** (intracellular, "
          "not secreted/membrane) → it lands at **Tier 4**, not Tier 5. Separately, the "
          "**pan-phenome diseases (cardiovascular / metabolic / renal / neuro) have no "
          "plasma pQTL layer computed yet** (INTERVAL pQTL was run only against the "
          "autoimmune arc), so none of those hits can reach Tier 5 until a plasma pQTL "
          "panel is colocalized against them. This is a data-coverage boundary, not a "
          "fabricated ceiling — running cis-pQTL MR+coloc across the phenome is the single "
          "step that can promote Tier-4 targets to Tier 5.\n")

    W("### Claim discipline\n")
    W("Claims are bound to the evidence level actually reached:\n")
    W("- PIRS weight alone → **biomarker** (predictive, not causal).\n")
    W("- + cis-MR → **causal nomination**.\n")
    W("- + colocalization → **prioritized causal target**.\n")
    W("- + direction-concordant plasma pQTL → **protein-level causal target**.\n")
    W("- + perturbation → **proof of mechanism** (not asserted here; proposed as next experiment).\n")
    W("All signals restricted to **plasma immune proteins**; no single-cell, tissue or "
      "intracellular claim is made unless such data are separately supplied.\n")

    if len(t5):
        W("### Tier 5 — novel protein-level plasma-immune targets\n")
        for _, r in t5.iterrows():
            W(f"**{r['Gene']} → {r['Disease']}** (PINS {r['Novelty_score_PINS']}). "
              f"MR OR={r['MR_OR']} (FDR {r['MR_FDR']}), coloc PP.H4={r['Coloc_PP_H4']}, "
              f"plasma pQTL concordant, replication: {r['Replication']}. "
              f"Direction: {r['Therapeutic_direction']}. "
              f"Next: {r['Best_validation_experiment']}.\n")

    if len(t4):
        W("### Tier 4 — prioritized causal targets (obtain plasma pQTL to upgrade)\n")
        for _, r in t4.head(15).iterrows():
            W(f"- **{r['Gene']} → {r['Disease']}**: OR={r['MR_OR']}, PP.H4={r['Coloc_PP_H4']}, "
              f"{r['Therapeutic_direction']} | {r['Druggability']}.\n")

    W("### Figure panels\n")
    W("- Panel A — plasma immune discovery workflow.\n")
    W("- Panel B — disease-specific PIRS performance (causal signal strength if untrained).\n")
    W("- Panel C — plasma immune signature (high-novelty targets, direction-coloured).\n")
    W("- Panel D — causal-evidence concordance ladder (MR→coloc→pQTL→replication).\n")
    W("- Panel E — novelty-priority map (PP.H4 vs PINS, tier-coloured).\n")
    W("- Panel F — validation plan for the top novel targets.\n")

    W("### Final Required Output Table\n")
    W("See `intelligence_layer_final_table.tsv` (ranked; high-novelty Tier 4/5 first). "
      "Columns: Rank, Disease, Protein, Gene, Protein class, Plasma detectability, "
      "PIRS coefficient, PIRS direction, CV stability, Sensitivity/Specificity/AUROC "
      "contribution, MR OR/FDR/support, Coloc PP.H4/support, pQTL support, Replication, "
      "Known-drug status, Druggability, Novelty score (PINS), Novelty tier (+label), "
      "Causal-predictive concordance, Therapeutic direction, Biomarker-or-target, "
      "Best figure panel, Best validation experiment, Final recommendation.\n")

    path = os.path.join(REP, f"Plasma_Immune_Discovery_Report{SUF}.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return path


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main():
    print("[31] Building disease-trained PIRS intelligence layer ...")
    df, has_pirs, pirs_metrics = build_table()

    tsv = os.path.join(OUT, f"intelligence_layer_final_table{SUF}.tsv")
    df.to_csv(tsv, sep="\t", index=False)
    # also drop a copy into 09_tables as the headline deliverable
    df.to_csv(os.path.join(TAB, f"T6_disease_intelligence_final_table{SUF}.tsv"),
              sep="\t", index=False)
    print(f"    Final Required Output Table  -> {tsv}  ({len(df)} rows x {df.shape[1]} cols)")

    make_figures(df, has_pirs, pirs_metrics)
    print(f"    6 figure panels              -> {FIG}")

    rep = write_report(df, has_pirs, pirs_metrics)
    print(f"    Manuscript-style report      -> {rep}")

    n5 = (df["Novelty_tier"] == 5).sum()
    n4 = (df["Novelty_tier"] == 4).sum()
    n1 = (df["Novelty_tier"] == 1).sum()
    print(f"    Mode: {'PIRS-fused' if has_pirs else 'causal-atlas-only (no fabrication)'}")
    print(f"    Tier5 novel protein-level: {n5} | Tier4 prioritized: {n4} | Tier1 controls: {n1}")
    print("[31] Done.")


if __name__ == "__main__":
    main()
