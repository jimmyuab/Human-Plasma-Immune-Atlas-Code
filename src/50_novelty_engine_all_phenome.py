#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
50 - NOVELTY PRIORITY ENGINE across the WHOLE FinnGen phenome
=============================================================
src/28 scored only the 28 deeply-annotated core diseases (176 pairs). This
rescores EVERY pan-phenome causal hit (cis_MR_ALL_finngen_results.tsv,
FDR<0.05) with the identical, auditable formula, and adds the two evidence
layers that now also cover the whole phenome:

  NoveltyPriority = s_causal + s_coloc + s_pleio + s_drug + s_cell
                    + s_protein + s_replication
                    - known_drug_penalty - MHC_LD_penalty

s_causal/s_coloc/s_pleio/s_drug/s_cell/penalties are IDENTICAL to src/28 so the
core-28 ranking is reproduced exactly when the two new terms are zero.
  s_protein     1.0 direction-concordant protein pQTL-MR (p<0.05)
                0.5 protein cis-pQTL colocalizes (PP.H4>=0.8) but MR not sig
  s_replication 0.5 two-population (FinnGen + UK Biobank) validated

Inputs:
  06_genetic_causality/cis_MR_ALL_finngen_results.tsv
  06_genetic_causality/coloc_ALL_finngen_results.tsv (+ legacy coloc files)
  06_genetic_causality/pqtl_MR_ALL_finngen_results.tsv
  06_genetic_causality/pqtl_coloc_ALL_finngen_results.tsv
  06_genetic_causality/uk_panphenome_concordance_ALL.tsv (+ legacy exact-code file)
  02_data_processed/plasma_immune_protein_annotation.tsv
Outputs:
  06_genetic_causality/novelty_engine_ranked_ALL.tsv
  09_tables/T7_novelty_priority_targets_ALL.tsv

Run:  python src/50_novelty_engine_all_phenome.py
"""
import os
import re

import numpy as np
import pandas as pd

ROOT = r"I:\Plasma immune atalas"
GEN = os.path.join(ROOT, "06_genetic_causality")
PROC = os.path.join(ROOT, "02_data_processed")
TAB = os.path.join(ROOT, "09_tables")
os.makedirs(TAB, exist_ok=True)

KNOWN = {"IL6ST", "CTLA4", "TNFRSF1A", "IL4", "TNFSF14", "TNFRSF14",
         "TNFRSF1B", "IL6R", "TYK2", "JAK2"}

SYSTEM = [
    ("Cardiovascular", r"circulatory|Cardiometabolic|hypertension"),
    ("Metabolic/Endocrine", r"Endocrine|Diabetes|metabolic"),
    ("Respiratory", r"respiratory|Asthma|COPD"),
    ("Musculoskeletal", r"musculoskeletal|Rheuma"),
    ("Dermatological", r"skin and subcutaneous"),
    ("Neurological/Psychiatric", r"nervous system|Mental|Neurolog"),
    ("Ophthalmological", r"eye and adnexa"),
    ("Gastrointestinal", r"digestive system|Gastrointestinal"),
    ("Renal/Genitourinary", r"genitourinary|Kidney|renal"),
    ("Haematological/Immune", r"blood and blood-forming|immune mechanism|autimmune|autoimmune"),
    ("Neoplasm", r"Neoplasm|cancer"),
    ("Infection", r"infectious|parasitic"),
    ("Medication/Other", r"Drug purchase|Operation|Comorbidities"),
]


def system_of(cat):
    c = str(cat)
    for name, pat in SYSTEM:
        if re.search(pat, c, flags=re.I):
            return name
    return "Other"


def load_coloc():
    frames = []
    for fn in ("coloc_ALL_finngen_results.tsv", "coloc_phenome_results.tsv",
               "coloc_results.tsv"):
        fp = os.path.join(GEN, fn)
        if os.path.exists(fp):
            d = pd.read_csv(fp, sep="\t")
            if {"gene", "disease", "PP_H4"} <= set(d.columns):
                frames.append(d[["gene", "disease", "PP_H4"]])
    col = pd.concat(frames, ignore_index=True).dropna(subset=["PP_H4"])
    col["gene"] = col.gene.astype(str).str.upper()
    col["disease"] = col.disease.astype(str).str.lower()
    return (col.sort_values("PP_H4", ascending=False)
               .drop_duplicates(["gene", "disease"])
               .set_index(["gene", "disease"]).PP_H4.to_dict())


def load_pqtl():
    """(gene,disease) -> dict(pqtl_OR, pqtl_p, pqtl_PPH4)."""
    out = {}
    for fn in ("pqtl_MR_ALL_finngen_results.tsv", "pqtl_MR_results.tsv"):
        fp = os.path.join(GEN, fn)
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp, sep="\t")
        d = d.sort_values("MR_p")
        for _, r in d.iterrows():
            k = (str(r.gene).upper(), str(r.disease).lower())
            if k in out and "pqtl_p" in out[k]:
                continue
            out.setdefault(k, {}).update(pqtl_OR=float(r.OR), pqtl_p=float(r.MR_p))
    for fn in ("pqtl_coloc_ALL_finngen_results.tsv", "pqtl_coloc_results.tsv"):
        fp = os.path.join(GEN, fn)
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp, sep="\t").dropna(subset=["PP_H4"])
        for _, r in d.sort_values("PP_H4", ascending=False).iterrows():
            k = (str(r.gene).upper(), str(r.disease).lower())
            e = out.setdefault(k, {})
            e["pqtl_PPH4"] = max(e.get("pqtl_PPH4", 0.0), float(r.PP_H4))
    return out


def load_repl():
    """Union of the exact-code (src/44) and expanded UKB (src/52) replications."""
    out = set()
    for fn in ("uk_panphenome_concordance_ALL.tsv", "uk_panphenome_concordance.tsv"):
        fp = os.path.join(GEN, fn)
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp, sep="\t")
        ok = d[d.two_population_validated.astype(str).str.lower().isin(("true", "1", "yes"))]
        out |= {(str(r.gene_symbol).upper(), str(r.disease).lower()) for _, r in ok.iterrows()}
    return out


def main():
    mr = pd.read_csv(os.path.join(GEN, "cis_MR_ALL_finngen_results.tsv"), sep="\t")
    df = mr[mr.FDR < 0.05].copy()
    df["gene_symbol"] = df.gene_symbol.astype(str).str.upper()
    df = df.sort_values("FDR").drop_duplicates(["gene_symbol", "phenocode"])
    df = df.rename(columns={"phenotype": "disease"})
    df["disease_system"] = df.category.map(system_of)
    print(f"[50] pan-phenome causal hits: {len(df)} | "
          f"{df.disease.nunique()} diseases | {df.gene_symbol.nunique()} proteins", flush=True)

    ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
    acls = ann.set_index("gene_symbol").hpa_protein_class.fillna("").astype(str)
    src = ann.set_index("gene_symbol").immune_source_cells.fillna("")

    def flag(g, key):
        return int(key.lower() in acls.get(g, "").lower())

    df["fda_target"] = df.gene_symbol.map(lambda g: flag(g, "FDA approved drug target"))
    df["is_secreted"] = df.gene_symbol.map(lambda g: flag(g, "secreted"))
    df["is_membrane"] = df.gene_symbol.map(lambda g: flag(g, "membrane"))
    df["druggability"] = (2 * df.fda_target + df.is_secreted + df.is_membrane).clip(upper=3)
    df["has_cellsource"] = df.gene_symbol.map(
        lambda g: int(str(src.get(g, "")).strip() != ""))

    gcat = df.groupby("gene_symbol").disease_system.nunique()
    df["n_categories"] = df.gene_symbol.map(gcat).fillna(1).astype(int)

    coloc = load_coloc()
    pq = load_pqtl()
    repl = load_repl()

    keys = list(zip(df.gene_symbol, df.disease.astype(str).str.lower()))
    df["PP_H4"] = [coloc.get(k, np.nan) for k in keys]
    df["pqtl_OR"] = [pq.get(k, {}).get("pqtl_OR", np.nan) for k in keys]
    df["pqtl_p"] = [pq.get(k, {}).get("pqtl_p", np.nan) for k in keys]
    df["pqtl_PPH4"] = [pq.get(k, {}).get("pqtl_PPH4", np.nan) for k in keys]
    df["two_population_validated"] = [k in repl for k in keys]

    same_dir = ((df.OR > 1) & (df.pqtl_OR > 1)) | ((df.OR < 1) & (df.pqtl_OR < 1))
    df["pqtl_concordant"] = same_dir & (df.pqtl_p < 0.05)
    df["pqtl_tested"] = df.pqtl_p.notna() | df.pqtl_PPH4.notna()

    df["known_axis"] = df.gene_symbol.isin(KNOWN).astype(int)
    df["mhc_penalty"] = (df.chr == 6).astype(int)

    df["s_causal"] = (-np.log10(df.FDR.clip(lower=1e-300))).clip(upper=10) / 5.0
    df["s_coloc"] = df.PP_H4.fillna(0)
    df["s_pleio"] = (df.n_categories - 1).clip(upper=3) * 0.5
    df["s_drug"] = df.druggability * 0.5
    df["s_cell"] = df.has_cellsource * 0.5
    df["s_protein"] = np.where(df.pqtl_concordant, 1.0,
                               np.where(df.pqtl_PPH4.fillna(0) >= 0.8, 0.5, 0.0))
    df["s_repl"] = df.two_population_validated.astype(int) * 0.5
    df["p_known"] = df.known_axis * 1.5
    df["p_mhc"] = df.mhc_penalty * 2.0
    df["novelty_priority"] = (df.s_causal + df.s_coloc + df.s_pleio + df.s_drug
                              + df.s_cell + df.s_protein + df.s_repl
                              - df.p_known - df.p_mhc)

    df["category_label"] = np.where(
        df.known_axis == 1, "recovered known axis",
        np.where(df.mhc_penalty == 1, "MHC-caution",
                 np.where((df.PP_H4.fillna(0) >= 0.8) & df.pqtl_concordant,
                          "NOVEL protein-confirmed",
                          np.where(df.PP_H4.fillna(0) >= 0.8, "NOVEL colocalized",
                                   "novel nomination"))))

    rank = df.sort_values("novelty_priority", ascending=False)
    keep = ["gene_symbol", "disease", "phenocode", "category", "disease_system",
            "immune_class", "OR", "FDR", "PP_H4", "pqtl_OR", "pqtl_p", "pqtl_PPH4",
            "pqtl_concordant", "two_population_validated", "n_categories",
            "druggability", "fda_target", "category_label", "novelty_priority",
            "s_causal", "s_coloc", "s_pleio", "s_drug", "s_cell", "s_protein",
            "s_repl", "p_known", "p_mhc"]
    out = os.path.join(GEN, "novelty_engine_ranked_ALL.tsv")
    rank[keep].to_csv(out, sep="\t", index=False)

    pub = rank[rank.category_label.str.startswith("NOVEL")]
    pub[["gene_symbol", "disease", "disease_system", "immune_class", "OR", "FDR",
         "PP_H4", "pqtl_p", "pqtl_PPH4", "n_categories", "druggability",
         "category_label", "novelty_priority"]].head(100).to_csv(
        os.path.join(TAB, "T7_novelty_priority_targets_ALL.tsv"), sep="\t", index=False)

    print("\n[50] === PAN-PHENOME NOVELTY ENGINE ===")
    vc = df.category_label.value_counts()
    print(vc.to_string())
    print(f"protein-layer tested pairs: {int(df.pqtl_tested.sum())} | "
          f"concordant: {int(df.pqtl_concordant.sum())} | "
          f"two-population validated: {int(df.two_population_validated.sum())}")
    show = ["gene_symbol", "disease", "disease_system", "OR", "FDR", "PP_H4",
            "pqtl_p", "category_label", "novelty_priority"]
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print("\nTop 30:")
        print(rank[show].head(30).to_string(index=False))
    print("\n[50] wrote", out)


if __name__ == "__main__":
    main()
