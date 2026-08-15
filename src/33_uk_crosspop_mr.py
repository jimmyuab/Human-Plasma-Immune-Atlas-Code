#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
33 - UK (England) cross-population Mendelian randomization
==========================================================

"Why only Finland?" -- this module answers it. It re-runs the SAME immune
cis-eQTL instruments used for the FinnGen discovery against *public UK Biobank
disease GWAS* (Neale lab / large curated UKB studies, served by OpenGWAS), so
every FinnGen causal nomination can be re-tested in an independent English
population. It then reports Finland <-> England directional concordance and the
subset that is significant in BOTH populations ("two-population validated").

Individual-level UKB data is gated; UKB *summary statistics* are public. This
uses only the summary stats, via the user's own OpenGWAS token, for its intended
purpose. No individual-level data, no fabrication -- diseases with no usable UK
GWAS are skipped and reported as such.

Run:
    python src/33_uk_crosspop_mr.py
"""

import os
import time
import json
import urllib.request
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GC   = os.path.join(ROOT, "06_genetic_causality")
SEC  = os.path.join(ROOT, ".secrets", "opengwas_token.txt")

API  = "https://api.opengwas.io/api"
FDR_SIG = 0.05

# --------------------------------------------------------------------------- #
#  Disease -> best public UK Biobank / UK-population GWAS (OpenGWAS id)
#  Chosen for the clearest disease definition + usable case count. Proxy or
#  weak definitions are flagged so claims stay honest.
# --------------------------------------------------------------------------- #
UK_GWAS = {
    "Rheumatoid arthritis":       ("ukb-d-M13_RHEUMA",  "Neale UKB, RA (ICD M13)", ""),
    "Coronary heart disease":     ("ukb-d-I9_CHD",       "Neale UKB, major CHD event", ""),
    "Heart failure":              ("ukb-d-I9_HEARTFAIL", "Neale UKB, heart failure", ""),
    "Venous thromboembolism":     ("ukb-d-I9_VTE",       "Neale UKB, VTE", ""),
    "Glaucoma":                   ("ukb-d-H40",          "Neale UKB, glaucoma (ICD H40)", ""),
    "Atrial fibrillation":        ("ukb-b-964",          "UKB, AF/flutter (ICD I48)", ""),
    "Hypertension":               ("ieu-b-5144",         "Tang UKB, hypertension", ""),
    "Psoriasis":                  ("ukb-b-10537",        "UKB, psoriasis (self-report)", "self-report"),
    "Coeliac disease":            ("ukb-b-8631",         "UKB, coeliac (self-report)", "self-report"),
    "Ankylosing spondylitis":     ("ukb-b-18194",        "UKB, AS (self-report)", "self-report"),
    "Type 2 diabetes":            ("ukb-b-13806",        "UKB, T2D (self-report)", "self-report"),
    "Osteoporosis":               ("ukb-b-12141",        "UKB, osteoporosis (self-report)", "self-report"),
    "Multiple sclerosis":         ("ukb-b-17670",        "UKB, MS (self-report)", "self-report"),
    "Crohn's disease":            ("ukb-a-552",          "Neale UKB, Crohn (ICD K50)", ""),
}


def load_token():
    if not os.path.exists(SEC):
        raise SystemExit("No OpenGWAS token at .secrets/opengwas_token.txt")
    return open(SEC).read().strip()


def post(path, payload, token, retries=3):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(API + path, data=data,
                                 headers={"Authorization": "Bearer " + token,
                                          "Content-Type": "application/json"})
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if i == retries - 1:
                print(f"    [warn] {path} failed: {type(e).__name__} {e}")
                return None
            time.sleep(3 * (i + 1))
    return None


# --------------------------------------------------------------------------- #
#  Build the immune instrument table (gene, SNP, alleles, eQTL beta)
# --------------------------------------------------------------------------- #
def build_instruments():
    inst = pd.read_csv(os.path.join(GC, "immune_cis_eqtl_instruments.tsv"), sep="\t")
    mr   = pd.read_csv(os.path.join(GC, "cis_MR_phenome_results.tsv"), sep="\t")
    # Restrict the cross-population test to genes that were FinnGen causal hits
    # (FDR < 0.05). Re-testing the discovery nominations in England is the point;
    # it also keeps the (slow) OpenGWAS lookups tractable.
    hit_genes = set(mr.loc[mr["FDR"] < FDR_SIG, "gene_symbol"])
    mr = mr[mr["gene_symbol"].isin(hit_genes)]
    inst = inst[inst["gene_symbol"].isin(hit_genes)]
    # one b_exp per gene (identical across diseases) from the phenome MR run
    bexp = (mr[["gene_symbol", "SNP", "eaf", "b_exp"]]
            .drop_duplicates("gene_symbol").set_index("gene_symbol"))
    rows = []
    for _, r in inst.iterrows():
        g = r["gene_symbol"]
        if g not in bexp.index:
            continue
        b = bexp.loc[g]
        if str(b["SNP"]) != str(r["SNP"]):
            continue  # use the exact instrument the discovery used
        rows.append({"gene_symbol": g, "SNP": r["SNP"],
                     "ea": str(r["effect_allele"]).upper(),
                     "oa": str(r["other_allele"]).upper(),
                     "eaf": float(b["eaf"]), "b_exp": float(b["b_exp"])})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
#  Query UK GWAS for the instrument SNPs and run Wald-ratio MR
# --------------------------------------------------------------------------- #
def uk_mr_for_disease(disease, uk_id, instruments, token):
    # unique SNPs only -- the OpenGWAS /associations endpoint is slow and times
    # out on large variant lists, so query in small chunks of unique rsIDs.
    rsids = sorted(set(str(s) for s in instruments["SNP"].tolist()))
    assoc = []
    for i in range(0, len(rsids), 8):
        chunk = rsids[i:i + 8]
        res = post("/associations", {"variant": chunk, "id": [uk_id]}, token)
        if isinstance(res, list):
            assoc.extend(res)
        time.sleep(0.3)
    if not assoc:
        return None
    a = pd.DataFrame(assoc)
    if "rsid" not in a.columns or "beta" not in a.columns:
        return None
    # UK association lookup keyed by rsid (first occurrence per SNP)
    uk_lu = {}
    for _, r in a.iterrows():
        snp = str(r["rsid"])
        if snp in uk_lu:
            continue
        uk_lu[snp] = r
    from math import erf, sqrt
    out = []
    for _, ins in instruments.iterrows():   # one row per gene -> keeps duplicates distinct
        snp = str(ins["SNP"])
        if snp not in uk_lu:
            continue
        r = uk_lu[snp]
        try:
            b_out = float(r["beta"]); se_out = float(r["se"])
            ea_uk = str(r.get("ea", "")).upper(); nea_uk = str(r.get("nea", "")).upper()
            p_out = float(r.get("p", np.nan))
        except (TypeError, ValueError):
            continue
        # harmonise UK effect to the eQTL effect allele
        if ea_uk == ins["ea"] and nea_uk == ins["oa"]:
            sign = 1.0
        elif ea_uk == ins["oa"] and nea_uk == ins["ea"]:
            sign = -1.0
        else:
            continue  # ambiguous / strand issue -> drop (honest)
        b_exp = float(ins["b_exp"])
        if b_exp == 0:
            continue
        mr_beta = (sign * b_out) / b_exp
        mr_se   = abs(se_out / b_exp)
        z = mr_beta / mr_se if mr_se > 0 else np.nan
        p = 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2)))) if z == z else np.nan
        out.append({"gene_symbol": ins["gene_symbol"], "SNP": snp,
                    "uk_beta": mr_beta, "uk_se": mr_se, "uk_z": z, "uk_p": p,
                    "uk_OR": float(np.exp(mr_beta)), "uk_snp_p": p_out})
    if not out:
        return None
    d = pd.DataFrame(out)
    d["disease"] = disease
    d["uk_id"] = uk_id
    return d


# --------------------------------------------------------------------------- #
def bh_fdr(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    idx = np.where(ok)[0]
    ps = p[idx]
    order = np.argsort(ps)
    m = len(ps)
    ranked = ps[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[idx[order]] = np.clip(ranked, 0, 1)
    return q


def main():
    print("[33] UK (England) cross-population MR ...")
    token = load_token()
    instruments = build_instruments()
    print(f"    immune instruments available: {len(instruments)} genes")
    print(f"    UK Biobank disease GWAS mapped: {len(UK_GWAS)}")

    all_uk = []
    for disease, (uk_id, label, caveat) in UK_GWAS.items():
        print(f"    querying {disease:28s} -> {uk_id} ...", end="", flush=True)
        d = uk_mr_for_disease(disease, uk_id, instruments, token)
        if d is None or not len(d):
            print(" no usable data (skipped)")
            continue
        d["uk_label"] = label
        d["uk_caveat"] = caveat
        all_uk.append(d)
        print(f" {len(d)} instruments tested")

    if not all_uk:
        print("[33] No UK data retrieved (token/rate limit?). Nothing written -- no fabrication.")
        return

    uk = pd.concat(all_uk, ignore_index=True)
    uk["uk_FDR"] = bh_fdr(uk["uk_p"].values)
    uk_cols = ["gene_symbol", "disease", "uk_id", "uk_label", "uk_caveat",
               "SNP", "uk_OR", "uk_beta", "uk_se", "uk_p", "uk_FDR", "uk_snp_p"]
    uk = uk[uk_cols].sort_values("uk_p")
    uk.to_csv(os.path.join(GC, "uk_cis_MR_results.tsv"), sep="\t", index=False)
    print(f"    UK MR results -> uk_cis_MR_results.tsv ({len(uk)} tests)")

    # -------- Finland <-> England concordance on the FinnGen hits ---------- #
    fin = pd.read_csv(os.path.join(GC, "cis_MR_phenome_results.tsv"), sep="\t")
    fin_hits = fin[fin["FDR"] < FDR_SIG][["gene_symbol", "disease", "OR", "FDR"]].copy()
    fin_hits = fin_hits.rename(columns={"OR": "fin_OR", "FDR": "fin_FDR"})
    merged = fin_hits.merge(uk, on=["gene_symbol", "disease"], how="inner")
    if len(merged):
        merged["fin_dir"] = np.where(merged["fin_OR"] > 1, "risk", "protective")
        merged["uk_dir"]  = np.where(merged["uk_OR"] > 1, "risk", "protective")
        merged["concordant"] = merged["fin_dir"] == merged["uk_dir"]
        merged["uk_replicates"] = merged["concordant"] & (merged["uk_p"] < FDR_SIG)
        merged["two_population_validated"] = (merged["concordant"]
                                              & (merged["uk_FDR"] < FDR_SIG))
        keep = ["gene_symbol", "disease", "fin_OR", "fin_FDR", "uk_OR", "uk_p",
                "uk_FDR", "fin_dir", "uk_dir", "concordant", "uk_replicates",
                "two_population_validated", "uk_id", "uk_label", "uk_caveat"]
        merged[keep].sort_values(["two_population_validated", "concordant", "uk_p"],
                                 ascending=[False, False, True]).to_csv(
            os.path.join(GC, "finngen_vs_uk_concordance.tsv"), sep="\t", index=False)
        n_cov = len(merged)
        n_conc = int(merged["concordant"].sum())
        n_rep = int(merged["uk_replicates"].sum())
        n_val = int(merged["two_population_validated"].sum())
        print(f"    Finland-hit coverage in UK: {n_cov} pairs")
        print(f"    directionally concordant : {n_conc}/{n_cov}")
        print(f"    UK-replicated (P<0.05)   : {n_rep}")
        print(f"    two-population validated (UK FDR<0.05): {n_val}")
        print(f"    -> finngen_vs_uk_concordance.tsv")
    else:
        print("    No FinnGen hits overlapped the mapped UK diseases.")
    print("[33] Done.")


if __name__ == "__main__":
    main()
