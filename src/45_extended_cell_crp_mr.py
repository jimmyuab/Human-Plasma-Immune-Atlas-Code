#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
45 - EXTENDED plasma immune protein -> blood immune-trait Mendelian randomization
=================================================================================
"No skip" expansion of the cellular arm (src/36). src/36 tested 6 immune cell
counts but basophil/eosinophil *count* GWAS are not in the OpenGWAS Neale set,
so they were dropped. This module (additive; writes a NEW output file, touches
nothing existing) covers EVERY public immune blood trait we can reliably reach:

  Counts     (Neale UKB _irnt) : leukocyte, lymphocyte, monocyte, neutrophil
  Counts     (Astle 2016, fills the Neale gap) : eosinophil, basophil
  Extra cell (Neale UKB _irnt) : platelet, reticulocyte, red blood cell
  Fractions  (Neale UKB _irnt) : lymphocyte%, monocyte%, neutrophil%,
                                 eosinophil%, basophil%
  Inflammation(Neale UKB _irnt): C-reactive protein (CRP)

Every trait uses public summary statistics + the user's own OpenGWAS token.
Same immune cis-eQTL instruments and identical Wald-ratio MR / allele
harmonisation as src/36. Nothing is fabricated; gated data is never used.

Inputs:
  06_genetic_causality/immune_cis_eqtl_instruments.tsv
  06_genetic_causality/cis_MR_phenome_results.tsv     (per-gene exposure beta)
  .secrets/opengwas_token.txt
Output:
  06_genetic_causality/extended_cell_crp_MR_results.tsv

Run:  python src/45_extended_cell_crp_mr.py
"""

import os
import time
import json
import urllib.request
from math import erf, sqrt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GC   = os.path.join(ROOT, "06_genetic_causality")
SEC  = os.path.join(ROOT, ".secrets", "opengwas_token.txt")

API  = "https://api.opengwas.io/api"
FDR_SIG = 0.05

# --------------------------------------------------------------------------- #
#  Immune blood trait -> public GWAS id (all real, all queryable via OpenGWAS)
#  group tag lets us summarise counts vs fractions vs inflammation separately.
# --------------------------------------------------------------------------- #
TRAIT_GWAS = {
    # ---- absolute cell counts (Neale UKB) --------------------------------- #
    "Leukocyte count":     ("ukb-d-30000_irnt",     "count", "Neale UKB white blood cell count"),
    "Lymphocyte count":    ("ukb-d-30120_irnt",     "count", "Neale UKB lymphocyte count"),
    "Monocyte count":      ("ukb-d-30130_irnt",     "count", "Neale UKB monocyte count"),
    "Neutrophil count":    ("ukb-d-30140_irnt",     "count", "Neale UKB neutrophil count"),
    # ---- counts the Neale set lacks -> Astle 2016 (European) -------------- #
    "Eosinophil count":    ("ebi-a-GCST004606",     "count", "Astle 2016 eosinophil count (n=172,275)"),
    "Basophil count":      ("ebi-a-GCST004618",     "count", "Astle 2016 basophil count (n=171,846)"),
    # ---- extra blood cell counts (Neale UKB) ------------------------------ #
    "Platelet count":      ("ukb-d-30080_irnt",     "count", "Neale UKB platelet count"),
    "Reticulocyte count":  ("ukb-d-30250_irnt",     "count", "Neale UKB reticulocyte count"),
    "Red blood cell count":("ukb-d-30010_irnt",     "count", "Neale UKB red blood cell count"),
    # ---- differential fractions / percentages (Neale UKB) ----------------- #
    "Lymphocyte %":        ("ukb-d-30180_irnt",     "fraction", "Neale UKB lymphocyte percentage"),
    "Monocyte %":          ("ukb-d-30190_irnt",     "fraction", "Neale UKB monocyte percentage"),
    "Neutrophil %":        ("ukb-d-30200_irnt",     "fraction", "Neale UKB neutrophil percentage"),
    "Eosinophil %":        ("ukb-d-30210_irnt",     "fraction", "Neale UKB eosinophil percentage"),
    "Basophil %":          ("ukb-d-30220_irnt",     "fraction", "Neale UKB basophil percentage"),
    # ---- systemic inflammation (Neale UKB) -------------------------------- #
    "C-reactive protein":  ("ukb-d-30710_irnt",     "inflammation", "Neale UKB C-reactive protein"),
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
        except Exception as e:  # noqa: BLE001
            if i == retries - 1:
                print(f"    [warn] {path} failed: {type(e).__name__} {e}")
                return None
            time.sleep(3 * (i + 1))
    return None


# --------------------------------------------------------------------------- #
#  ALL instrumented immune proteins (one row per gene, with eQTL exposure beta)
# --------------------------------------------------------------------------- #
def build_instruments():
    inst = pd.read_csv(os.path.join(GC, "immune_cis_eqtl_instruments.tsv"), sep="\t")
    mr   = pd.read_csv(os.path.join(GC, "cis_MR_phenome_results.tsv"), sep="\t")
    bexp = (mr[["gene_symbol", "SNP", "eaf", "b_exp"]]
            .drop_duplicates("gene_symbol").set_index("gene_symbol"))
    rows = []
    for _, r in inst.iterrows():
        g = r["gene_symbol"]
        if g not in bexp.index:
            continue
        b = bexp.loc[g]
        if str(b["SNP"]) != str(r["SNP"]):
            continue
        rows.append({"gene_symbol": g, "SNP": r["SNP"],
                     "ea": str(r["effect_allele"]).upper(),
                     "oa": str(r["other_allele"]).upper(),
                     "eaf": float(b["eaf"]), "b_exp": float(b["b_exp"])})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def uk_mr_for_trait(trait, uk_id, instruments, token):
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
    uk_lu = {}
    for _, r in a.iterrows():
        snp = str(r["rsid"])
        if snp not in uk_lu:
            uk_lu[snp] = r
    out = []
    for _, ins in instruments.iterrows():
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
        if ea_uk == ins["ea"] and nea_uk == ins["oa"]:
            sign = 1.0
        elif ea_uk == ins["oa"] and nea_uk == ins["ea"]:
            sign = -1.0
        else:
            continue
        b_exp = float(ins["b_exp"])
        if b_exp == 0:
            continue
        mr_beta = (sign * b_out) / b_exp
        mr_se   = abs(se_out / b_exp)
        z = mr_beta / mr_se if mr_se > 0 else np.nan
        p = 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2)))) if z == z else np.nan
        out.append({"gene_symbol": ins["gene_symbol"], "SNP": snp,
                    "trait_beta": mr_beta, "trait_se": mr_se, "trait_z": z,
                    "trait_p": p, "trait_snp_p": p_out})
    if not out:
        return None
    d = pd.DataFrame(out)
    d["trait"] = trait
    d["uk_id"] = uk_id
    return d


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
    print("[45] Extended immune protein -> blood immune-trait MR (no skip) ...")
    token = load_token()
    instruments = build_instruments()
    print(f"    immune instruments available: {len(instruments)} genes")
    print(f"    immune blood traits mapped:   {len(TRAIT_GWAS)}")

    groups = {t: g for t, (uid, g, lab) in TRAIT_GWAS.items()}
    labels = {t: lab for t, (uid, g, lab) in TRAIT_GWAS.items()}

    allc = []
    for trait, (uk_id, group, label) in TRAIT_GWAS.items():
        print(f"    querying {trait:22s} -> {uk_id:20s} ...", end="", flush=True)
        d = uk_mr_for_trait(trait, uk_id, instruments, token)
        if d is None or not len(d):
            print(" no usable data (skipped)")
            continue
        d["trait_group"] = group
        d["trait_label"] = label
        allc.append(d)
        print(f" {len(d)} instruments tested")

    if not allc:
        print("[45] No trait data retrieved. Nothing written -- no fabrication.")
        return

    cc = pd.concat(allc, ignore_index=True)
    cc["trait_FDR"] = bh_fdr(cc["trait_p"].values)
    cols = ["gene_symbol", "trait", "trait_group", "uk_id", "trait_label", "SNP",
            "trait_beta", "trait_se", "trait_z", "trait_p", "trait_FDR", "trait_snp_p"]
    cc = cc[cols].sort_values("trait_p")
    out = os.path.join(GC, "extended_cell_crp_MR_results.tsv")
    cc.to_csv(out, sep="\t", index=False)
    n_hit = int((cc["trait_FDR"] < FDR_SIG).sum())
    print(f"\n[45] extended blood-trait MR -> extended_cell_crp_MR_results.tsv "
          f"({len(cc)} tests, {cc['trait'].nunique()} traits, {n_hit} at FDR<0.05)")

    # per-trait summary
    summ = (cc[cc["trait_FDR"] < FDR_SIG]
            .groupby("trait").size().sort_values(ascending=False))
    for t, n in summ.items():
        print(f"      {t:22s}: {n} causal immune proteins (FDR<0.05)")
    print("[45] Done.")


if __name__ == "__main__":
    main()
