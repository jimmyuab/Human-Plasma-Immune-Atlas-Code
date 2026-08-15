#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
36 - Plasma immune protein -> circulating immune-CELL-COUNT Mendelian randomization
===================================================================================

The atlas is a plasma *protein* resource, but immune cell counts are a real,
public immune phenotype. This module adds a cellular layer that mirrors the
"trait" arm of large phenome atlases (metabolite->trait): it asks which plasma
immune proteins CAUSALLY drive the circulating immune cell populations, by
running the SAME immune cis-eQTL instruments used everywhere else against public
UK Biobank blood immune-cell-count GWAS (Neale lab, served by OpenGWAS).

Scope is ALL 673 cis-eQTL-instrumented immune proteins (every protein with a
usable instrument), tested against 6 immune cell counts. Only public summary
statistics and the user's own OpenGWAS token are used; nothing is fabricated,
and no existing file is modified (new outputs only).

Run:
    python src/36_immune_cell_count_mr.py
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
#  Immune cell count -> public UK Biobank (Neale) GWAS id
# --------------------------------------------------------------------------- #
CELL_GWAS = {
    "Lymphocyte count": ("ukb-d-30120_irnt", "Neale UKB, lymphocyte count"),
    "Monocyte count":   ("ukb-d-30130_irnt", "Neale UKB, monocyte count"),
    "Neutrophil count": ("ukb-d-30140_irnt", "Neale UKB, neutrophil count"),
    "Eosinophil count": ("ukb-d-30150_irnt", "Neale UKB, eosinophil count"),
    "Basophil count":   ("ukb-d-30160_irnt", "Neale UKB, basophil count"),
    "Leukocyte count":  ("ukb-d-30000_irnt", "Neale UKB, white blood cell count"),
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
                    "cell_beta": mr_beta, "cell_se": mr_se, "cell_z": z, "cell_p": p,
                    "cell_snp_p": p_out})
    if not out:
        return None
    d = pd.DataFrame(out)
    d["cell_trait"] = trait
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
    print("[36] Immune protein -> immune-cell-count MR ...")
    token = load_token()
    instruments = build_instruments()
    print(f"    immune instruments available: {len(instruments)} genes")
    print(f"    immune cell-count GWAS mapped: {len(CELL_GWAS)}")

    allc = []
    for trait, (uk_id, label) in CELL_GWAS.items():
        print(f"    querying {trait:18s} -> {uk_id} ...", end="", flush=True)
        d = uk_mr_for_trait(trait, uk_id, instruments, token)
        if d is None or not len(d):
            print(" no usable data (skipped)")
            continue
        d["cell_label"] = label
        allc.append(d)
        print(f" {len(d)} instruments tested")

    if not allc:
        print("[36] No cell-count data retrieved. Nothing written -- no fabrication.")
        return

    cc = pd.concat(allc, ignore_index=True)
    cc["cell_FDR"] = bh_fdr(cc["cell_p"].values)
    cols = ["gene_symbol", "cell_trait", "uk_id", "cell_label", "SNP",
            "cell_beta", "cell_se", "cell_z", "cell_p", "cell_FDR", "cell_snp_p"]
    cc = cc[cols].sort_values("cell_p")
    cc.to_csv(os.path.join(GC, "immune_cell_count_MR_results.tsv"), sep="\t", index=False)
    n_hit = int((cc["cell_FDR"] < FDR_SIG).sum())
    print(f"    immune-cell-count MR -> immune_cell_count_MR_results.tsv "
          f"({len(cc)} tests, {n_hit} at FDR<0.05)")

    # short per-trait summary
    summ = (cc[cc["cell_FDR"] < FDR_SIG]
            .groupby("cell_trait").size().sort_values(ascending=False))
    for t, n in summ.items():
        print(f"      {t:18s}: {n} causal immune proteins (FDR<0.05)")
    print("[36] Done.")


if __name__ == "__main__":
    main()
