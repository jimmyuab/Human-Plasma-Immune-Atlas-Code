#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
44 - Expanded, RELIABLE UK Biobank cross-population replication
===============================================================
The original UK replication (src/33) used a hand-written map of 14 diseases.
This expands it to EVERY pan-phenome FinnGen R12 hit-disease that has an
EXACT-CODE Neale UKB counterpart on OpenGWAS (id == 'ukb-d-<FinnGen phenocode>').
Exact-code matching is used deliberately (not fuzzy name matching) so each
FinnGen<->UKB pair is the SAME harmonised endpoint definition -- reliable, not
approximate.

For each such disease we re-test only the genes that were FinnGen-causal
(FDR<0.05) in that disease -- replicating the discovery nominations in an
independent (English) population, using the identical Wald-ratio MR and allele
harmonisation as src/33. FinnGen(Finland) vs UKB(England) same-direction and
UKB nominal significance define two-population validation.

Inputs:
  06_genetic_causality/cis_MR_ALL_finngen_results.tsv     (pan-phenome hits)
  06_genetic_causality/immune_instruments_hg38.tsv        (instrument alleles)
  .secrets/opengwas_token.txt
Output:
  06_genetic_causality/uk_panphenome_concordance.tsv

Run:  python src/44_uk_panphenome_replication.py [--limit N]
"""
import os
import sys
import json
import time
import math
import urllib.request

import numpy as np
import pandas as pd

ROOT = r"I:\Plasma immune atalas"
GC = os.path.join(ROOT, "06_genetic_causality")
SEC = os.path.join(ROOT, ".secrets", "opengwas_token.txt")
API = "https://api.opengwas.io/api"
FDR_SIG = 0.05


def token():
    return open(SEC).read().strip()


def get(path, tok, retries=3):
    req = urllib.request.Request(API + path,
                                 headers={"Authorization": "Bearer " + tok})
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if i == retries - 1:
                print(f"  [warn] GET {path}: {type(e).__name__}")
                return None
            time.sleep(3 * (i + 1))


def post(path, payload, tok, retries=4):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(API + path, data=data,
                                 headers={"Authorization": "Bearer " + tok,
                                          "Content-Type": "application/json"})
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if i == retries - 1:
                return None
            time.sleep(3 * (i + 1))


def main():
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])
    tok = token()

    # pan-phenome hits + per-gene exposure (b_exp, alleles from instruments)
    mr = pd.read_csv(os.path.join(GC, "cis_MR_ALL_finngen_results.tsv"), sep="\t")
    sig = mr[mr.FDR < FDR_SIG].copy()
    inst = pd.read_csv(os.path.join(GC, "immune_instruments_hg38.tsv"), sep="\t")
    inst["SNP"] = inst["SNP"].astype(str)
    inst_lu = inst.drop_duplicates("gene_symbol").set_index("gene_symbol")

    # available UKB datasets on OpenGWAS
    info = get("/gwasinfo", tok)
    ukb_ids = {d["id"] for d in info.values() if str(d.get("id", "")).startswith("ukb-")}

    # exact-code matches: ukb-d-<phenocode>
    hit_codes = sig[["phenocode", "phenotype"]].drop_duplicates()
    matches = []
    for _, r in hit_codes.iterrows():
        uid = "ukb-d-" + r.phenocode
        if uid in ukb_ids:
            matches.append((r.phenocode, r.phenotype, uid))
    if limit:
        matches = matches[:limit]
    print(f"[44] exact ukb-d matches among hit diseases: {len(matches)}", flush=True)

    rows = []
    for k, (code, pheno, uid) in enumerate(matches, 1):
        genes = sig.loc[sig.phenocode == code, "gene_symbol"].unique().tolist()
        # collect this disease's hit-gene instruments
        recs = []
        for g in genes:
            if g not in inst_lu.index:
                continue
            gr = inst_lu.loc[g]
            frow = sig[(sig.phenocode == code) & (sig.gene_symbol == g)].iloc[0]
            recs.append(dict(gene=g, SNP=str(gr["SNP"]),
                             ea=str(gr["effect_allele"]).upper(),
                             oa=str(gr["other_allele"]).upper(),
                             b_exp=float(frow["b_exp"]),
                             fin_OR=float(frow["OR"]), fin_FDR=float(frow["FDR"])))
        if not recs:
            continue
        rsids = sorted({r["SNP"] for r in recs})
        assoc = []
        for i in range(0, len(rsids), 8):
            res = post("/associations", {"variant": rsids[i:i + 8], "id": [uid]}, tok)
            if isinstance(res, list):
                assoc.extend(res)
            time.sleep(0.25)
        uk_lu = {}
        for a in assoc:
            s = str(a.get("rsid", ""))
            if s and s not in uk_lu:
                uk_lu[s] = a
        n_hit = 0
        for rec in recs:
            a = uk_lu.get(rec["SNP"])
            if not a:
                continue
            try:
                b_out = float(a["beta"]); se_out = float(a["se"])
                ea_uk = str(a.get("ea", "")).upper(); nea_uk = str(a.get("nea", "")).upper()
                p_snp = float(a.get("p", np.nan))
            except (TypeError, ValueError, KeyError):
                continue
            if ea_uk == rec["ea"] and nea_uk == rec["oa"]:
                sign = 1.0
            elif ea_uk == rec["oa"] and nea_uk == rec["ea"]:
                sign = -1.0
            else:
                continue
            if rec["b_exp"] == 0:
                continue
            mr_beta = (sign * b_out) / rec["b_exp"]
            mr_se = abs(se_out / rec["b_exp"])
            z = mr_beta / mr_se if mr_se > 0 else np.nan
            uk_p = math.erfc(abs(z) / math.sqrt(2)) if z == z else np.nan
            fin_dir = "+" if rec["fin_OR"] > 1 else "-"
            uk_OR = math.exp(mr_beta)
            uk_dir = "+" if uk_OR > 1 else "-"
            rows.append(dict(gene_symbol=rec["gene"], phenocode=code, disease=pheno,
                             uk_id=uid, fin_OR=rec["fin_OR"], fin_FDR=rec["fin_FDR"],
                             uk_OR=uk_OR, uk_p=uk_p, fin_dir=fin_dir, uk_dir=uk_dir,
                             concordant=(fin_dir == uk_dir),
                             uk_replicates=(fin_dir == uk_dir and uk_p < 0.05),
                             two_population_validated=(fin_dir == uk_dir and uk_p < 0.05)))
            n_hit += 1
        print(f"[44] {k}/{len(matches)} {code:22s} {uid:26s} tested={n_hit}", flush=True)

    d = pd.DataFrame(rows)
    out = os.path.join(GC, "uk_panphenome_concordance.tsv")
    d.to_csv(out, sep="\t", index=False)
    if len(d):
        print(f"\n[44] === EXPANDED UK REPLICATION ===")
        print(f"pairs tested: {len(d)} | diseases: {d.phenocode.nunique()} | "
              f"genes: {d.gene_symbol.nunique()}")
        print(f"concordant: {int(d.concordant.sum())} | "
              f"two-population validated: {int(d.two_population_validated.sum())}")
    print("[44] wrote", out)


if __name__ == "__main__":
    main()
