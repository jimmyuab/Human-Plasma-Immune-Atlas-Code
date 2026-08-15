#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
52 - Expanded two-population (FinnGen -> UK Biobank) replication
================================================================
src/44 replicated pan-phenome hits only where a FinnGen phenocode had an
EXACT-code Neale counterpart (id == 'ukb-d-<phenocode>'), covering 48 diseases.
This adds a second, still-conservative matching route -- normalised trait-name
matching against UK Biobank datasets on OpenGWAS -- and records how each pair
was matched so every pairing is auditable.

INDEPENDENCE RULE (why only UK Biobank):
Many OpenGWAS disease GWAS are meta-analyses that silently INCLUDE FinnGen
(e.g. Sakaue 2021, ebi-a-GCST90018*, is BBJ + UKB + FinnGen), which would make
"replication" circular. Only ukb-d / ukb-b / ukb-a / ukb-e datasets are used:
single-cohort UK Biobank, unambiguously FinnGen-free, matching the Finland-vs-
England two-population design. Name matches additionally require >=200 cases
and a European/unspecified ancestry label.

Statistics are identical to src/44 (Wald ratio on the same discovery
instrument, strand-aware harmonisation, direction + nominal significance).

Output: 06_genetic_causality/uk_panphenome_concordance_ALL.tsv
Run:    python src/52_uk_replication_expanded.py [--limit N]
"""
import os
import sys
import re
import json
import math
import time
import urllib.request

import numpy as np
import pandas as pd

ROOT = r"I:\Plasma immune atalas"
GC = os.path.join(ROOT, "06_genetic_causality")
SEC = os.path.join(ROOT, ".secrets", "opengwas_token.txt")
API = "https://api.opengwas.io/api"
FDR_SIG = 0.05
MIN_CASES = 200
UKB_PREFIX = ("ukb-d-", "ukb-b-", "ukb-a-", "ukb-e-")


def token():
    return open(SEC).read().strip()


def get(path, tok, retries=3):
    req = urllib.request.Request(API + path,
                                 headers={"Authorization": "Bearer " + tok})
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
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
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except Exception:  # noqa: BLE001
            if i == retries - 1:
                return None
            time.sleep(3 * (i + 1))


def norm(s):
    s = str(s).lower()
    s = re.sub(r"\(.*?\)", "", s)
    s = re.sub(r"^(diagnoses - main icd10:|non-cancer illness code, self-reported:)", "", s)
    s = re.sub(r"\b(fg|mode|more control exclusions|definitions combined|"
               r"main-diagnosis|only as)\b", "", s)
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    return " ".join(s.split())


def build_map(info, sig):
    """FinnGen phenocode -> (ukb_id, match_method, trait, ncase)."""
    ids = set(info)
    by_name = {}
    for d in info.values():
        i = str(d.get("id", ""))
        if not i.startswith(UKB_PREFIX):
            continue
        by_name.setdefault(norm(d.get("trait", "")), []).append(d)
    out = {}
    for code, pheno in zip(sig.phenocode, sig.phenotype):
        uid = "ukb-d-" + code
        if uid in ids:
            out[code] = (uid, "exact-code", info[uid].get("trait"),
                         info[uid].get("ncase") or 0)
            continue
        cands = [d for d in by_name.get(norm(pheno), [])
                 if (d.get("ncase") or 0) >= MIN_CASES
                 and str(d.get("population")) in ("European", "None", "NA", "nan")]
        if cands:
            b = max(cands, key=lambda d: d.get("ncase") or 0)
            out[code] = (b["id"], "name-match", b.get("trait"), b.get("ncase") or 0)
    return out


def main():
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])
    tok = token()

    mr = pd.read_csv(os.path.join(GC, "cis_MR_ALL_finngen_results.tsv"), sep="\t")
    sig = mr[mr.FDR < FDR_SIG].copy()
    dis = sig[["phenocode", "phenotype"]].drop_duplicates()
    inst = pd.read_csv(os.path.join(GC, "immune_instruments_hg38.tsv"), sep="\t")
    inst["SNP"] = inst["SNP"].astype(str)
    inst_lu = inst.drop_duplicates("gene_symbol").set_index("gene_symbol")

    cache = os.path.join(GC, "opengwas_gwasinfo.json")
    info = json.load(open(cache)) if os.path.exists(cache) else get("/gwasinfo", tok)
    if not os.path.exists(cache):
        json.dump(info, open(cache, "w"))

    dmap = build_map(info, dis)
    codes = [c for c in dis.phenocode if c in dmap]
    if limit:
        codes = codes[:limit]
    ex = sum(1 for c in codes if dmap[c][1] == "exact-code")
    print(f"[52] hit diseases: {len(dis)} | UKB-replicable: {len(codes)} "
          f"(exact-code {ex}, name-match {len(codes) - ex})", flush=True)

    rows = []
    for k, code in enumerate(codes, 1):
        uid, method, utrait, ncase = dmap[code]
        pheno = dis.loc[dis.phenocode == code, "phenotype"].iloc[0]
        sub = sig[sig.phenocode == code]
        recs = []
        for g in sub.gene_symbol.unique():
            if g not in inst_lu.index:
                continue
            gr = inst_lu.loc[g]
            fr = sub[sub.gene_symbol == g].iloc[0]
            recs.append(dict(gene=g, SNP=str(gr["SNP"]),
                             ea=str(gr["effect_allele"]).upper(),
                             oa=str(gr["other_allele"]).upper(),
                             b_exp=float(fr["b_exp"]),
                             fin_OR=float(fr["OR"]), fin_FDR=float(fr["FDR"])))
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
        n = 0
        for rec in recs:
            a = uk_lu.get(rec["SNP"])
            if not a:
                continue
            try:
                b_out = float(a["beta"]); se_out = float(a["se"])
                ea_uk = str(a.get("ea", "")).upper()
                nea_uk = str(a.get("nea", "")).upper()
            except (TypeError, ValueError, KeyError):
                continue
            if ea_uk == rec["ea"] and nea_uk == rec["oa"]:
                sign = 1.0
            elif ea_uk == rec["oa"] and nea_uk == rec["ea"]:
                sign = -1.0
            else:
                continue
            if rec["b_exp"] == 0 or se_out <= 0:
                continue
            mr_beta = (sign * b_out) / rec["b_exp"]
            mr_se = abs(se_out / rec["b_exp"])
            z = mr_beta / mr_se
            uk_p = math.erfc(abs(z) / math.sqrt(2))
            fin_dir = "+" if rec["fin_OR"] > 1 else "-"
            uk_OR = math.exp(mr_beta)
            uk_dir = "+" if uk_OR > 1 else "-"
            conc = fin_dir == uk_dir
            rows.append(dict(gene_symbol=rec["gene"], phenocode=code, disease=pheno,
                             uk_id=uid, uk_trait=utrait, match_method=method,
                             uk_ncase=ncase, fin_OR=rec["fin_OR"], fin_FDR=rec["fin_FDR"],
                             uk_OR=uk_OR, uk_p=uk_p, fin_dir=fin_dir, uk_dir=uk_dir,
                             concordant=conc,
                             uk_replicates=(conc and uk_p < 0.05),
                             two_population_validated=(conc and uk_p < 0.05)))
            n += 1
        print(f"[52] {k}/{len(codes)} {code:24s} {uid:26s} {method:11s} tested={n}",
              flush=True)

    d = pd.DataFrame(rows)
    out = os.path.join(GC, "uk_panphenome_concordance_ALL.tsv")
    d.to_csv(out, sep="\t", index=False)
    if len(d):
        print("\n[52] === EXPANDED TWO-POPULATION REPLICATION ===")
        print(f"pairs tested: {len(d)} | diseases: {d.phenocode.nunique()} | "
              f"genes: {d.gene_symbol.nunique()}")
        print(f"concordant: {int(d.concordant.sum())} "
              f"({100 * d.concordant.mean():.0f}%) | "
              f"two-population validated: {int(d.two_population_validated.sum())}")
        print(d.match_method.value_counts().to_string())
    print("[52] wrote", out)


if __name__ == "__main__":
    main()
