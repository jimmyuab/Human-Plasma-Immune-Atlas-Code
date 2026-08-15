#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
41 - PAN-PHENOME immune cis-MR across ALL FinnGen R12 endpoints
================================================================
Extends the 28-disease cis-eQTL(eQTLGen) -> disease(FinnGen R12) Wald-ratio MR
to EVERY endpoint in the FinnGen R12 manifest (~2,469 diseases), using remote
tabix point-queries (src/remote_tabix.py) so nothing is downloaded in full.

Same, already-validated statistics as src/24_run_mr_phenome.py:
  * Zhu (2016) Z -> beta/se reconstruction for the eQTL exposure
  * strand-aware allele harmonisation, palindromic SNPs dropped
  * Wald ratio MR, one cis instrument per gene

Purely additive. Instruments matched to FinnGen by rsID after hg19->hg38
liftover of the instrument position (only used to locate the tabix region).

Outputs (all inside project):
  06_genetic_causality/panphenome/<phenocode>.tsv   (per-endpoint, resumable)
  06_genetic_causality/cis_MR_ALL_finngen_results.tsv  (aggregated + FDR)

Run:
    python src/41_panphenome_mr.py            # full run (parallel, ~hours)
    python src/41_panphenome_mr.py --limit 5  # quick smoke test
"""
import os
import csv
import sys
import math
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remote_tabix import RemoteTabix  # noqa: E402

import pandas as pd  # noqa: E402

ROOT = r"I:\Plasma immune atalas"
GEN = os.path.join(ROOT, "06_genetic_causality")
MAN = os.path.join(ROOT, "01_data_raw", "FinnGen_GWAS", "finngen_R12_manifest.tsv")
PART = os.path.join(GEN, "panphenome")
os.makedirs(PART, exist_ok=True)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
N_WORKERS = 16      # endpoints processed concurrently
INNER_WORKERS = 12  # SNP queries per endpoint processed concurrently
# -> up to ~192 concurrent range requests

# FinnGen R12 sumstat columns (0-based)
C_CHROM, C_POS, C_REF, C_ALT, C_RSID = 0, 1, 2, 3, 4
C_PVAL, C_BETA, C_SEBETA, C_AF = 6, 8, 9, 12  # pval, beta, sebeta, af_alt_controls


def ambiguous(a, b):
    return COMP.get(a) == b


def zhu_beta_se(z, n, p):
    denom = 2 * p * (1 - p) * (n + z * z)
    if denom <= 0:
        return None, None
    s = math.sqrt(denom)
    return z / s, 1.0 / s


def load_instruments():
    ins = pd.read_csv(os.path.join(GEN, "immune_instruments_hg38.tsv"), sep="\t")
    ins["SNP"] = ins["SNP"].astype(str)
    recs = []
    for _, r in ins.iterrows():
        recs.append((str(r["chr38"]), int(r["pos38"]), str(r["SNP"]),
                     str(r["effect_allele"]).upper(), str(r["other_allele"]).upper(),
                     float(r["Zscore"]), float(r["NrSamples"]), str(r["gene_symbol"])))
    return recs


def _mr_one_snp(rt, inst):
    chrom, pos, rs, ea, oa, z, n, gene = inst
    if ambiguous(ea, oa):
        return None
    try:
        rows = rt.query(chrom, pos, fetch=49152)
    except Exception:
        return None
    hit = None
    for f in rows:  # match by rsID (robust to coordinate quirks)
        if len(f) > C_AF and rs in f[C_RSID]:
            hit = f
            break
    if hit is None:
        return None
    ref = hit[C_REF].upper()
    alt = hit[C_ALT].upper()
    try:
        beta = float(hit[C_BETA])
        se_out = float(hit[C_SEBETA])
        af = float(hit[C_AF])
    except (ValueError, IndexError):
        return None
    if alt == ea and ref == oa:
        bo, pe = beta, af
    elif ref == ea and alt == oa:
        bo, pe = -beta, 1 - af
    else:
        return None
    if not (0 < pe < 1):
        return None
    b_exp, _ = zhu_beta_se(z, n, pe)
    if b_exp is None or b_exp == 0:
        return None
    theta = bo / b_exp
    se_theta = se_out / abs(b_exp)
    zt = theta / se_theta
    pval = math.erfc(abs(zt) / math.sqrt(2))
    return {
        "gene_symbol": gene, "SNP": rs, "chr": chrom,
        "eqtl_Z": z, "eqtl_N": int(n), "eaf": round(pe, 4),
        "b_exp": b_exp, "b_out_oriented": bo, "se_out": se_out,
        "MR_beta": theta, "MR_se": se_theta, "MR_z": zt, "MR_p": pval,
        "OR": math.exp(theta),
        "OR_l95": math.exp(theta - 1.96 * se_theta),
        "OR_u95": math.exp(theta + 1.96 * se_theta),
    }


def mr_for_endpoint(url, instruments, inner=1):
    """Query all instruments in one FinnGen endpoint; return list of MR rows."""
    rt = RemoteTabix(url)
    out = []
    if inner <= 1:
        for inst in instruments:
            r = _mr_one_snp(rt, inst)
            if r:
                out.append(r)
    else:
        with ThreadPoolExecutor(max_workers=inner) as ex:
            for r in ex.map(lambda i: _mr_one_snp(rt, i), instruments):
                if r:
                    out.append(r)
    return out


_print_lock = threading.Lock()


def worker(ep, instruments, i, total):
    code = ep["phenocode"]
    dest = os.path.join(PART, code + ".tsv")
    if os.path.exists(dest):
        return code, -1  # already done
    t0 = time.time()
    try:
        rows = mr_for_endpoint(ep["path_https"], instruments, inner=INNER_WORKERS)
    except Exception as e:  # noqa: BLE001
        with _print_lock:
            print(f"[{i}/{total}] {code:28s} FAILED {type(e).__name__}", flush=True)
        return code, -2
    for r in rows:
        r["phenocode"] = code
        r["phenotype"] = ep["phenotype"]
        r["num_cases"] = ep["num_cases"]
        r["num_controls"] = ep["num_controls"]
        r["category"] = ep["category"]
    pd.DataFrame(rows).to_csv(dest, sep="\t", index=False)
    with _print_lock:
        print(f"[{i}/{total}] {code:28s} hits={len(rows):4d}  {time.time()-t0:5.1f}s", flush=True)
    return code, len(rows)


def main():
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])

    instruments = load_instruments()
    eps = list(csv.DictReader(open(MAN, encoding="utf-8"), delimiter="\t"))
    if "--reverse" in sys.argv:
        eps = eps[::-1]
    if limit:
        eps = eps[:limit]
    total = len(eps)
    print(f"[41] instruments={len(instruments)} | endpoints={total} | workers={N_WORKERS}", flush=True)

    done = 0
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(worker, ep, instruments, i, total): ep
                for i, ep in enumerate(eps, 1)}
        for fut in as_completed(futs):
            fut.result()
            done += 1

    # aggregate
    frames = []
    for fn in os.listdir(PART):
        if fn.endswith(".tsv"):
            try:
                df = pd.read_csv(os.path.join(PART, fn), sep="\t")
                if len(df):
                    frames.append(df)
            except Exception:
                pass
    if not frames:
        print("[41] no results")
        return
    res = pd.concat(frames, ignore_index=True).sort_values("MR_p").reset_index(drop=True)
    m = len(res)
    res["FDR"] = (res["MR_p"] * m / (res.index + 1)).clip(upper=1.0)
    res["FDR"] = res["FDR"][::-1].cummin()[::-1]
    ann = pd.read_csv(os.path.join(ROOT, "02_data_processed",
                                   "plasma_immune_protein_annotation.tsv"), sep="\t")
    res = res.merge(ann[["gene_symbol", "immune_class"]], on="gene_symbol", how="left")
    out = os.path.join(GEN, "cis_MR_ALL_finngen_results.tsv")
    res.to_csv(out, sep="\t", index=False)
    print(f"\n[41] DONE: {res['phenocode'].nunique()} endpoints with data | "
          f"{m:,} tests | FDR<0.05: {int((res.FDR < 0.05).sum())} | "
          f"genes: {res['gene_symbol'].nunique()}")
    print("wrote", out)


if __name__ == "__main__":
    main()
