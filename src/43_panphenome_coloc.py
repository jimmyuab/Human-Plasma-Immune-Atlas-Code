#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
43 - REAL phenome-wide colocalization
=====================================
Runs coloc.abf for EVERY pan-phenome cis-MR hit (FDR<0.05, non-MHC) from
src/41 (cis_MR_ALL_finngen_results.tsv), across all FinnGen R12 diseases --
not just the 28 core. Uses the IDENTICAL corrected per-SNP Wakefield ABF
coloc.abf as src/11 / src/27 (same priors, same W). The only change is that
the disease cis-window is fetched by REMOTE TABIX region query
(src/remote_tabix.py) instead of scanning a full local file, so it scales to
the whole phenome without downloading ~2 TB.

Inputs:
  06_genetic_causality/cis_MR_ALL_finngen_results.tsv   (pan-phenome hits)
  01_data_raw/eQTLGen/cis-eQTLs-full.txt.gz             (eQTL cis windows)
  01_data_raw/liftover/hg19ToHg38.over.chain.gz         (eQTL hg19 -> FinnGen hg38)
Output:
  06_genetic_causality/coloc_ALL_finngen_results.tsv
  06_genetic_causality/coloc_eqtl_cis_cache_pos.pkl     (gene -> cis SNPs w/ pos)

Run:  python src/43_panphenome_coloc.py [--limit N]
"""
import os
import sys
import gzip
import math
import pickle
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remote_tabix import RemoteTabix  # noqa: E402
from pyliftover import LiftOver  # noqa: E402

ROOT = r"I:\Plasma immune atalas"
RAW = os.path.join(ROOT, "01_data_raw")
GEN = os.path.join(ROOT, "06_genetic_causality")
EQTL_FULL = os.path.join(RAW, "eQTLGen", "cis-eQTLs-full.txt.gz")
CHAIN = os.path.join(RAW, "liftover", "hg19ToHg38.over.chain.gz")
CACHE = os.path.join(GEN, "coloc_eqtl_cis_cache_pos.pkl")
MAN = os.path.join(RAW, "FinnGen_GWAS", "finngen_R12_manifest.tsv")
FG_URL = ("https://storage.googleapis.com/finngen-public-data-r12/"
          "summary_stats/release/finngen_R12_{code}.gz")

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
W_EQTL = 0.15 ** 2
W_GWAS = 0.20 ** 2
N_WORKERS = 12
# FinnGen columns
FG_CHROM, FG_POS, FG_REF, FG_ALT, FG_RSID = 0, 1, 2, 3, 4
FG_BETA, FG_SEBETA, FG_AF = 8, 9, 12


def ambiguous(a, b):
    return COMP.get(a) == b


def lABF(Z, V, W):
    r = W / (W + V)
    return 0.5 * (math.log(1 - r) + r * Z * Z)


def logsumexp(a):
    m = np.max(a)
    return m + math.log(np.sum(np.exp(a - m)))


def build_cache(genes):
    """gene -> {rsid: (Z, a1, a2, N, chr, pos_hg19)} for the requested genes."""
    eq = {}
    if os.path.exists(CACHE):
        eq = pickle.load(open(CACHE, "rb"))
    missing = set(genes) - set(eq)
    if not missing:
        return eq
    print(f"[43] scanning eQTLGen full for {len(missing)} genes ...", flush=True)
    add = defaultdict(dict)
    with gzip.open(EQTL_FULL, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        gi, si, zi = ix["GeneSymbol"], ix["SNP"], ix["Zscore"]
        a1i, a2i, ni = ix["AssessedAllele"], ix["OtherAllele"], ix["NrSamples"]
        ci, pi = ix["SNPChr"], ix["SNPPos"]
        nrow = 0
        for line in f:
            nrow += 1
            p = line.rstrip("\n").split("\t")
            if p[gi] not in missing:
                continue
            add[p[gi]][p[si]] = (float(p[zi]), p[a1i].upper(), p[a2i].upper(),
                                 float(p[ni]), p[ci], int(p[pi]))
        print(f"[43]   scanned {nrow:,} rows", flush=True)
    for g, v in add.items():
        eq[g] = v
    pickle.dump(eq, open(CACHE, "wb"))
    return eq


def coloc_one(code, disease, genes_here, eq, lo):
    """Fetch each gene's cis-window from one disease (remote) and run coloc."""
    url = FG_URL.format(code=code)
    try:
        rt = RemoteTabix(url)
    except Exception:
        return []
    res = []
    for g in genes_here:
        if g not in eq:
            continue
        snps = eq[g]
        chrom = next(iter(snps.values()))[4]
        poss19 = [v[5] for v in snps.values()]
        # liftover window bounds hg19 -> hg38
        lo_b = lo.convert_coordinate("chr" + str(chrom), min(poss19))
        hi_b = lo.convert_coordinate("chr" + str(chrom), max(poss19))
        if not lo_b or not hi_b:
            continue
        p_lo = min(lo_b[0][1], hi_b[0][1])
        p_hi = max(lo_b[0][1], hi_b[0][1])
        try:
            rows = rt.region(str(chrom), p_lo - 1000, p_hi + 1000)
        except Exception:
            continue
        out = {}
        for f in rows:
            if len(f) <= FG_AF:
                continue
            rs = f[FG_RSID]
            try:
                b = float(f[FG_BETA]); se = float(f[FG_SEBETA]); af = float(f[FG_AF])
            except (ValueError, IndexError):
                continue
            if se <= 0 or not (0 < af < 1):
                continue
            # rsids field may be comma-joined
            for r in rs.split(","):
                out[r] = (b, se, f[FG_REF].upper(), f[FG_ALT].upper(), af)
        labf = []
        for rs, (ze, a1, a2, N, _c, _p) in snps.items():
            if rs not in out:
                continue
            b_o, se_o, ref, alt, af = out[rs]
            if ambiguous(a1, a2):
                continue
            if alt == a1 and ref == a2:
                sign = 1.0
            elif ref == a1 and alt == a2:
                sign = -1.0
            else:
                continue
            Zg = b_o / se_o; Vg = se_o * se_o
            denom = 2.0 * af * (1.0 - af) * (N + ze * ze)
            if denom <= 0:
                continue
            se_e = 1.0 / math.sqrt(denom); Ve = se_e * se_e; Ze = ze * sign
            labf.append((lABF(Ze, Ve, W_EQTL), lABF(Zg, Vg, W_GWAS)))
        nsnp = len(labf)
        rec = dict(gene=g, disease=disease, disease_code=code, nsnp=nsnp,
                   PP_H0=np.nan, PP_H1=np.nan, PP_H2=np.nan, PP_H3=np.nan, PP_H4=np.nan)
        if nsnp >= 5:
            le = np.array([x[0] for x in labf]); lg = np.array([x[1] for x in labf])
            p1 = p2 = 1e-4; p12 = 1e-5
            l1 = logsumexp(le); l2 = logsumexp(lg); l4 = logsumexp(le + lg)
            l12 = l1 + l2; diff = 1.0 - math.exp(l4 - l12)
            lH3 = (-np.inf if diff <= 0 else math.log(p1) + math.log(p2) + l12 + math.log(diff))
            arr = np.array([0.0, math.log(p1) + l1, math.log(p2) + l2, lH3, math.log(p12) + l4])
            pp = np.exp(arr - logsumexp(arr))
            rec.update(PP_H0=pp[0], PP_H1=pp[1], PP_H2=pp[2], PP_H3=pp[3], PP_H4=pp[4])
        res.append(rec)
    return res


def main():
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])

    mr = pd.read_csv(os.path.join(GEN, "cis_MR_ALL_finngen_results.tsv"), sep="\t")
    sig = mr[(mr.FDR < 0.05) & (mr.chr != 6)].copy()  # non-MHC
    hits = sig[["gene_symbol", "phenocode", "phenotype"]].drop_duplicates()
    by_code = hits.groupby("phenocode")
    codes = list(by_code.groups.keys())
    if limit:
        codes = codes[:limit]
    print(f"[43] non-MHC pan-phenome coloc: {len(hits)} hits | "
          f"{hits.gene_symbol.nunique()} genes | {len(codes)} diseases", flush=True)

    eq = build_cache(set(hits.gene_symbol))
    print(f"[43] eQTL cis cache: {len(eq)} genes", flush=True)
    lo = LiftOver(CHAIN)

    all_rows = []
    done = 0

    def run_code(code):
        sub = by_code.get_group(code)
        disease = sub.phenotype.iloc[0]
        genes_here = list(sub.gene_symbol.unique())
        return coloc_one(code, disease, genes_here, eq, lo)

    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(run_code, c): c for c in codes}
        for fut in as_completed(futs):
            rows = fut.result()
            all_rows.extend(rows)
            done += 1
            if done % 20 == 0:
                print(f"[43] {done}/{len(codes)} diseases done | rows={len(all_rows)}", flush=True)

    cd = pd.DataFrame(all_rows).sort_values("PP_H4", ascending=False)
    out = os.path.join(GEN, "coloc_ALL_finngen_results.tsv")
    cd.to_csv(out, sep="\t", index=False)
    print("\n[43] === PAN-PHENOME COLOC (all diseases) ===")
    print(cd[["gene", "disease", "nsnp", "PP_H4", "PP_H3"]].head(25).to_string(index=False))
    print(f"\n[43] PP.H4>=0.8: {int((cd.PP_H4 >= 0.8).sum())} | "
          f">=0.5: {int((cd.PP_H4 >= 0.5).sum())} | rows: {len(cd)}")
    print("[43] wrote", out)


if __name__ == "__main__":
    main()
