#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
49 - PAN-PHENOME protein-level pQTL MR + colocalization
=======================================================
Lifts the protein (INTERVAL Sun 2018 plasma pQTL) layer from the 13-disease
autoimmune core to EVERY pan-phenome cis-MR hit, so the pQTL layer matches the
coverage of the cis-MR (2,465 diseases) and coloc (all hits) layers.

Statistics are IDENTICAL to src/16 (Wald-ratio MR with real beta+SE on both
sides; per-SNP-variance Wakefield coloc.abf, same priors and W). The only
change is that the disease cis-window is fetched by REMOTE TABIX region query
(src/remote_tabix.py) rather than from a locally downloaded sumstat file, so it
scales to the whole phenome.

Exposure : 01_data_raw/INTERVAL_pQTL/cis/<GENE>__<ACC>.tsv   (GRCh37, real SE)
Outcome  : FinnGen R12 public sumstats, remote byte-range   (GRCh38)
Matching : by rsid; the query window is the aptamer cis window lifted hg19->hg38

Outputs:
  06_genetic_causality/pqtl_MR_ALL_finngen_results.tsv
  06_genetic_causality/pqtl_coloc_ALL_finngen_results.tsv

Run:  python src/49_pqtl_all_phenome_mr_coloc.py [--limit N]
"""
import os
import sys
import math
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remote_tabix import RemoteTabix  # noqa: E402
from pyliftover import LiftOver  # noqa: E402

ROOT = r"I:\Plasma immune atalas"
RAW = os.path.join(ROOT, "01_data_raw")
GEN = os.path.join(ROOT, "06_genetic_causality")
CIS = os.path.join(RAW, "INTERVAL_pQTL", "cis")
CHAIN = os.path.join(RAW, "liftover", "hg19ToHg38.over.chain.gz")
FG_URL = ("https://storage.googleapis.com/finngen-public-data-r12/"
          "summary_stats/release/finngen_R12_{code}.gz")

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
W_PQTL = 0.15 ** 2
W_GWAS = 0.20 ** 2
N_WORKERS = 10
FG_REF, FG_ALT, FG_RSID = 2, 3, 4
FG_BETA, FG_SEBETA, FG_AF = 8, 9, 12


def ambiguous(a, b):
    return COMP.get(a) == b


def lABF(Z, V, W):
    r = W / (W + V)
    return 0.5 * (math.log(1 - r) + r * Z * Z)


def logsumexp(a):
    m = np.max(a)
    return m + math.log(np.sum(np.exp(a - m)))


def load_cis(fp):
    df = pd.read_csv(fp, sep="\t", dtype=str, low_memory=False)
    df = df.rename(columns={"hm_rsid": "rsid"})
    for c in ["beta", "standard_error", "p_value", "base_pair_location"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["rsid", "beta", "standard_error", "p_value",
                           "effect_allele", "other_allele", "base_pair_location"])
    df = df[df.standard_error > 0]
    df["effect_allele"] = df["effect_allele"].str.upper()
    df["other_allele"] = df["other_allele"].str.upper()
    df = df[df.rsid.str.startswith("rs")]
    df["chromosome"] = df["chromosome"].astype(str)
    return df[["rsid", "chromosome", "base_pair_location", "effect_allele",
               "other_allele", "beta", "standard_error", "p_value"]]


def fetch_fg(rt, chrom, p_lo, p_hi):
    try:
        rows = rt.region(str(chrom), p_lo, p_hi)
    except Exception:
        return {}
    out = {}
    for f in rows:
        if len(f) <= FG_AF:
            continue
        try:
            b = float(f[FG_BETA]); se = float(f[FG_SEBETA]); af = float(f[FG_AF])
        except (ValueError, IndexError):
            continue
        if se <= 0 or not (0 < af < 1):
            continue
        rec = (b, se, f[FG_REF].upper(), f[FG_ALT].upper(), af)
        for r in f[FG_RSID].split(","):
            out[r] = rec
    return out


def run_disease(code, disease, aptamers, lo):
    """aptamers: list of (gene, acc, cis_df). Returns (mr_rows, coloc_rows)."""
    try:
        rt = RemoteTabix(FG_URL.format(code=code))
    except Exception:
        return [], []
    mr_rows, cl_rows = [], []
    for gene, acc, df in aptamers:
        chrom = df.chromosome.iloc[0]
        b19 = lo.convert_coordinate("chr" + str(chrom), int(df.base_pair_location.min()))
        e19 = lo.convert_coordinate("chr" + str(chrom), int(df.base_pair_location.max()))
        if not b19 or not e19:
            continue
        p_lo = min(b19[0][1], e19[0][1]) - 1000
        p_hi = max(b19[0][1], e19[0][1]) + 1000
        fg = fetch_fg(rt, chrom, p_lo, p_hi)
        if not fg:
            continue

        # ---- MR: strongest cis-pQTL present in FinnGen, non-ambiguous ----
        cand = df[df.rsid.isin(fg.keys())].sort_values("p_value")
        for _, inst in cand.iterrows():
            if ambiguous(inst.effect_allele, inst.other_allele):
                continue
            b_o, se_o, ref, alt, af = fg[inst.rsid]
            if alt == inst.effect_allele and ref == inst.other_allele:
                sign = 1.0
            elif ref == inst.effect_allele and alt == inst.other_allele:
                sign = -1.0
            else:
                continue
            b_exp = float(inst.beta)
            if b_exp == 0:
                continue
            b_out = b_o * sign
            mr_beta = b_out / b_exp
            mr_se = abs(se_o / b_exp)
            z = mr_beta / mr_se
            p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
            mr_rows.append(dict(
                gene=gene, acc=acc, disease=disease, disease_code=code,
                chr=chrom, rsid=inst.rsid, pqtl_p=float(inst.p_value),
                b_exp=b_exp, b_out=b_out, se_out=se_o,
                MR_beta=mr_beta, MR_se=mr_se, MR_p=p, OR=math.exp(mr_beta),
                OR_l95=math.exp(mr_beta - 1.96 * mr_se),
                OR_u95=math.exp(mr_beta + 1.96 * mr_se)))
            break

        # ---- coloc over shared cis SNPs ----
        labf = []
        for rs, ea, oa, b, se in zip(df.rsid, df.effect_allele, df.other_allele,
                                     df.beta, df.standard_error):
            r = fg.get(rs)
            if r is None or ambiguous(ea, oa):
                continue
            b_o, se_o, ref, alt, af = r
            if not ((alt == ea and ref == oa) or (ref == ea and alt == oa)):
                continue
            labf.append((lABF(b / se, se * se, W_PQTL),
                         lABF(b_o / se_o, se_o * se_o, W_GWAS)))
        rec = dict(gene=gene, acc=acc, disease=disease, disease_code=code,
                   nsnp=len(labf), PP_H0=np.nan, PP_H1=np.nan, PP_H2=np.nan,
                   PP_H3=np.nan, PP_H4=np.nan)
        if len(labf) >= 5:
            le = np.array([x[0] for x in labf]); lg = np.array([x[1] for x in labf])
            p1 = p2 = 1e-4; p12 = 1e-5
            l1 = logsumexp(le); l2 = logsumexp(lg); l4 = logsumexp(le + lg)
            l12 = l1 + l2; diff = 1.0 - math.exp(l4 - l12)
            lH3 = (-np.inf if diff <= 0 else
                   math.log(p1) + math.log(p2) + l12 + math.log(diff))
            arr = np.array([0.0, math.log(p1) + l1, math.log(p2) + l2,
                            lH3, math.log(p12) + l4])
            pp = np.exp(arr - logsumexp(arr))
            rec.update(PP_H0=pp[0], PP_H1=pp[1], PP_H2=pp[2],
                       PP_H3=pp[3], PP_H4=pp[4])
        cl_rows.append(rec)
    return mr_rows, cl_rows


def main():
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])

    mr = pd.read_csv(os.path.join(GEN, "cis_MR_ALL_finngen_results.tsv"), sep="\t",
                     usecols=["gene_symbol", "phenocode", "phenotype", "FDR"])
    hits = mr[mr.FDR < 0.05][["gene_symbol", "phenocode", "phenotype"]].drop_duplicates()
    hits["gene_u"] = hits.gene_symbol.astype(str).str.upper()

    files = {}
    for fp in sorted(glob.glob(os.path.join(CIS, "*.tsv"))):
        gene, acc = os.path.basename(fp)[:-4].split("__")
        files.setdefault(gene.upper(), []).append((acc, fp))
    print(f"[49] aptamer cis files: {sum(len(v) for v in files.values())} "
          f"over {len(files)} genes", flush=True)

    hits = hits[hits.gene_u.isin(files)]
    codes = sorted(hits.phenocode.unique())
    if limit:
        codes = codes[:limit]
    print(f"[49] protein-layer coverage: {len(hits)} gene-disease hits | "
          f"{hits.gene_u.nunique()} genes | {len(codes)} diseases", flush=True)

    cache = {}

    def cis_for(gene_u):
        if gene_u not in cache:
            cache[gene_u] = [(gene_u, acc, load_cis(fp)) for acc, fp in files[gene_u]]
        return cache[gene_u]

    for g in hits.gene_u.unique():
        cis_for(g)
    print(f"[49] loaded cis tables for {len(cache)} genes", flush=True)

    by_code = hits.groupby("phenocode")
    all_mr, all_cl = [], []
    done = 0

    def job(code):
        sub = by_code.get_group(code)
        apt = [a for g in sub.gene_u.unique() for a in cache[g]]
        return run_disease(code, sub.phenotype.iloc[0], apt, LO)

    global LO
    LO = LiftOver(CHAIN)

    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(job, c): c for c in codes}
        for fut in as_completed(futs):
            try:
                m, c = fut.result()
            except Exception as e:
                print(f"[49] disease {futs[fut]} FAILED: {e}", flush=True)
                continue
            all_mr.extend(m); all_cl.extend(c)
            done += 1
            if done % 10 == 0:
                print(f"[49] {done}/{len(codes)} diseases | MR={len(all_mr)} "
                      f"coloc={len(all_cl)}", flush=True)

    mrd = pd.DataFrame(all_mr)
    if len(mrd):
        mrd["FDR"] = multipletests(mrd.MR_p, method="fdr_bh")[1]
        mrd = mrd.sort_values("MR_p")
    mrd.to_csv(os.path.join(GEN, "pqtl_MR_ALL_finngen_results.tsv"), sep="\t", index=False)

    cld = pd.DataFrame(all_cl).sort_values("PP_H4", ascending=False)
    cld.to_csv(os.path.join(GEN, "pqtl_coloc_ALL_finngen_results.tsv"), sep="\t", index=False)

    print("\n[49] === PAN-PHENOME pQTL-MR (top) ===")
    if len(mrd):
        print(mrd.head(25)[["gene", "disease", "OR", "OR_l95", "OR_u95",
                            "MR_p", "FDR"]].to_string(index=False))
        print(f"[49] tests={len(mrd)} FDR<0.05={int((mrd.FDR < 0.05).sum())} "
              f"| diseases={mrd.disease.nunique()} genes={mrd.gene.nunique()}")
    print("\n[49] === PAN-PHENOME pQTL coloc PP.H4>=0.8 ===")
    strong = cld[cld.PP_H4 >= 0.8]
    print(strong.head(25)[["gene", "disease", "nsnp", "PP_H4"]].to_string(index=False))
    print(f"[49] coloc loci={len(cld)} strong={len(strong)} "
          f"| diseases={cld.disease.nunique()}")


if __name__ == "__main__":
    main()
