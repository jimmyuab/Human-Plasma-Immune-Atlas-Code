#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
48 - Expand the protein-level (INTERVAL pQTL) layer from the 13-disease
     autoimmune core to the WHOLE FinnGen phenome hit set.

The cis-MR and colocalization layers already cover all 2,465 FinnGen R12
endpoints (1,016 hits / 380 diseases / 220 proteins). The pQTL layer was the
only layer still restricted to the autoimmune core. This script selects every
pan-phenome hit gene that has a public INTERVAL (Sun 2018) aptamer, resolves
GRCh37 gene coordinates, downloads the harmonised genome-wide sumstats from
the EBI GWAS Catalog FTP, extracts the +/-500 kb cis window and deletes the
genome-wide file again (bounded disk use).

Inputs :  06_genetic_causality/cis_MR_ALL_finngen_results.tsv
          06_genetic_causality/interval_pqtl_study_map.tsv
Outputs:  06_genetic_causality/interval_pqtl_hit_studies_ALL.tsv
          06_genetic_causality/interval_gene_coords_grch37.json   (extended)
          01_data_raw/INTERVAL_pQTL/cis/<GENE>__<ACC>.tsv

Run:  python src/48_expand_pqtl_all_phenome.py [--limit N]
"""
import os
import sys
import json
import gzip
import time
import shutil
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

N_DL = 4

ROOT = r"I:\Plasma immune atalas"
GEN = os.path.join(ROOT, "06_genetic_causality")
RAW = os.path.join(ROOT, "01_data_raw", "INTERVAL_pQTL")
CIS = os.path.join(RAW, "cis")
os.makedirs(CIS, exist_ok=True)

WIN = 500_000
FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
ENS = "https://grch37.rest.ensembl.org"
COORDS_F = os.path.join(GEN, "interval_gene_coords_grch37.json")


def ens_coords(sym):
    url = f"{ENS}/lookup/symbol/homo_sapiens/{sym}?content-type=application/json"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        d = json.load(r)
    chrom = str(d["seq_region_name"])
    if chrom not in [str(i) for i in range(1, 23)] + ["X", "Y"]:
        raise ValueError("non-primary contig " + chrom)
    return [chrom, int(d["start"]), int(d["end"]), int(d.get("strand", 1))]


def ftp_dir(acc):
    n = int(acc.replace("GCST", ""))
    lo = ((n - 1) // 1000) * 1000 + 1
    return f"{FTP}/GCST{lo:08d}-GCST{lo + 999:08d}/{acc}"


def download(url, dest, tries=5):
    """Resumable download; a stalled EBI connection is retried with HTTP Range."""
    tmp = dest + ".part"
    for i in range(tries):
        have = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        req = urllib.request.Request(url)
        if have:
            req.add_header("Range", f"bytes={have}-")
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                mode = "ab" if (have and r.status == 206) else "wb"
                if mode == "wb":
                    have = 0
                with open(tmp, mode) as f:
                    while True:
                        chunk = r.read(4 * 1024 * 1024)
                        if not chunk:
                            break
                        f.write(chunk)
            os.replace(tmp, dest)
            return True
        except Exception as e:
            print(f"    retry {i} ({os.path.basename(dest)}): {e}", flush=True)
            time.sleep(5)
    return False


def extract_cis(gz_path, chrom, lo, hi, out_path):
    keep = []
    with gzip.open(gz_path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        ci, pi = ix["chromosome"], ix["base_pair_location"]
        cs = str(chrom)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ci] != cs:
                continue
            try:
                pos = int(p[pi])
            except ValueError:
                continue
            if lo <= pos <= hi:
                keep.append(p)
    df = pd.DataFrame(keep, columns=hdr)
    df.to_csv(out_path, sep="\t", index=False)
    return len(df)


def main():
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])

    mr = pd.read_csv(os.path.join(GEN, "cis_MR_ALL_finngen_results.tsv"),
                     sep="\t", usecols=["gene_symbol", "FDR"])
    hit_genes = sorted(mr.loc[mr.FDR < 0.05, "gene_symbol"].astype(str).str.upper().unique())
    print(f"pan-phenome hit genes: {len(hit_genes)}", flush=True)

    smap = pd.read_csv(os.path.join(GEN, "interval_pqtl_study_map.tsv"), sep="\t")
    smap = smap[smap.full_pvalue.astype(str) == "True"].copy()
    smap["gene_token_u"] = smap.gene_token.astype(str).str.upper()
    sel = smap[smap.gene_token_u.isin(hit_genes)].copy()
    sel = sel.sort_values(["gene_token_u", "accession"])
    sel.to_csv(os.path.join(GEN, "interval_pqtl_hit_studies_ALL.tsv"), sep="\t", index=False)
    print(f"aptamer studies for hit genes: {len(sel)} across "
          f"{sel.gene_token_u.nunique()} genes", flush=True)

    coords = json.load(open(COORDS_F))
    coords_u = {k.upper(): v for k, v in coords.items()}
    need = sorted(set(sel.gene_token_u) - set(coords_u))
    print(f"resolving GRCh37 coords for {len(need)} new genes ...", flush=True)
    for g in need:
        try:
            coords[g] = ens_coords(g)
            coords_u[g] = coords[g]
            print(f"  {g} -> {coords[g]}", flush=True)
        except Exception as e:
            print(f"  {g} COORD FAIL: {e}", flush=True)
        time.sleep(0.12)
    json.dump(coords, open(COORDS_F, "w"), indent=1)

    todo = []
    for _, r in sel.iterrows():
        g, acc = r.gene_token_u, r.accession
        out = os.path.join(CIS, f"{g}__{acc}.tsv")
        if os.path.exists(out) and os.path.getsize(out) > 200:
            continue
        if g not in coords_u:
            continue
        todo.append((g, acc, out))
    if limit:
        todo = todo[:limit]
    print(f"to download: {len(todo)} aptamer files", flush=True)

    log = []
    done = [0]

    def one(item):
        g, acc, out = item
        chrom, start, end, _ = coords_u[g]
        url = f"{ftp_dir(acc)}/harmonised/{acc}.h.tsv.gz"
        raw = os.path.join(RAW, f"{acc}.h.tsv.gz")
        if not (os.path.exists(raw) and os.path.getsize(raw) > 1_000_000):
            if not download(url, raw):
                return dict(gene=g, accession=acc, n_cis=0, status="dl_fail")
        try:
            n = extract_cis(raw, chrom, start - WIN, end + WIN, out)
            return dict(gene=g, accession=acc, n_cis=n, status="ok")
        except Exception as e:
            return dict(gene=g, accession=acc, n_cis=0, status=f"extract_fail:{e}")
        finally:
            if os.path.exists(raw):
                os.remove(raw)

    with ThreadPoolExecutor(max_workers=N_DL) as ex:
        futs = {ex.submit(one, t): t for t in todo}
        for fut in as_completed(futs):
            r = fut.result()
            log.append(r)
            done[0] += 1
            print(f"[{done[0]}/{len(todo)}] {r['gene']} {r['accession']} "
                  f"{r['status']} cis={r['n_cis']}", flush=True)

    if log:
        pd.DataFrame(log).to_csv(os.path.join(GEN, "interval_pqtl_download_log.tsv"),
                                 sep="\t", index=False)
    have = sorted({f.split("__")[0] for f in os.listdir(CIS) if f.endswith(".tsv")})
    print(f"\nDONE. cis files now cover {len(have)} genes", flush=True)


if __name__ == "__main__":
    main()
