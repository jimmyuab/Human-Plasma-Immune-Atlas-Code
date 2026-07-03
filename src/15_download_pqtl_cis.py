#!/usr/bin/env python
"""
HDDM Layer 54 - Step 15
Download INTERVAL (Sun 2018) harmonised plasma pQTL summary stats for the
significant-hit genes and extract the cis window (+/-500 kb of the gene body)
into small per-aptamer tables. Fully public GWAS Catalog FTP (no login/form).

For each aptamer study:
  - build FTP path from accession (GCSTnnn ranges of 1000)
  - stream-download the harmonised .h.tsv.gz (~300 MB) if not cached
  - keep only rows on the gene's chromosome within +/-500 kb of the gene body
  - write cis table with real beta + standard_error (no Zhu needed)
Output: 01_data_raw/INTERVAL_pQTL/cis/<gene>__<accession>.tsv
"""
import os, json, gzip, urllib.request, shutil
import pandas as pd

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality")
RAW  = os.path.join(ROOT, "01_data_raw", "INTERVAL_pQTL")
CIS  = os.path.join(RAW, "cis")
os.makedirs(CIS, exist_ok=True)

WIN = 500_000
FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

coords = json.load(open(os.path.join(GEN, "interval_gene_coords_grch37.json")))
hit = pd.read_csv(os.path.join(GEN, "interval_pqtl_hit_studies.tsv"), sep="\t")

def ftp_dir(acc):
    n = int(acc.replace("GCST", ""))
    lo = ((n - 1) // 1000) * 1000 + 1
    hi = lo + 999
    return f"{FTP}/GCST{lo:08d}-GCST{hi:08d}/{acc}"

def download(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1_000_000:
        return dest
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f, length=1024*1024)
    os.replace(tmp, dest)
    return dest

summary = []
for _, row in hit.iterrows():
    g = row.gene_token.upper()
    if g not in coords:
        # coords keyed by original symbol; try direct
        alt = [k for k in coords if k.upper()==g]
        if not alt:
            print("no coords for", g); continue
        g_key = alt[0]
    else:
        g_key = g
    chrom, start, end, strand = coords[g_key]
    acc = row.accession
    d = ftp_dir(acc)
    fname = f"{acc}.h.tsv.gz"
    url = f"{d}/harmonised/{fname}"
    raw = os.path.join(RAW, fname)
    cis_out = os.path.join(CIS, f"{g}__{acc}.tsv")
    if os.path.exists(cis_out):
        n = sum(1 for _ in open(cis_out)) - 1
        print(f"{g} {acc}: cached cis ({n} SNPs)"); summary.append((g,acc,n)); continue
    print(f"{g} {acc}: downloading {url} ...")
    try:
        download(url, raw)
    except Exception as e:
        print("  DL FAIL", e); continue
    lo, hi = start - WIN, end + WIN
    keep = []
    with gzip.open(raw, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h:i for i,h in enumerate(hdr)}
        ci, pi = ix["chromosome"], ix["base_pair_location"]
        for line in f:
            p = line.rstrip("\n").split("\t")
            try:
                if p[ci] != str(chrom): continue
                bp = int(p[pi])
            except (ValueError, IndexError): continue
            if lo <= bp <= hi:
                keep.append(line)
    with open(cis_out, "w") as o:
        o.write("\t".join(hdr) + "\n"); o.writelines(keep)
    print(f"  {g}: {len(keep)} cis SNPs -> {os.path.basename(cis_out)}")
    summary.append((g, acc, len(keep)))
    # remove the big genome-wide file to save space
    try: os.remove(raw)
    except OSError: pass

print("\n=== cis extraction summary ===")
for g,a,n in summary: print(f"  {g:10s} {a}  {n} cis SNPs")
print("done:", len(summary), "aptamers")
