#!/usr/bin/env python
"""
HDDM Layer 54 - optional bulk downloader
Download FinnGen R12 summary statistics for IMMUNE / INFLAMMATORY / INFECTIOUS /
AUTOIMMUNE endpoints from the manifest (each ~50-150 MB). Full set is ~250 GB; this
filters to disease chapters most relevant to the plasma immunome for MR/coloc.
Run:  python src/05_download_finngen_immune.py [--all]
"""
import os, sys, csv, urllib.request

ROOT = r"I:\Plasma immune atalas"
MAN  = os.path.join(ROOT, "01_data_raw", "FinnGen_GWAS", "finngen_R12_manifest.tsv")
DEST = os.path.join(ROOT, "01_data_raw", "FinnGen_GWAS", "sumstats")
os.makedirs(DEST, exist_ok=True)

KEEP = ("immune", "inflamm", "infectious", "autoimmun", "rheumat", "AB1_", "III ",
        "blood and blood-forming")
download_all = "--all" in sys.argv

rows = list(csv.DictReader(open(MAN, encoding="utf-8"), delimiter="\t"))
sel = rows if download_all else [
    r for r in rows if any(k.lower() in (r["category"] + r["phenotype"]).lower() for k in KEEP)]
print(f"Manifest endpoints: {len(rows)} | selected: {len(sel)} "
      f"({'ALL' if download_all else 'immune-relevant'})")

for i, r in enumerate(sel, 1):
    url = r["path_https"]; fn = os.path.join(DEST, os.path.basename(url))
    if os.path.exists(fn) and os.path.getsize(fn) > 0:
        print(f"[{i}/{len(sel)}] skip {os.path.basename(fn)}"); continue
    try:
        print(f"[{i}/{len(sel)}] {r['phenocode']} ...", end=" ", flush=True)
        urllib.request.urlretrieve(url, fn)
        print(f"{os.path.getsize(fn)//1_000_000} MB")
    except Exception as e:
        print("FAILED", e)
print("Done ->", DEST)
