#!/usr/bin/env python
"""
HDDM Layer 54 - Step 23
Download an EXPANDED, well-powered NON-immune FinnGen R12 phenome so the plasma
immunome cis-MR can be extended beyond the 13 autoimmune diseases into
cardiovascular / metabolic / renal / neuro / aging endpoints. All public, from
the local manifest (path_https). Each ~700-800 MB.

Selected endpoints (phenocode -> readable disease, curated for case count > ~4k):
"""
import os, sys, csv, urllib.request

ROOT = r"I:\Plasma immune atalas"
MAN  = os.path.join(ROOT, "01_data_raw", "FinnGen_GWAS", "finngen_R12_manifest.tsv")
DEST = os.path.join(ROOT, "01_data_raw", "FinnGen_GWAS", "sumstats")
os.makedirs(DEST, exist_ok=True)

# curated expansion phenome (non-immune; complements the 13 autoimmune already present)
WANT = {
 "I9_HYPTENS":"Hypertension", "T2D":"Type 2 diabetes", "I9_AF":"Atrial fibrillation",
 "I9_CHD":"Coronary heart disease", "I9_HEARTFAIL":"Heart failure",
 "I9_STR":"Stroke", "E4_OBESITY":"Obesity", "I9_VTE":"Venous thromboembolism",
 "F5_DEMENTIA":"Dementia", "G6_EPLEPSY":"Epilepsy",
 "N14_CHRONKIDNEYDIS":"Chronic kidney disease", "G6_AD_WIDE":"Alzheimer's disease",
 "M13_OSTEOPOROSIS":"Osteoporosis", "T1D":"Type 1 diabetes",
 "H7_GLAUCOMA":"Glaucoma",
}

rows = list(csv.DictReader(open(MAN, encoding="utf-8"), delimiter="\t"))
by_code = {r["phenocode"]: r for r in rows}
sel = [by_code[c] for c in WANT if c in by_code]
missing = [c for c in WANT if c not in by_code]
if missing:
    print("WARNING: phenocodes not in manifest:", missing)
print(f"Expanded phenome: {len(sel)} endpoints to fetch")

for i, r in enumerate(sel, 1):
    url = r["path_https"]; fn = os.path.join(DEST, os.path.basename(url))
    if os.path.exists(fn) and os.path.getsize(fn) > 100_000_000:
        print(f"[{i}/{len(sel)}] skip {os.path.basename(fn)} (exists)", flush=True); continue
    try:
        print(f"[{i}/{len(sel)}] {r['phenocode']} ({WANT[r['phenocode']]}) ...", end=" ", flush=True)
        urllib.request.urlretrieve(url, fn)
        print(f"{os.path.getsize(fn)//1_000_000} MB", flush=True)
    except Exception as e:
        print("FAILED", e, flush=True)
print("Done ->", DEST, flush=True)
