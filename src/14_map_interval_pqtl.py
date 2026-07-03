#!/usr/bin/env python
"""
HDDM Layer 54 - Step 14
Build a gene_symbol -> INTERVAL (Sun 2018) pQTL study-accession map from the
GWAS Catalog, so protein-level cis-pQTL instruments can replace / confirm the
transcript-level cis-eQTL instruments. Fully public (no login, no form).

INTERVAL traits are named like  "Beta-defensin 119 levels (DEFB119.10689.5.3)"
-> the token inside the parentheses starts with the SomaScan target gene symbol.
We also map via the UniProt/gene token before the first dot.
Output: 06_genetic_causality/interval_pqtl_study_map.tsv
"""
import os, time, json, urllib.request, re
import pandas as pd

ROOT = r"I:\Plasma immune atalas"
GEN  = os.path.join(ROOT, "06_genetic_causality")
PMID = 29875488  # Sun et al. 2018 Nature, INTERVAL plasma proteome

def fetch(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            print("  retry", i, e); time.sleep(3)
    raise RuntimeError("failed "+url)

base = ("https://www.ebi.ac.uk/gwas/rest/api/studies/search/"
        f"findByPublicationIdPubmedId?pubmedId={PMID}&size=500")
rows = []
page = 0
while True:
    d = fetch(base + f"&page={page}")
    studies = d.get("_embedded", {}).get("studies", [])
    if not studies: break
    for s in studies:
        acc = s["accessionId"]
        trait = (s.get("diseaseTrait") or {}).get("trait", "")
        full = s.get("fullPvalueSet", False)
        # token inside parentheses e.g. (DEFB119.10689.5.3)
        m = re.search(r"\(([^)]+)\)", trait)
        tok = m.group(1) if m else ""
        gene = tok.split(".")[0] if tok else ""
        rows.append(dict(accession=acc, trait=trait, soma_target=tok,
                         gene_token=gene, full_pvalue=full))
    print(f"page {page}: {len(studies)} studies (total {len(rows)})")
    tp = d.get("page", {}).get("totalPages", 1)
    page += 1
    if page >= tp: break

m = pd.DataFrame(rows)
m.to_csv(os.path.join(GEN, "interval_pqtl_study_map.tsv"), sep="\t", index=False)
print("total INTERVAL studies:", len(m), "| with full sumstats:", int(m.full_pvalue.sum()))

# match to our significant-hit genes
tier = pd.read_csv(os.path.join(GEN, "evidence_tiered_targets.tsv"), sep="\t")
hit_genes = sorted(tier.gene_symbol.unique())
hit = m[m.gene_token.str.upper().isin([g.upper() for g in hit_genes])]
print("\nsig-hit genes with an INTERVAL pQTL study:",
      sorted(hit.gene_token.str.upper().unique()))
missing = set(g.upper() for g in hit_genes) - set(hit.gene_token.str.upper())
print("NOT found in INTERVAL:", sorted(missing))
hit.to_csv(os.path.join(GEN, "interval_pqtl_hit_studies.tsv"), sep="\t", index=False)
print("wrote map + hit-study list")
