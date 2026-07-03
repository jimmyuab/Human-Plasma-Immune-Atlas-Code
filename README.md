# An open, genetics-anchored plasma immunome atlas

**Recovers established therapeutic immune axes and nominates colocalized, replicated autoimmune targets — entirely from public data.**

This repository is the fully reproducible, open resource behind the manuscript
*"An open, genetics-anchored plasma immunome atlas recovers established therapeutic immune axes and nominates colocalized autoimmune targets."*
It maps the Olink plasma immune proteome to disease through a calibrated, four-layer causal-evidence arc and a novelty/translation layer — with **no credentialed or individual-level data** required.

---

## Headline results

- **1,007 plasma immune proteins** curated from the 2,923-analyte Olink Explore universe (HPA + MSigDB annotation).
- **cis-eQTL Mendelian randomization** (eQTLGen → 13 FinnGen R12 immune diseases): 8,749 tests, **32 gene–disease hits at FDR < 0.05**.
- Recovers **known drug-target biology from genetics alone**: IL6ST→RA (tocilizumab), CTLA4→autoimmune thyroid/RA (abatacept), TNFRSF1A→ankylosing spondylitis (etanercept), IL4→psoriasis (dupilumab) — and recovers the *correct pharmacological direction* of each.
- **Transcript colocalization** (coloc.abf, per-SNP Wakefield ABF): 17 loci PP.H4 ≥ 0.8.
- **Protein-level pQTL** validation (INTERVAL plasma cis-pQTLs): **2 tier-5 protein-level causal targets** — TNFSF14→multiple sclerosis and SWAP70→rheumatoid arthritis (colocalized at both transcript and protein).
- **Independent replication** (OpenGWAS, non-FinnGen consortium GWAS): 19/19 covered hits directionally concordant, **17 replicated** (P < 0.05), incl. TNFSF14→MS P=1.5×10⁻¹².
- **Novelty layer**: immune-class enrichment, cross-disease pleiotropy map, genetic support for therapeutic direction, and a novel-vs-known target prioritization.

Every claim is bound to an explicit **evidence tier** (T2 MHC-caution → T5 protein-level causal); MHC signals are held at nomination.

---

## The evidence arc

```
discovery              transcript            protein-level          independent           translation
cis-eQTL MR    ─►      colocalization  ─►    pQTL (INTERVAL)  ─►    replication     ─►    novelty /
(eQTLGen →              (coloc.abf,           MR + coloc            (OpenGWAS,            drug-direction
 FinnGen R12)           PP.H4 ≥ 0.8)                                non-FinnGen)          & prioritization
```

## Figures

- `figures/main/` — 9 multi-panel main figures (Figure1–9).
- `figures/supplementary/` — 70 Extended Data / Supplementary figures (per-disease MR volcano, QQ + genomic inflation λ, forest plots, per-gene INTERVAL cis-pQTL regional plots, discovery-vs-replication).

## Manuscript

- `manuscript/Plasma_Immunome_Phenome_Atlas_Nature.docx` — 9 figures, 7 results sections, evidence-tiered Table 1.

---

## Reproduce the pipeline

```bash
pip install -r requirements.txt
# steps are numbered; run in order
python src/01_build_plasma_immune_annotation.py
python src/06_build_mr_instruments.py
python src/07_run_cis_mr.py
python src/11_run_coloc.py
python src/12_claim_gate.py
python src/16_run_pqtl_mr_coloc.py
python src/17_integrate_pqtl_tiers.py
python src/18_opengwas_replication.py      # needs a free OpenGWAS JWT (see below)
python src/19_integrate_replication.py
python src/20_supplementary_figures.py     # 70 supplementary figures
python src/21_novelty_analysis.py          # novelty layer + Figure 9
python src/13_write_manuscript.py          # assembles the .docx
```

Scripts expect the project root layout (`01_data_raw/`, `02_data_processed/`, `06_genetic_causality/`, `08_figures/`). The pre-computed **result tables are included** under `results/` so figures and the manuscript can be regenerated without re-downloading multi-GB sumstats.

## Training a predictive model (PIRS)

The atlas itself is a **summary-statistics causal resource**, not a per-individual training set. A supervised **Plasma Immune Risk Score (PIRS)** — predicting future disease from a person's Olink profile — requires **individual-level UK Biobank Olink NPX + phenotypes**, which are **controlled-access** and deliberately **not** shipped here.

`src/22_train_pirs.py` is a ready-to-run, honest-by-design scaffold:

- With **no UKB data present** it prints the exact input schema it needs and exits cleanly — it never fabricates data.
- Once approved UKB data is placed under `01_data_raw/UKB_Olink_NPX/` (`npx_matrix.tsv`, `outcomes.tsv`, optional `covariates.tsv`), one command trains a **cross-validated elastic-net survival model** (scikit-survival Coxnet → lifelines → logistic fallback), restricted to the 1,007 immune proteins, and writes per-protein PIRS weights, a CV C-index/AUROC table, fitted model pickles, and a performance figure.

```bash
python src/22_train_pirs.py     # prints schema + exits if UKB data absent
```

To obtain the data, apply via the UK Biobank Access Management System (Olink Field 30900 + disease phenotypes).

## Data sources (all public)

| Layer | Source | Access |
|---|---|---|
| Protein universe | Olink Explore (UKB coding 143) | public |
| Annotation | Human Protein Atlas, MSigDB C7/C8 | public |
| eQTL instruments | eQTLGen cis-eQTLs | public |
| Disease GWAS (discovery) | FinnGen R12 | public |
| Protein pQTL | INTERVAL (Sun 2018) via EBI GWAS Catalog FTP | public |
| Replication GWAS | IMSGC, Okada, IGAS, Stuart, Fischer, Sakaue via **OpenGWAS** | free JWT token |

> **OpenGWAS token:** replication (`src/18`) needs a free personal JWT from https://api.opengwas.io — register, place the token in a local file, and never commit it. The included `results/genetic_causality/opengwas_replication.tsv` already contains the saved replication output, so re-running `src/18` is optional.

## Not included (gated, by design)

Individual-level UK Biobank Olink NPX + phenotypes (controlled access) and UKB-PPP / deCODE pQTL (Synapse / DUA-gated) are **not** in this repository and are **not required** — the atlas is built to be fully reproducible from open data.

---

## One-click publish

To publish this folder as a public GitHub repository and cut a release:

- **Windows:** double-click `publish.bat`
- **macOS/Linux:** `bash publish.sh`

The script initializes git, commits, creates the GitHub repo, pushes, and tags a `v1.0.0` release. See the top of the script for the two variables (repo name / visibility) you can edit. Requires the [GitHub CLI](https://cli.github.com) (`gh`) authenticated, or an existing `origin` remote.

## License

- **Code** (`src/`): MIT — see `LICENSE`.
- **Data tables & figures** (`results/`, `figures/`): CC-BY-4.0.

## Citation

See `CITATION.cff`.
