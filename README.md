# Human Plasma Immune Atlas

**An open, genetics-anchored causal map from the plasma immunome to the human disease phenome.**

This repository is the fully reproducible, open resource behind the manuscript
*"Human Plasma Immune Atlas: an open, genetics-anchored causal map of the plasma immunome across the human disease phenome."*
It maps the Olink plasma immune proteome to disease through a calibrated, **seven-layer** causal-evidence arc and a novelty/translation engine — with **no credentialed or individual-level data** required.

👉 **Live interactive app — click and run, no login, no account, no upload:**
https://huggingface.co/spaces/jianlizhao/Human-Plasma-Immune-Atlas

---

## Headline results

- **1,007 plasma immune proteins** curated from the 2,923-analyte Olink Explore universe (HPA + MSigDB annotation).
- **Phenome-wide cis-eQTL Mendelian randomization**: **672 instrumentable plasma immune proteins × 2,466 FinnGen R12 endpoints = 1,656,872 MR tests**, of which **1,016 gene–disease pairs are causal at FDR < 5 %** (220 proteins, 379 diseases). *No curated disease list is used anywhere — every layer runs across the whole phenome.*
- **Bayesian colocalization** (`coloc.abf`, per-SNP Wakefield ABF): 900 loci tested, **417 colocalize at PP.H4 ≥ 0.8**.
- **Protein-level pQTL validation** (INTERVAL plasma cis-pQTLs, Sun 2018): **473 protein-level MR tests across 219 diseases and 89 proteins, 149 significant at FDR < 5 %**, plus **595 protein colocalization loci**.
- **Two-population replication (Finland → England)**: **189 like-for-like FinnGen-vs-UK-Biobank tests over 62 diseases and 90 genes; 162 (86 %) agree in causal direction and 89 replicate at P < 0.05** (133 exact-code + 56 name-matched endpoints). Only *single-cohort* UK Biobank datasets are eligible — several large public meta-analyses of the same endpoints silently include FinnGen, which would make "replication" circular.
- **Immune-cell / inflammation layer**: the same instruments run against blood cell counts, CRP and cytokines, to show how each protein perturbs the immune system itself.
- **Novelty & therapeutic-direction engine**: causal strength + coloc + pleiotropy + druggability + protein confirmation + replication, penalised for known drug axes and the MHC region. Classifies all 1,016 causal pairs into **329 NOVEL colocalized**, **35 NOVEL protein-confirmed**, **447 novel nominations**, **90 recovered known drug axes** and **115 MHC-caution**.
- **Recovers known drug-target biology from genetics alone** — IL6ST→RA (tocilizumab), CTLA4→autoimmune thyroid/RA (abatacept), TNFRSF1A→ankylosing spondylitis (etanercept), IL4→psoriasis (dupilumab) — and recovers the *correct pharmacological direction* of each.
- **Disease intelligence layer** (`src/31`, `src/50`): the **Final Required Output Table** (`results/genetic_causality/intelligence_layer_final_table_ALL.tsv`, **1,016 gene–disease pairs × 33 columns**): plasma detectability, PIRS coefficient, MR/coloc/pQTL/replication support, a Plasma Immune Novelty Score, a 1–5 novelty tier, therapeutic direction, and the best next validation experiment. Tier distribution: **5 protein-level novel targets, 359 prioritized transcript-level targets, 447 causal nominations, 115 MHC-held, 90 known-drug positive controls.**

Every claim is bound to an explicit **evidence tier** (T2 MHC-caution → T5 protein-level causal); MHC signals are held at nomination.

---

## The seven evidence layers

```
1 cis-eQTL MR        2 colocalization    3 protein pQTL      4 independent
(eQTLGen n≈31,684 ─► (coloc.abf,      ─► MR + coloc       ─► replication
 → FinnGen R12,       PP.H4 ≥ 0.8)       (INTERVAL)          (OpenGWAS,
 whole phenome)                                              non-FinnGen)
        │
        ▼
5 two-population     6 immune-cell /     7 novelty &
  replication     ─► inflammation     ─► therapeutic-direction
  (Finland→England)  (counts, CRP,       engine (tiering,
                      cytokines)          drug direction)
```

Every layer runs across the whole FinnGen R12 phenome — **none is restricted to a curated disease list.**

## Figures

- `figures/main/` (= `figures/nature/`) — 21 manuscript main figures (Figure 1–22): the calibrated autoimmune arc (1–9), the pan-phenome layer (11–16), the layer-coverage / protein-level / novelty panels (17–22).
- `figures/curation/` — 16 protein-universe curation and annotation panels (Fig01–Fig16).
- `figures/paper_style/` — 32 deep-layer panels, including `DL2b_uk_replication_expanded.png` (the widened Finland→England arm).
- `figures/heart/` — 4 cardiovascular-validation panels.
- `figures/intelligence_layer/` — 13 disease-intelligence panels.
- `figures/phenome/` — 34 per-disease / per-category phenome figures.
- `figures/supplementary/` — 70 Extended Data / Supplementary figures (per-disease MR volcano, QQ + genomic inflation λ, forest plots, per-gene INTERVAL cis-pQTL regional plots, discovery-vs-replication).

**211 figures total.**

## Manuscript & methodology

- `manuscript/Plasma_Immunome_Phenome_Atlas_Nature.docx` — the manuscript with figures embedded.
- `manuscript/PIRS_and_Atlas_Methodology.docx` — the complete step-by-step methodology, every parameter and every decision, with the figure that each step produces.
- `manuscript/Heart_Cardiovascular_Validation_Methodology.docx` — the cardiovascular deep-dive validation.
- `manuscript/Plasma_Immune_Discovery_Report_ALL.md` — the whole-phenome discovery report.

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
python src/18_opengwas_replication.py       # needs a free OpenGWAS JWT (see below)
python src/19_integrate_replication.py
python src/20_supplementary_figures.py      # 70 supplementary figures
python src/21_novelty_analysis.py
# --- whole-phenome expansion (all 2,466 FinnGen R12 endpoints) ---
python src/23_download_finngen_phenome.py
python src/41_panphenome_mr.py              # immune cis-MR across the FULL phenome
python src/43_panphenome_coloc.py           # coloc on every phenome-wide hit
python src/48_expand_pqtl_all_phenome.py    # INTERVAL pQTL over the full phenome
python src/49_pqtl_all_phenome_mr_coloc.py
python src/44_uk_panphenome_replication.py  # two-population, exact endpoint codes
python src/52_uk_replication_expanded.py    # two-population, + name-matched endpoints
python src/45_extended_cell_crp_mr.py       # immune-cell / CRP / cytokine layer
python src/50_novelty_engine_all_phenome.py # integrated novelty-priority engine
python src/42_panphenome_figures.py
python src/46_deep_layer_figures.py
python src/51_all_phenome_protein_figures.py
python src/31_disease_intelligence_layer.py # Final Required Output Table + panels
python src/13_write_manuscript.py           # assembles the .docx
python src/30_write_methodology.py          # assembles the methodology .docx
python src/47_heart_methodology.py          # cardiovascular validation .docx
```

Scripts expect the project root layout (`01_data_raw/`, `02_data_processed/`, `06_genetic_causality/`, `08_figures/`). The pre-computed **result tables are included** under `results/` so figures and the manuscript can be regenerated without re-downloading multi-GB sumstats.

> **Note:** the full phenome-wide MR output `cis_MR_ALL_finngen_results.tsv` (1,656,872 rows, ~615 MB) exceeds GitHub's file-size limit and is **not** committed. Every downstream table derived from it *is* included, and `src/41_panphenome_mr.py` regenerates it.

## Train your own model (PIRS) — bring your own data

The atlas is a **summary-statistics causal resource**, but you can also train a supervised
**Plasma Immune Risk Score (PIRS)** on **any cohort you are authorised to use** (UK Biobank,
a clinical cohort, your own Olink/other proteomic run). `src/22_train_pirs.py` ships **no
individual-level data** and never fabricates any.

```bash
python src/22_train_pirs.py --write-templates           # blank input templates
python src/22_train_pirs.py --npx npx.tsv --outcomes outcomes.tsv
```

You provide three simple tables — `npx` (participant × protein, columns named by Olink ID or
gene symbol), `outcomes` (`id, disease, event, time_years`), and optional `covariates`. The
trainer restricts features to the 1,007 curated immune proteins, runs **cross-validated
elastic-net survival** (scikit-survival Coxnet → lifelines → logistic fallback, auto-selected),
and writes per-protein PIRS weights, a CV C-index/AUROC table, model pickles, and a
performance figure.

**➡ Full step-by-step guide: [`docs/TRAINING.md`](docs/TRAINING.md)** — input schema, options,
how to score new individuals, and validation notes.

> UK Biobank Olink NPX + phenotypes are controlled-access (apply via the UKB Access Management
> System, Olink Field 30900). PIRS is not shipped pre-trained because doing so would require
> individual-level data.

## Data sources (all public)

| Layer | Source | Access |
|---|---|---|
| Protein universe | Olink Explore (UKB coding 143) | public |
| Annotation | Human Protein Atlas, MSigDB C7/C8 | public |
| eQTL instruments | eQTLGen cis-eQTLs (n ≈ 31,684) | public |
| Disease GWAS (discovery) | FinnGen R12, all 2,466 endpoints | public |
| Protein pQTL | INTERVAL (Sun 2018) via EBI GWAS Catalog FTP | public |
| Replication GWAS | non-FinnGen consortium GWAS via **OpenGWAS** | free JWT token |
| Two-population replication | UK Biobank single-cohort GWAS (`ukb-d/b/a/e-*`) | free JWT token |

> **OpenGWAS token:** replication (`src/18`, `src/44`, `src/52`) needs a free personal JWT from https://api.opengwas.io — register, place the token in a local file, and never commit it. The included replication tables already contain the saved output, so re-running those steps is optional.

## Limits

* MR estimates are **lifelong genetically-proxied** effects, not the effect of a drug course.
* MHC/HLA-region hits are held at *nomination* — long-range LD defeats colocalisation there.
* eQTL instruments proxy transcript, not always circulating protein; layer 3 is the arbiter.
* 110 of the 220 causal proteins have no public SomaScan aptamer, so they cannot reach a protein-level tier from login-free data.
* The two-population layer covers 62 diseases, not the full phenome, because it requires an independent single-cohort UK Biobank GWAS of the same endpoint to exist at all.
* Research resource only — **not** clinical advice or a validated diagnostic.

## Not included (gated, by design)

Individual-level UK Biobank Olink NPX + phenotypes (controlled access) and UKB-PPP / deCODE pQTL (Synapse / DUA-gated) are **not** in this repository and are **not required** — the atlas is built to be fully reproducible from open data.

## License

- **Code** (`src/`): MIT — see `LICENSE`.
- **Data tables & figures** (`results/`, `figures/`): CC-BY-4.0.

## Citation

See `CITATION.cff`.
