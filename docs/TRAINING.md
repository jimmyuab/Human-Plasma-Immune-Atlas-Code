# Training a Plasma Immune Risk Score (PIRS) on your own data

`src/22_train_pirs.py` lets you train a **cross-validated, elastic-net survival model**
that maps an individual's plasma immune protein profile to future disease risk — using
**your own cohort**. Features are automatically restricted to the 1,007 curated
plasma-immune proteins in this atlas.

The trainer ships **no individual-level data** and never fabricates any. You bring a cohort
you are authorised to use (UK Biobank, a hospital/clinical cohort, your own Olink or other
proteomic run).

---

## 1. Install

```bash
pip install -r requirements.txt
# optional, for a proper Cox survival model (recommended):
pip install scikit-survival        # OR:  pip install lifelines
```

Without a survival library the trainer automatically falls back to a logistic
"incident-within-horizon" model — it still runs, just with AUROC instead of C-index.

## 2. Prepare three tables

All files are TSV (or CSV). The participant-id column is `id` by default
(change with `--id-col`).

**`npx.tsv`** — one row per participant, one column per protein.
Protein columns may be named by **Olink assay ID or by gene symbol**; both are matched
against the immune panel. Values are normalised expression (e.g. Olink NPX, log2).

| id | IL6 | CTLA4 | TNFSF14 | ... |
|----|-----|-------|---------|-----|
| 1001 | 2.13 | 0.88 | 1.42 | ... |
| 1002 | 1.05 | 1.30 | 0.77 | ... |

**`outcomes.tsv`** — one row per participant **per disease** you want to model.
`event` = 1 if the disease occurred during follow-up else 0; `time_years` = time from
baseline to the event or to censoring.

| id | disease | event | time_years |
|----|---------|-------|-----------|
| 1001 | Rheumatoid arthritis | 1 | 4.2 |
| 1001 | Multiple sclerosis   | 0 | 11.0 |
| 1002 | Rheumatoid arthritis | 0 | 11.0 |

**`covariates.tsv`** *(optional)* — confounders to include.

| id | age | sex |
|----|-----|-----|
| 1001 | 58 | 1 |

> Tip: generate blank templates with the exact headers:
> ```bash
> python src/22_train_pirs.py --write-templates
> ```

## 3. Train

```bash
python src/22_train_pirs.py --npx npx.tsv --outcomes outcomes.tsv
# with covariates and a custom id column:
python src/22_train_pirs.py --npx npx.tsv --outcomes outcomes.tsv \
       --covariates covariates.tsv --id-col participant_id
```

Useful options:

| flag | default | meaning |
|------|---------|---------|
| `--folds` | 5 | cross-validation folds |
| `--l1-ratio` | 0.5 | elastic-net mix (0 = ridge, 1 = lasso) |
| `--min-events` | 50 | skip a disease with fewer incident cases |
| `--horizon` | 10 | years, used only by the logistic fallback |
| `--outdir` | `05_machine_learning/` | where results are written |
| `--annotation` | auto | override the immune-panel annotation path |

## 4. Outputs (`--outdir`)

| file | contents |
|------|----------|
| `pirs_<disease>_weights.tsv` | per-protein PIRS coefficients — **the score** |
| `pirs_cv_metrics.tsv` | cross-validated C-index (or AUROC) per disease |
| `pirs_<disease>_model.pkl` | fitted pipeline (median-impute → scale → model) |
| `pirs_performance.png` | discrimination across all modelled diseases |

## 5. Score new people with a trained model

```python
import pandas as pd, pickle
m = pickle.load(open("05_machine_learning/pirs_Rheumatoid_arthritis_model.pkl", "rb"))
new = pd.read_csv("new_npx.tsv", sep="\t")          # same protein columns
X = m["pre"].transform(new[m["features"]].values)   # impute + scale
risk = (m["model"].predict(X) if m["backend"] == "coxnet"
        else m["model"].predict_partial_hazard(pd.DataFrame(X, columns=m["features"])).values
        if m["backend"] == "lifelines"
        else m["model"].predict_proba(X)[:, 1])       # higher = higher risk
```

## Notes on rigour

- **Discrimination is cross-validated** (stratified k-fold); the reported C-index/AUROC is
  out-of-fold. The saved weights are refit on the full cohort.
- **Validate externally.** A model trained on one cohort should be tested on an independent
  cohort before any claim of clinical utility.
- **Immune-only by design.** Features are the atlas's 1,007 plasma-immune proteins, so PIRS
  is an *immune* risk score; add non-immune proteins or clinical covariates only if that
  matches your question.
- **Ethics/consent.** Use only data you are authorised and consented to analyse. UK Biobank
  Olink NPX + phenotypes require an approved application (Olink Field 30900).
