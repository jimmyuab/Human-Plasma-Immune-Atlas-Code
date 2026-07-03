#!/usr/bin/env python
"""
HDDM Layer 54 - Step 22
Plasma Immune Risk Score (PIRS) trainer  --  READY-TO-RUN SCAFFOLD.

This module trains a cross-validated, elastic-net risk model that maps an
individual's plasma immune Olink NPX profile to future disease risk, restricted
to the 1,007 curated plasma-immune proteins.

IMPORTANT — HONEST-BY-DESIGN:
  Individual-level UK Biobank Olink NPX + phenotypes are CONTROLLED-ACCESS and
  are NOT shipped with this repository. This script therefore does NOT fabricate
  data. If the expected input files are absent it prints the exact schema it
  needs and exits cleanly (code 0). The moment approved UKB data is dropped into
  the paths below, `python src/22_train_pirs.py` trains PIRS end-to-end.

EXPECTED INPUTS (place under 01_data_raw/UKB_Olink_NPX/):
  1. npx_matrix.tsv        rows = participants, cols = ['eid', <olink_id ...>]
                           values = Olink NPX (log2, one column per assay).
  2. outcomes.tsv          rows = participants, cols =
                           ['eid', 'disease', 'event', 'time_years']
                           event = 1 incident / 0 censored ; time_years = follow-up
                           (baseline -> first event or censoring). One row per
                           participant x modelled disease.
  (optional) covariates.tsv  ['eid','age','sex', ... ]  adjusted-for confounders.

OUTPUTS (05_machine_learning/):
  pirs_<disease>_weights.tsv     per-protein PIRS coefficients (the score).
  pirs_cv_metrics.tsv            per-disease cross-validated discrimination.
  pirs_<disease>_model.pkl       fitted pipeline (impute+scale+model).
  Figure10_pirs_performance.png  CV performance across modelled diseases.
"""
import os, sys, glob, pickle
import numpy as np, pandas as pd

ROOT   = r"I:\Plasma immune atalas"
PROC   = os.path.join(ROOT, "02_data_processed")
UKB    = os.path.join(ROOT, "01_data_raw", "UKB_Olink_NPX")
OUT    = os.path.join(ROOT, "05_machine_learning")
FIG    = os.path.join(ROOT, "08_figures", "nature")
os.makedirs(OUT, exist_ok=True)

NPX_FP  = os.path.join(UKB, "npx_matrix.tsv")
OUT_FP  = os.path.join(UKB, "outcomes.tsv")
COV_FP  = os.path.join(UKB, "covariates.tsv")

# --- config ---
N_FOLDS       = 5
SEED          = 42
L1_RATIO      = 0.5          # elastic-net mix (0=ridge, 1=lasso)
HORIZON_YEARS = 10          # used only by the logistic fallback
MIN_EVENTS    = 50          # skip a disease with too few incident cases

# ============================================================
# 0. GUARD: no fabricated data. If inputs are missing, instruct and exit.
# ============================================================
def _instruct_and_exit():
    print("="*72)
    print("PIRS trainer: UK Biobank individual-level inputs NOT found.")
    print("This is expected until controlled-access UKB data is approved and")
    print("placed locally. No data is fabricated. Provide these files:\n")
    print(f"  {NPX_FP}")
    print("     rows=participants  cols=['eid', <olink_id...>]  values=NPX (log2)")
    print(f"  {OUT_FP}")
    print("     cols=['eid','disease','event','time_years']  (one row / eid / disease)")
    print(f"  {COV_FP}  (optional) cols=['eid','age','sex',...]")
    print("\nApply via the UK Biobank Access Management System (Olink Field 30900")
    print("+ disease phenotypes). Then re-run:  python src/22_train_pirs.py")
    print("="*72)
    sys.exit(0)

if not (os.path.exists(NPX_FP) and os.path.exists(OUT_FP)):
    _instruct_and_exit()

# ============================================================
# 1. Load + restrict features to the curated plasma immunome
# ============================================================
ann = pd.read_csv(os.path.join(PROC, "plasma_immune_protein_annotation.tsv"), sep="\t")
imm = ann[ann.is_plasma_immune == 1]
immune_assays = imm["olink_id"].dropna().astype(str).tolist()

npx = pd.read_csv(NPX_FP, sep="\t")
out = pd.read_csv(OUT_FP, sep="\t")
cov = pd.read_csv(COV_FP, sep="\t") if os.path.exists(COV_FP) else None

feat_cols = [c for c in npx.columns if c in set(immune_assays)]
if not feat_cols:
    print("ERROR: no immune-protein columns matched annotation olink_id. "
          "Check that npx_matrix.tsv columns use Olink assay IDs.")
    sys.exit(1)
print(f"features: {len(feat_cols)} immune-protein assays "
      f"(of {len(immune_assays)} curated); participants: {npx.eid.nunique()}")

# ============================================================
# 2. Model backend selection (proper survival if available)
# ============================================================
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold

BACKEND = None
try:
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.metrics import concordance_index_censored
    from sksurv.util import Surv
    BACKEND = "coxnet"
except Exception:
    try:
        from lifelines import CoxPHFitter
        from lifelines.utils import concordance_index
        BACKEND = "lifelines"
    except Exception:
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import roc_auc_score
        BACKEND = "logistic"
print(f"survival backend: {BACKEND}")

def make_pre():
    return Pipeline([("impute", SimpleImputer(strategy="median")),
                     ("scale",  StandardScaler())])

# ============================================================
# 3. Per-disease cross-validated training
# ============================================================
metrics_rows = []
diseases = sorted(out.disease.unique())

for dis in diseases:
    od = out[out.disease == dis][["eid", "event", "time_years"]]
    df = npx[["eid"] + feat_cols].merge(od, on="eid", how="inner").dropna(subset=["event", "time_years"])
    if cov is not None:
        df = df.merge(cov, on="eid", how="left")
    n_ev = int(df.event.sum())
    if n_ev < MIN_EVENTS:
        print(f"skip {dis}: only {n_ev} incident events (<{MIN_EVENTS})")
        continue

    X = df[feat_cols].values
    y_event = df.event.astype(int).values
    y_time  = df.time_years.astype(float).values
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)

    fold_c = []
    for tr, te in skf.split(X, y_event):
        pre = make_pre().fit(X[tr])
        Xtr, Xte = pre.transform(X[tr]), pre.transform(X[te])
        if BACKEND == "coxnet":
            ytr = Surv.from_arrays(event=y_event[tr].astype(bool), time=y_time[tr])
            mdl = CoxnetSurvivalAnalysis(l1_ratio=L1_RATIO, alpha_min_ratio=0.01,
                                         max_iter=100000)
            mdl.fit(Xtr, ytr)
            risk = mdl.predict(Xte)
            c = concordance_index_censored(y_event[te].astype(bool), y_time[te], risk)[0]
        elif BACKEND == "lifelines":
            tr_df = pd.DataFrame(Xtr, columns=feat_cols)
            tr_df["T"], tr_df["E"] = y_time[tr], y_event[tr]
            cph = CoxPHFitter(penalizer=0.1, l1_ratio=L1_RATIO)
            cph.fit(tr_df, "T", "E")
            risk = cph.predict_partial_hazard(pd.DataFrame(Xte, columns=feat_cols)).values.ravel()
            c = concordance_index(y_time[te], -risk, y_event[te])
        else:  # logistic fallback: incident within HORIZON
            lab_tr = ((y_event[tr] == 1) & (y_time[tr] <= HORIZON_YEARS)).astype(int)
            lab_te = ((y_event[te] == 1) & (y_time[te] <= HORIZON_YEARS)).astype(int)
            if lab_tr.sum() < 5 or len(set(lab_te)) < 2:
                continue
            mdl = LogisticRegression(penalty="elasticnet", solver="saga",
                                     l1_ratio=L1_RATIO, C=1.0, max_iter=5000)
            mdl.fit(Xtr, lab_tr)
            risk = mdl.predict_proba(Xte)[:, 1]
            c = roc_auc_score(lab_te, risk)
        fold_c.append(c)

    if not fold_c:
        continue
    cmean, cstd = float(np.mean(fold_c)), float(np.std(fold_c))
    metric = "C-index" if BACKEND != "logistic" else "AUROC"
    print(f"{dis:32s} {metric}={cmean:.3f} +/- {cstd:.3f}  (n_ev={n_ev})")
    metrics_rows.append(dict(disease=dis, backend=BACKEND, metric=metric,
                             score_mean=cmean, score_std=cstd, n_events=n_ev,
                             n_participants=len(df)))

    # --- refit on ALL data -> save PIRS weights + model ---
    pre = make_pre().fit(X)
    Xall = pre.transform(X)
    if BACKEND == "coxnet":
        yall = Surv.from_arrays(event=y_event.astype(bool), time=y_time)
        final = CoxnetSurvivalAnalysis(l1_ratio=L1_RATIO, alpha_min_ratio=0.01, max_iter=100000).fit(Xall, yall)
        coef = final.coef_[:, -1]
    elif BACKEND == "lifelines":
        alldf = pd.DataFrame(Xall, columns=feat_cols); alldf["T"], alldf["E"] = y_time, y_event
        final = CoxPHFitter(penalizer=0.1, l1_ratio=L1_RATIO).fit(alldf, "T", "E")
        coef = final.params_.values
    else:
        lab = ((y_event == 1) & (y_time <= HORIZON_YEARS)).astype(int)
        final = LogisticRegression(penalty="elasticnet", solver="saga",
                                   l1_ratio=L1_RATIO, C=1.0, max_iter=5000).fit(Xall, lab)
        coef = final.coef_.ravel()

    id2gene = dict(zip(imm.olink_id.astype(str), imm.gene_symbol))
    w = pd.DataFrame({"olink_id": feat_cols,
                      "gene_symbol": [id2gene.get(a, a) for a in feat_cols],
                      "pirs_weight": coef}).sort_values("pirs_weight", key=np.abs, ascending=False)
    safe = dis.replace(" ", "_").replace("/", "_")[:24]
    w.to_csv(os.path.join(OUT, f"pirs_{safe}_weights.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT, f"pirs_{safe}_model.pkl"), "wb") as fh:
        pickle.dump({"pre": pre, "model": final, "features": feat_cols, "backend": BACKEND}, fh)

# ============================================================
# 4. Metrics table + performance figure
# ============================================================
if not metrics_rows:
    print("No disease had enough incident events to train PIRS.")
    sys.exit(0)

mt = pd.DataFrame(metrics_rows).sort_values("score_mean", ascending=False)
mt.to_csv(os.path.join(OUT, "pirs_cv_metrics.tsv"), sep="\t", index=False)
print("\nwrote pirs_cv_metrics.tsv")

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, max(3, 0.4*len(mt)+1)))
y = np.arange(len(mt))
ax.barh(y, mt.score_mean, xerr=mt.score_std, color="#1f4e79", ecolor="#888", capsize=3)
ax.axvline(0.5, ls="--", c="red", lw=1, label="no discrimination (0.5)")
ax.set_yticks(y); ax.set_yticklabels(mt.disease, fontsize=8)
ax.set_xlabel(f"cross-validated {mt.metric.iloc[0]} (mean +/- SD, {N_FOLDS}-fold)")
ax.set_xlim(0.4, 1.0); ax.legend(fontsize=8, loc="lower right")
for yi, (_, r) in enumerate(mt.iterrows()):
    ax.text(r.score_mean + 0.005, yi, f"{r.score_mean:.3f}", va="center", fontsize=7)
fig.suptitle("Figure 10  |  Plasma Immune Risk Score (PIRS) discrimination across diseases",
             fontsize=12.5, fontweight="bold", x=0.02, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(FIG, "Figure10_pirs_performance.png"), dpi=150, bbox_inches="tight")
print("wrote Figure10_pirs_performance.png")
