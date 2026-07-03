#!/usr/bin/env python
"""
HDDM Layer 54 - Step 22
Plasma Immune Risk Score (PIRS) trainer.

Trains a cross-validated, elastic-net survival/risk model that maps an
individual's plasma immune protein profile to future disease risk, restricted
to the 1,007 curated plasma-immune proteins.

BRING YOUR OWN DATA
-------------------
This trainer is data-source-agnostic: point it at ANY cohort you are authorised
to use (UK Biobank, a hospital cohort, your own Olink/other proteomic run).
It does NOT ship or fabricate individual-level data. If the input files are
absent it prints the exact schema it needs (and can write blank templates) and
exits cleanly.

INPUTS (tab- or comma-separated; paths configurable via CLI):
  --npx        proteomic matrix : rows = participants
                 cols = ['id', <protein columns...>]  values = normalised
                 expression (e.g. Olink NPX log2, or any per-protein level).
                 Protein columns may be named by Olink assay ID OR gene symbol;
                 both are matched against the curated immune panel.
  --outcomes   survival labels  : cols = ['id','disease','event','time_years']
                 event = 1 incident / 0 censored ; time_years = follow-up from
                 baseline to first event or censoring. One row per id x disease.
  --covariates (optional)       : cols = ['id','age','sex', ...] confounders.
  --id-col     name of the participant-id column (default: 'id').

OUTPUTS (--outdir, default 05_machine_learning/):
  pirs_<disease>_weights.tsv     per-protein PIRS coefficients (the score).
  pirs_cv_metrics.tsv            per-disease cross-validated discrimination.
  pirs_<disease>_model.pkl       fitted pipeline (impute+scale+model).
  pirs_performance.png           CV performance across modelled diseases.

QUICK START
-----------
  # 1. see the schema / write blank templates you can fill in
  python src/22_train_pirs.py --write-templates
  # 2. train on your data
  python src/22_train_pirs.py --npx my_npx.tsv --outcomes my_outcomes.tsv
"""
import os, sys, argparse, pickle
import numpy as np, pandas as pd

# ---- portable project root (works wherever the repo is cloned) ----
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

def find_annotation():
    """Locate the curated immune-protein annotation across repo layouts."""
    for c in [os.path.join(ROOT, "02_data_processed", "plasma_immune_protein_annotation.tsv"),
              os.path.join(ROOT, "results", "annotation", "plasma_immune_protein_annotation.tsv"),
              os.path.join(ROOT, "plasma_immune_protein_annotation.tsv")]:
        if os.path.exists(c):
            return c
    return None

def read_table(fp):
    sep = "\t" if fp.lower().endswith((".tsv", ".txt")) else ","
    return pd.read_csv(fp, sep=sep)

# ============================================================
# CLI
# ============================================================
ap = argparse.ArgumentParser(description="Train a Plasma Immune Risk Score (PIRS) on your cohort.")
ap.add_argument("--npx",        default=None, help="proteomic matrix (rows=participants, cols=id+proteins)")
ap.add_argument("--outcomes",   default=None, help="survival labels (id,disease,event,time_years)")
ap.add_argument("--covariates", default=None, help="optional confounders (id,age,sex,...)")
ap.add_argument("--id-col",     default="id",  help="participant id column name (default: id)")
ap.add_argument("--annotation", default=None, help="override path to immune-protein annotation")
ap.add_argument("--outdir",     default=os.path.join(ROOT, "05_machine_learning"))
ap.add_argument("--figdir",     default=os.path.join(ROOT, "08_figures", "nature"))
ap.add_argument("--folds",      type=int,   default=5)
ap.add_argument("--l1-ratio",   type=float, default=0.5, help="elastic-net mix (0=ridge,1=lasso)")
ap.add_argument("--horizon",    type=float, default=10,  help="years, logistic fallback only")
ap.add_argument("--min-events", type=int,   default=50,  help="skip disease with fewer incident cases")
ap.add_argument("--seed",       type=int,   default=42)
ap.add_argument("--write-templates", action="store_true",
                help="write blank input templates next to --outdir and exit")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# ============================================================
# 0. templates / guard — never fabricate data
# ============================================================
SCHEMA = f"""
INPUT SCHEMA (id column = '{args.id_col}')
  npx        : {args.id_col}, <protein_1>, <protein_2>, ...   (values = NPX/expression)
  outcomes   : {args.id_col}, disease, event, time_years      (event 1/0; time in years)
  covariates : {args.id_col}, age, sex, ...                   (optional)
Protein columns may be Olink assay IDs or gene symbols; both are matched to the
curated immune panel."""

def write_templates():
    d = args.outdir
    idc = args.id_col
    pd.DataFrame(columns=[idc, "PROTEIN_A", "PROTEIN_B"]).to_csv(
        os.path.join(d, "TEMPLATE_npx.tsv"), sep="\t", index=False)
    pd.DataFrame(columns=[idc, "disease", "event", "time_years"]).to_csv(
        os.path.join(d, "TEMPLATE_outcomes.tsv"), sep="\t", index=False)
    pd.DataFrame(columns=[idc, "age", "sex"]).to_csv(
        os.path.join(d, "TEMPLATE_covariates.tsv"), sep="\t", index=False)
    print("Wrote blank templates to", d)
    print(SCHEMA)

if args.write_templates:
    write_templates(); sys.exit(0)

if not args.npx or not args.outcomes or not os.path.exists(args.npx) or not os.path.exists(args.outcomes):
    print("="*72)
    print("PIRS trainer: no input data provided.")
    print("This trainer needs an individual-level cohort you are authorised to use.")
    print("It does not ship or fabricate data.\n")
    print("  python src/22_train_pirs.py --write-templates       # blank input files")
    print("  python src/22_train_pirs.py --npx NPX.tsv --outcomes OUT.tsv")
    print(SCHEMA)
    print("="*72)
    sys.exit(0)

# ============================================================
# 1. Load + restrict features to the curated plasma immunome
# ============================================================
ann_fp = args.annotation or find_annotation()
if ann_fp is None:
    print("ERROR: could not find plasma_immune_protein_annotation.tsv. "
          "Pass it with --annotation."); sys.exit(1)
ann = read_table(ann_fp)
imm = ann[ann.is_plasma_immune == 1]
immune_ids   = set(imm["olink_id"].dropna().astype(str))
immune_genes = set(imm["gene_symbol"].dropna().astype(str))
id2gene = dict(zip(imm.olink_id.astype(str), imm.gene_symbol))

npx = read_table(args.npx)
out = read_table(args.outcomes)
cov = read_table(args.covariates) if (args.covariates and os.path.exists(args.covariates)) else None
idc = args.id_col
for nm, df in [("npx", npx), ("outcomes", out)]:
    if idc not in df.columns:
        print(f"ERROR: id column '{idc}' not found in {nm} (use --id-col)."); sys.exit(1)

# match protein columns by Olink id OR gene symbol
feat_cols = [c for c in npx.columns if c != idc and (str(c) in immune_ids or str(c) in immune_genes)]
if not feat_cols:
    print("ERROR: none of the proteomic columns matched the curated immune panel.\n"
          "Columns should be Olink assay IDs or gene symbols. Example immune genes: "
          + ", ".join(list(immune_genes)[:8])); sys.exit(1)
print(f"features: {len(feat_cols)} immune proteins matched "
      f"(of {len(immune_ids)} curated); participants: {npx[idc].nunique()}")

# ============================================================
# 2. Survival backend (proper Cox if available, else logistic)
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
for dis in sorted(out.disease.unique()):
    od = out[out.disease == dis][[idc, "event", "time_years"]]
    df = npx[[idc] + feat_cols].merge(od, on=idc, how="inner").dropna(subset=["event", "time_years"])
    if cov is not None:
        df = df.merge(cov, on=idc, how="left")
    n_ev = int(df.event.sum())
    if n_ev < args.min_events:
        print(f"skip {dis}: only {n_ev} incident events (<{args.min_events})")
        continue

    X = df[feat_cols].values
    y_event = df.event.astype(int).values
    y_time  = df.time_years.astype(float).values
    skf = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=args.seed)

    fold_c = []
    for tr, te in skf.split(X, y_event):
        pre = make_pre().fit(X[tr]); Xtr, Xte = pre.transform(X[tr]), pre.transform(X[te])
        if BACKEND == "coxnet":
            ytr = Surv.from_arrays(event=y_event[tr].astype(bool), time=y_time[tr])
            mdl = CoxnetSurvivalAnalysis(l1_ratio=args.l1_ratio, alpha_min_ratio=0.01, max_iter=100000)
            mdl.fit(Xtr, ytr); risk = mdl.predict(Xte)
            c = concordance_index_censored(y_event[te].astype(bool), y_time[te], risk)[0]
        elif BACKEND == "lifelines":
            tr_df = pd.DataFrame(Xtr, columns=feat_cols); tr_df["T"], tr_df["E"] = y_time[tr], y_event[tr]
            cph = CoxPHFitter(penalizer=0.1, l1_ratio=args.l1_ratio); cph.fit(tr_df, "T", "E")
            risk = cph.predict_partial_hazard(pd.DataFrame(Xte, columns=feat_cols)).values.ravel()
            c = concordance_index(y_time[te], -risk, y_event[te])
        else:
            lab_tr = ((y_event[tr] == 1) & (y_time[tr] <= args.horizon)).astype(int)
            lab_te = ((y_event[te] == 1) & (y_time[te] <= args.horizon)).astype(int)
            if lab_tr.sum() < 5 or len(set(lab_te)) < 2: continue
            mdl = LogisticRegression(penalty="elasticnet", solver="saga",
                                     l1_ratio=args.l1_ratio, C=1.0, max_iter=5000)
            mdl.fit(Xtr, lab_tr); risk = mdl.predict_proba(Xte)[:, 1]
            c = roc_auc_score(lab_te, risk)
        fold_c.append(c)
    if not fold_c: continue

    cmean, cstd = float(np.mean(fold_c)), float(np.std(fold_c))
    metric = "C-index" if BACKEND != "logistic" else "AUROC"
    print(f"{dis:32s} {metric}={cmean:.3f} +/- {cstd:.3f}  (n_ev={n_ev})")
    metrics_rows.append(dict(disease=dis, backend=BACKEND, metric=metric,
                             score_mean=cmean, score_std=cstd, n_events=n_ev, n_participants=len(df)))

    # refit on all data -> PIRS weights + saved model
    pre = make_pre().fit(X); Xall = pre.transform(X)
    if BACKEND == "coxnet":
        yall = Surv.from_arrays(event=y_event.astype(bool), time=y_time)
        final = CoxnetSurvivalAnalysis(l1_ratio=args.l1_ratio, alpha_min_ratio=0.01, max_iter=100000).fit(Xall, yall)
        coef = final.coef_[:, -1]
    elif BACKEND == "lifelines":
        alldf = pd.DataFrame(Xall, columns=feat_cols); alldf["T"], alldf["E"] = y_time, y_event
        final = CoxPHFitter(penalizer=0.1, l1_ratio=args.l1_ratio).fit(alldf, "T", "E"); coef = final.params_.values
    else:
        lab = ((y_event == 1) & (y_time <= args.horizon)).astype(int)
        final = LogisticRegression(penalty="elasticnet", solver="saga",
                                   l1_ratio=args.l1_ratio, C=1.0, max_iter=5000).fit(Xall, lab); coef = final.coef_.ravel()

    w = pd.DataFrame({"protein": feat_cols,
                      "gene_symbol": [id2gene.get(str(a), a) for a in feat_cols],
                      "pirs_weight": coef}).sort_values("pirs_weight", key=np.abs, ascending=False)
    safe = dis.replace(" ", "_").replace("/", "_")[:24]
    w.to_csv(os.path.join(args.outdir, f"pirs_{safe}_weights.tsv"), sep="\t", index=False)
    with open(os.path.join(args.outdir, f"pirs_{safe}_model.pkl"), "wb") as fh:
        pickle.dump({"pre": pre, "model": final, "features": feat_cols, "backend": BACKEND}, fh)

# ============================================================
# 4. Metrics table + performance figure
# ============================================================
if not metrics_rows:
    print("No disease had enough incident events to train PIRS."); sys.exit(0)

mt = pd.DataFrame(metrics_rows).sort_values("score_mean", ascending=False)
mt.to_csv(os.path.join(args.outdir, "pirs_cv_metrics.tsv"), sep="\t", index=False)
print("\nwrote pirs_cv_metrics.tsv")

os.makedirs(args.figdir, exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, max(3, 0.4*len(mt)+1)))
y = np.arange(len(mt))
ax.barh(y, mt.score_mean, xerr=mt.score_std, color="#1f4e79", ecolor="#888", capsize=3)
ax.axvline(0.5, ls="--", c="red", lw=1, label="no discrimination (0.5)")
ax.set_yticks(y); ax.set_yticklabels(mt.disease, fontsize=8)
ax.set_xlabel(f"cross-validated {mt.metric.iloc[0]} (mean +/- SD, {args.folds}-fold)")
ax.set_xlim(0.4, 1.0); ax.legend(fontsize=8, loc="lower right")
for yi, (_, r) in enumerate(mt.iterrows()):
    ax.text(r.score_mean + 0.005, yi, f"{r.score_mean:.3f}", va="center", fontsize=7)
fig.suptitle("Plasma Immune Risk Score (PIRS) discrimination across diseases",
             fontsize=12.5, fontweight="bold", x=0.02, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(args.figdir, "pirs_performance.png"), dpi=150, bbox_inches="tight")
print("wrote pirs_performance.png")
