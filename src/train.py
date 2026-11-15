"""Train and evaluate per-drug resistance models with leakage-aware validation.

For each drug and model, out-of-fold predictions come from two schemes:
  grouped : StratifiedGroupKFold on Mash lineage groups (the honest estimate;
            genomes from one lineage never sit in both train and test)
  random  : StratifiedKFold ignoring lineage (shows how much leakage inflates scores)
Feature filtering (dropping rare/constant elements) is fitted inside each fold.
Metrics come with 95 % confidence intervals from a bootstrap that resamples
whole lineage groups. A leave-one-BioProject-out test checks transfer to an
unseen study.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

DRUGS = ["meropenem", "ciprofloxacin", "gentamicin", "trimethoprim/sulfamethoxazole"]
SEED = 42
N_SPLITS = 5


class MinPrevalence(BaseEstimator, TransformerMixin):
    """Keep features present in at least min_frac and at most 1 - min_frac of the
    TRAINING rows. Fitted per fold, so test rows never influence the selection."""

    def __init__(self, min_frac: float = 0.01):
        self.min_frac = min_frac

    def fit(self, X, y=None):
        prev = np.asarray(X).mean(axis=0)
        self.keep_ = (prev >= self.min_frac) & (prev <= 1 - self.min_frac)
        if not self.keep_.any():
            self.keep_ = prev > 0
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.keep_]


def make_model(name: str, y_train: np.ndarray):
    """Fresh, untrained model. Class imbalance is handled by class weights."""
    if name == "logreg":
        return LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
    if name == "rf":
        return RandomForestClassifier(n_estimators=500, min_samples_leaf=2,
                                      class_weight="balanced_subsample",
                                      n_jobs=-1, random_state=SEED)
    if name == "xgb":
        from xgboost import XGBClassifier
        pos = max(int(y_train.sum()), 1)
        return XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.05,
                             subsample=0.8, colsample_bytree=0.8,
                             scale_pos_weight=(len(y_train) - pos) / pos,
                             eval_metric="logloss", n_jobs=4, random_state=SEED)
    raise ValueError(name)


def fit_predict(model_name, X_tr, y_tr, X_te):
    filt = MinPrevalence().fit(X_tr)
    model = make_model(model_name, y_tr).fit(filt.transform(X_tr), y_tr)
    return model.predict_proba(filt.transform(X_te))[:, 1]


def oof_predictions(model_name, X, y, groups, scheme):
    """Out-of-fold probabilities for every genome."""
    if scheme == "grouped":
        splits = StratifiedGroupKFold(N_SPLITS, shuffle=True, random_state=SEED).split(X, y, groups)
    else:
        splits = StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED).split(X, y)
    proba = np.full(len(y), np.nan)
    for tr, te in splits:
        proba[te] = fit_predict(model_name, X[tr], y[tr], X[te])
    return proba


def metrics(y, proba, threshold=0.5) -> dict:
    pred = (proba >= threshold).astype(int)
    out = {"auroc": np.nan, "f1": f1_score(y, pred, zero_division=0),
           "mcc": matthews_corrcoef(y, pred)}
    if len(np.unique(y)) == 2:
        out["auroc"] = roc_auc_score(y, proba)
    return out


def group_bootstrap_ci(y, proba, groups, n_boot=1000, seed=SEED) -> dict:
    """95 % CI per metric, resampling whole lineage groups with replacement."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    idx_by_group = {g: np.flatnonzero(groups == g) for g in uniq}
    draws = {k: [] for k in ("auroc", "f1", "mcc")}
    for _ in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([idx_by_group[g] for g in pick])
        if len(np.unique(y[idx])) < 2:
            continue
        for k, v in metrics(y[idx], proba[idx]).items():
            draws[k].append(v)
    return {f"{k}_lo": np.percentile(v, 2.5) for k, v in draws.items()} | \
           {f"{k}_hi": np.percentile(v, 97.5) for k, v in draws.items()}


def load_inputs(matrix, groups, labels):
    X = pd.read_csv(matrix, sep="\t", index_col="genome_id", dtype={"genome_id": str})
    X.index = X.index.astype(str)
    g = pd.read_csv(groups, sep="\t", dtype={"genome_id": str}).set_index("genome_id")["lineage"]
    lab = pd.read_csv(labels, sep="\t", dtype={"genome_id": str}).set_index("genome_id")
    common = X.index.intersection(g.index).intersection(lab.index)
    return X.loc[common], g.loc[common], lab.loc[common]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", default="data/features/amr_matrix.tsv")
    p.add_argument("--groups", default="data/features/lineage_groups.tsv")
    p.add_argument("--labels", default="data/metadata/selected_500.tsv")
    p.add_argument("--models", default="logreg,rf,xgb")
    p.add_argument("--holdout-project", default="PRJEB27256")
    p.add_argument("--out", default="reports/per_drug_results.csv")
    p.add_argument("--pred-out", default="results/oof_predictions.tsv")
    p.add_argument("--n-boot", type=int, default=1000)
    args = p.parse_args()

    X_df, groups_s, lab = load_inputs(args.matrix, args.groups, args.labels)
    X_all = X_df.values.astype(np.int8)
    print(f"genomes: {len(X_df)}   features: {X_df.shape[1]}   "
          f"lineage groups: {groups_s.nunique()}")

    rows, preds = [], []
    for drug in DRUGS:
        y = lab[drug].values.astype(int)
        g = groups_s.values
        for model_name in args.models.split(","):
            for scheme in ("grouped", "random"):
                proba = oof_predictions(model_name, X_all, y, g, scheme)
                row = {"drug": drug, "model": model_name, "validation": scheme,
                       "n": len(y), "n_resistant": int(y.sum())}
                row |= metrics(y, proba) | group_bootstrap_ci(y, proba, g, args.n_boot)
                rows.append(row)
                preds.append(pd.DataFrame({"genome_id": X_df.index, "drug": drug,
                                           "model": model_name, "validation": scheme,
                                           "label": y, "proba": proba}))
                print(f"{drug:31s} {model_name:6s} {scheme:7s} "
                      f"AUROC {row['auroc']:.3f} [{row['auroc_lo']:.3f}-{row['auroc_hi']:.3f}]  "
                      f"F1 {row['f1']:.3f}  MCC {row['mcc']:.3f}")

        # Leave-one-BioProject-out: train on every other study, test on this one.
        proj = lab["bioproject_accession"].fillna("").values
        te = proj == args.holdout_project
        if te.sum() and len(np.unique(y[te])) == 2:
            model_name = args.models.split(",")[-1]
            proba = fit_predict(model_name, X_all[~te], y[~te], X_all[te])
            row = {"drug": drug, "model": model_name,
                   "validation": f"holdout:{args.holdout_project}",
                   "n": int(te.sum()), "n_resistant": int(y[te].sum())}
            row |= metrics(y[te], proba) | group_bootstrap_ci(y[te], proba, g[te], args.n_boot)
            rows.append(row)
            print(f"{drug:31s} {model_name:6s} holdout AUROC {row['auroc']:.3f}  "
                  f"(n={row['n']}, R={row['n_resistant']})")

    res = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    res.round(4).to_csv(args.out, index=False)
    Path(args.pred_out).parent.mkdir(parents=True, exist_ok=True)
    pd.concat(preds).to_csv(args.pred_out, sep="\t", index=False)
    print(f"\nwrote {args.out} and {args.pred_out}")

    gap = res[res["validation"].isin(["grouped", "random"])].pivot_table(
        index=["drug", "model"], columns="validation", values="auroc")
    gap["inflation"] = gap["random"] - gap["grouped"]
    print("\n# AUROC: random split minus lineage-grouped split")
    print(gap.round(3).to_string())


if __name__ == "__main__":
    main()
