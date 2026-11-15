"""SHAP analysis: which AMR elements drive each drug's XGBoost model, and do they
match the drug class that NCBI's AMRFinderPlus database assigns to them?

For each drug an XGBoost model is fitted on all genomes (after the same
prevalence filter used in training). SHAP values come from shap.TreeExplainer;
if that fails on this xgboost/shap version pair, XGBoost's built-in
pred_contribs (the same TreeSHAP algorithm) is used instead.

Reading presence/absence SHAP values: a positive mean SHAP among genomes that
CARRY an element means its presence pushes the prediction towards resistant.
`ncbi_class_match` only says whether AMRFinderPlus's own class/subclass
annotation names the drug's class; every gene-drug link still has to be checked
by hand in CARD and the NCBI Reference Gene Catalog.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.train import DRUGS, MinPrevalence, load_inputs, make_model

# Keywords that AMRFinderPlus class/subclass annotations use for each drug's class.
DRUG_KEYWORDS = {
    "meropenem": ["CARBAPENEM"],
    "ciprofloxacin": ["QUINOLONE"],
    "gentamicin": ["GENTAMICIN"],
    "trimethoprim/sulfamethoxazole": ["TRIMETHOPRIM", "SULFONAMIDE"],
}


def shap_matrix(model, X: np.ndarray) -> tuple[np.ndarray, str]:
    """(n_genomes, n_features) SHAP values in log-odds units, plus the method used."""
    try:
        import shap
        vals = shap.TreeExplainer(model).shap_values(X)
        vals = np.asarray(vals[1] if isinstance(vals, list) else vals)
        if vals.ndim == 3:
            vals = vals[:, :, 1]
        return vals, f"shap {shap.__version__} TreeExplainer"
    except Exception as err:  # version mismatch between shap and xgboost
        import xgboost as xgb
        contribs = model.get_booster().predict(xgb.DMatrix(X), pred_contribs=True)
        return contribs[:, :-1], f"xgboost pred_contribs (shap failed: {type(err).__name__})"


def drug_table(drug, X_df, y, ann, top) -> tuple[pd.DataFrame, str]:
    X = X_df.values.astype(np.int8)
    filt = MinPrevalence().fit(X)
    Xf = filt.transform(X)
    names = X_df.columns[filt.keep_]
    model = make_model("xgb", y).fit(Xf, y)
    vals, method = shap_matrix(model, Xf)

    present = Xf.astype(bool)
    with np.errstate(invalid="ignore"):
        mean_when_present = np.nanmean(np.where(present, vals, np.nan), axis=0)
    tab = pd.DataFrame({
        "drug": drug,
        "element": names,
        "mean_abs_shap": np.abs(vals).mean(axis=0),
        "mean_shap_when_present": mean_when_present,
        "prev_in_R": Xf[y == 1].mean(axis=0),
        "prev_in_S": Xf[y == 0].mean(axis=0),
    })
    tab = tab.merge(ann[["subtype", "amr_class", "amr_subclass"]],
                    left_on="element", right_index=True, how="left")
    text = (tab["amr_class"].fillna("") + "/" + tab["amr_subclass"].fillna("")).str.upper()
    tab["ncbi_class_match"] = text.apply(
        lambda s: any(k in s for k in DRUG_KEYWORDS[drug]))
    tab["direction"] = np.where(tab["mean_shap_when_present"] > 0,
                                "towards R", "towards S")
    tab["verified_in_CARD"] = ""     # fill in by hand
    tab = tab.sort_values("mean_abs_shap", ascending=False).head(top)
    tab.insert(1, "rank", range(1, len(tab) + 1))
    return tab, method


def plot(tables: pd.DataFrame, path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, drug in zip(axes.flat, DRUGS):
        t = tables[tables["drug"] == drug].head(12).iloc[::-1]
        colors = ["#2a6f97" if m else "#bbbbbb" for m in t["ncbi_class_match"]]
        ax.barh(t["element"], t["mean_abs_shap"], color=colors)
        ax.set_title(drug)
        ax.set_xlabel("mean |SHAP| (log-odds)")
    fig.suptitle("Top XGBoost features per drug (blue = AMRFinderPlus annotates "
                 "the drug's class)")
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", default="data/features/amr_matrix.tsv")
    p.add_argument("--groups", default="data/features/lineage_groups.tsv")
    p.add_argument("--labels", default="data/metadata/selected_500.tsv")
    p.add_argument("--annotation", default="data/features/element_annotation.tsv")
    p.add_argument("--top", type=int, default=15)
    p.add_argument("--out", default="reports/shap_top_features.csv")
    p.add_argument("--figure", default="reports/figures/shap_top_features.png")
    args = p.parse_args()

    X_df, _, lab = load_inputs(args.matrix, args.groups, args.labels)
    ann = pd.read_csv(args.annotation, sep="\t", index_col="element")
    tables = []
    for drug in DRUGS:
        tab, method = drug_table(drug, X_df, lab[drug].values.astype(int), ann, args.top)
        tables.append(tab)
        k = min(10, len(tab))
        hit = int(tab["ncbi_class_match"].head(k).sum())
        print(f"\n### {drug}  ({method})  top-{k} with matching NCBI class: {hit}/{k}")
        print(tab[["rank", "element", "mean_abs_shap", "direction", "prev_in_R",
                   "prev_in_S", "amr_subclass", "ncbi_class_match"]]
              .head(10).round(3).to_string(index=False))
    out = pd.concat(tables, ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.round(4).to_csv(args.out, index=False)
    plot(out, args.figure)
    print(f"\nwrote {args.out} and {args.figure}")


if __name__ == "__main__":
    main()
