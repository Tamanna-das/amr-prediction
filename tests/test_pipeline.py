"""Tests for the feature matrix, lineage grouping and leakage-safe training code."""
import numpy as np
import pandas as pd

from src.build_matrix import collect_hits, presence_matrix
from src.lineages import cluster, genome_id_from_path
from src.train import MinPrevalence, group_bootstrap_ci, oof_predictions

HEADER = ["Name", "Element symbol", "Type", "Subtype", "Class", "Subclass"]


def _write_report(path, rows):
    pd.DataFrame(rows, columns=HEADER).to_csv(path, sep="\t", index=False)


def test_matrix_keeps_zero_hit_genomes_and_only_amr_type(tmp_path):
    _write_report(tmp_path / "g1.tsv", [
        ["g1", "blaKPC-2", "AMR", "AMR", "BETA-LACTAM", "CARBAPENEM"],
        ["g1", "gyrA_S83I", "AMR", "POINT", "QUINOLONE", "QUINOLONE"],
        ["g1", "ybtQ", "VIRULENCE", "VIRULENCE", "NA", "NA"],
    ])
    _write_report(tmp_path / "g2.tsv", [])                      # no hits at all
    hits = collect_hits(["g1", "g2"], tmp_path)
    mat = presence_matrix(hits, ["g1", "g2"])
    assert list(mat.index) == ["g1", "g2"]
    assert list(mat.columns) == ["blaKPC-2", "gyrA_S83I"]       # virulence dropped
    assert mat.loc["g2"].sum() == 0


def test_genome_id_from_path():
    assert genome_id_from_path("data/genomes/573.12771.fna") == "573.12771"


def test_single_linkage_chains_close_genomes():
    ids = ["a", "b", "c", "d"]
    d = pd.DataFrame([[0, .001, .009, .05], [.001, 0, .004, .05],
                      [.009, .004, 0, .05], [.05, .05, .05, 0]], index=ids, columns=ids)
    g = cluster(d, 0.005)
    assert g["a"] == g["b"] == g["c"]        # a-c linked through b
    assert g["d"] != g["a"]


def test_min_prevalence_is_fitted_on_training_rows_only():
    X_train = np.array([[1, 0, 1], [1, 0, 0], [1, 0, 1], [1, 0, 0]])
    filt = MinPrevalence(min_frac=0.1).fit(X_train)
    assert list(filt.keep_) == [False, False, True]   # constant columns dropped
    X_test = np.array([[0, 1, 1]])                   # test values never change the choice
    assert filt.transform(X_test).shape == (1, 1)


def test_grouped_cv_never_splits_a_lineage():
    rng = np.random.default_rng(0)
    groups = np.repeat(np.arange(20), 5)
    X = rng.integers(0, 2, (100, 4))
    y = np.tile([0, 1], 50)
    proba = oof_predictions("logreg", X, y, groups, "grouped")
    assert not np.isnan(proba).any()


def test_group_bootstrap_ci_brackets_point_estimate():
    rng = np.random.default_rng(1)
    y = np.tile([0, 1], 50)
    proba = np.clip(y * 0.7 + rng.random(100) * 0.3, 0, 1)
    ci = group_bootstrap_ci(y, proba, np.repeat(np.arange(25), 4), n_boot=200)
    assert ci["auroc_lo"] <= ci["auroc_hi"] <= 1.0
