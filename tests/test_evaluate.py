"""The nested CV against a near-verbatim copy of the published procedure (main_s2p.py -f_ml +
result_analysis.py extract_best_estimator) on synthetic data, with small grids so it runs in
seconds. The grids themselves are checked separately against hyper_range.py's values."""

import copy
import itertools

import numpy as np
import pandas as pd
import pytest
from sklearn import preprocessing
from sklearn.metrics import f1_score
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline

from s2g2p import presets
from s2g2p.ml import grids
from s2g2p.ml.evaluate import evaluate, select
from s2g2p.ml.model import predict, train

SMALL = {
    "svm": [{"cl__kernel": ["rbf"], "cl__gamma": [1e-2, 1e-3], "cl__C": [0.1, 10]},
            {"cl__kernel": ["linear"], "cl__C": [0.1, 10]}],
    "lr": [{"cl__C": [0.01, 1, 10], "cl__max_iter": [1000]}],
    "rf": [{"cl__n_estimators": [20], "cl__min_samples_split": [2, 5]}],
    "lsvm": [{"cl__C": [0.01, 1], "cl__dual": [True]}],
}


@pytest.fixture
def small_grids(monkeypatch):
    monkeypatch.setitem(grids.GRIDS, "hu2024", SMALL)


def synthetic(seed=0, n=80, n_folds=5):
    rng = np.random.default_rng(seed)
    ids = [f"g{i:03d}" for i in range(n)]
    y = rng.integers(0, 2, n)
    kmer = rng.poisson(20, size=(n, 12)).astype(float)
    kmer[:, 0] += 8 * y
    gpa = rng.integers(0, 2, size=(n, 15))
    gpa[:, 3] = y ^ (rng.random(n) < 0.15)
    gpa[:, 4] = gpa[:, 3]  # a duplicated pattern
    indel = rng.integers(0, 2, size=(n, 6))
    perm = rng.permutation(n)
    folds = [[ids[i] for i in perm[j::n_folds]] for j in range(n_folds)]
    labels = pd.Series(y, index=ids, name="drug")
    blocks = {"kmer": pd.DataFrame(kmer, index=ids, columns=[f"k{i}" for i in range(12)]),
              "gpa": pd.DataFrame(gpa, index=ids, columns=[f"fam{i}" for i in range(15)]),
              "indel": pd.DataFrame(indel, index=ids, columns=[f"fam{i}_indel" for i in range(6)])}
    return labels, folds, blocks


def published(labels, folds, blocks, cl_list):
    """main_s2p.py's -f_ml loop and result_analysis.py's per-fold choice."""
    id_all = np.array(labels.index)
    y_all = labels.to_numpy()
    X_all = blocks["kmer"].loc[id_all]
    scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(X_all)
    X_all = pd.DataFrame(scaler.transform(X_all), index=X_all.index, columns=X_all.columns)
    X_all = pd.concat([X_all, blocks["gpa"].reindex(X_all.index), blocks["indel"].reindex(X_all.index)], axis=1)
    p = {s: i for i, s in enumerate(id_all)}
    folders_index = [[p[s] for s in f] for f in folds]

    def create_generator(nFolds):
        for idx in range(len(nFolds)):
            test = nFolds[idx]
            train = list(itertools.chain(*[fold for idy, fold in enumerate(nFolds) if idy != idx]))
            yield train, test

    score = {}
    for chosen_cl in cl_list:
        cl, hyper_space = grids.estimator(chosen_cl), copy.deepcopy(SMALL[chosen_cl])
        per = []
        for out_cv in range(len(folds)):
            test_idx = folders_index[out_cv]
            tvt = folders_index[:out_cv] + folders_index[out_cv + 1:]
            id_vt = id_all[list(itertools.chain.from_iterable(tvt))]
            y_vt = y_all[list(itertools.chain.from_iterable(tvt))]
            search = GridSearchCV(estimator=Pipeline(steps=[("cl", cl)]), param_grid=hyper_space, n_jobs=1,
                                  scoring="f1_macro", cv=create_generator(tvt), refit=True)
            search.fit(X_all, y_all)
            pipe = search.best_estimator_
            pipe.fit(X_all.loc[id_vt, :], y_vt)
            y_pred = pipe.predict(X_all.loc[id_all[test_idx], :])
            per.append({"mean": search.best_score_,
                        "std": search.cv_results_["std_test_score"][search.best_index_],
                        "f1": f1_score(y_all[test_idx], y_pred, average="macro"),
                        "pred": y_pred.tolist()})
        score[chosen_cl] = per
    chosen = []
    for out_cv in range(len(folds)):
        MEAN = [score[c][out_cv]["mean"] for c in cl_list]
        STD = [score[c][out_cv]["std"] for c in cl_list]
        combined = sorted(zip(MEAN, STD), key=lambda x: (-x[0], x[1]))
        cl = cl_list[MEAN.index(combined[0][0])]
        chosen.append((cl, score[cl][out_cv]["f1"], score[cl][out_cv]["pred"]))
    return chosen


def test_hu2024_procedure_matches_published(small_grids):
    labels, folds, blocks = synthetic()
    cfg = presets.load("hu2024")
    res = evaluate(lambda k: blocks, labels, folds, cfg["ml"], log=lambda m: None)
    ref = published(labels, folds, blocks, cfg["ml"]["classifiers"])
    got = res.folds.sort_values("fold")
    assert got["classifier"].tolist() == [c for c, _, _ in ref]
    assert np.allclose(got["f1_macro"].to_numpy(), [f for _, f, _ in ref])
    for k, (_, _, pred) in enumerate(ref):
        p = res.predictions[res.predictions.fold == k]
        assert p["y_pred"].tolist() == pred


def test_noleak_runs_and_dedups_in_cv(small_grids):
    labels, folds, blocks = synthetic(1)
    cfg = presets.load("noleak")
    res = evaluate(lambda k: blocks, labels, folds, cfg["ml"], log=lambda m: None)
    assert len(res.folds) == len(folds)
    assert (res.candidates["n_features_used"] < res.candidates["n_features"]).all()
    assert set(res.summary["metric"]) >= {"f1_macro", "accuracy", "mcc"}


def test_tie_break():
    c = [{"inner_mean": 0.8, "inner_std": 0.2, "n": 1}, {"inner_mean": 0.8, "inner_std": 0.1, "n": 2},
         {"inner_mean": 0.7, "inner_std": 0.0, "n": 3}]
    assert select(c, "first")["n"] == 1
    assert select(c, "lowest_std")["n"] == 2


def test_train_predict_round_trip(small_grids, tmp_path):
    labels, folds, blocks = synthetic(2)
    cfg = presets.load("noleak")
    bundle = train(blocks, labels, folds, cfg, log=lambda m: None)
    out = predict(bundle, blocks, list(labels.index[:5]))
    assert out["prediction"].isin([0, 1]).all() and len(out) == 5
    bad = dict(blocks)
    bad["gpa"] = blocks["gpa"].iloc[:, 1:]
    with pytest.raises(ValueError):
        predict(bundle, bad)
    with pytest.raises(ValueError):
        train(blocks, labels, folds, presets.load("hu2024"))


def test_published_grid_values():
    g = grids.HU2024
    assert g["svm"][0]["cl__C"] == [0.001, 0.10, 0.1, 10, 25, 50, 100, 1000]
    assert g["lr"][0]["cl__tol"] == [1e-06, 1e-03, 1e-04]
    assert g["rf"][0]["cl__n_estimators"] == [100, 200, 500, 1000]
    assert sum(len(list(itertools.product(*d.values()))) for d in g["svm"]) == 72
    for clf in ("svm", "lr", "rf", "lsvm"):
        grids.grid(clf, "hu2024")  # every estimator accepts its grid on this scikit-learn


def test_presets_validate():
    for name in presets.available():
        presets.load(name)
