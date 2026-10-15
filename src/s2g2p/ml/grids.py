"""Classifiers and hyperparameter grids.

``hu2024`` is ``AMR_software/seq2geno/hyper_range.py`` of Hu et al. 2024 (hzi-bifo/
AMR_benchmarking @ 2850078), copied value for value, duplicates included (``0.10`` and ``0.1``
are the same C; GridSearchCV fits both, which changes nothing but cost). Two translations are
needed to run on a current scikit-learn and are the only differences:

* random forest ``max_features="auto"`` was removed in scikit-learn 1.3; for classifiers it
  meant ``"sqrt"``, which is what is passed.
* grid keys an installed scikit-learn no longer accepts are dropped when they only restate a
  default (recorded in the run's ``run.json``).

``g2p`` is NOT Geno2Pheno's grid: Geno2Pheno's server-side grids are not in its public client.
It is the small grid used for this project's Geno2Pheno-protocol rerun, kept so that preset is
runnable; replace it with the real one before reporting it as Geno2Pheno.
"""

from __future__ import annotations

import copy

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, LinearSVC

_C_SVM = [0.001, 0.10, 0.1, 10, 25, 50, 100, 1000]
_C_LIN = [0.001, 0.01, 0.015, 0.025, 0.10, 0.1, 10, 25, 50, 100, 1000]
_GAMMA = [1e-2, 1e-3, 1e-4, 1e-5]

HU2024 = {
    "svm": [
        {"cl__kernel": ["rbf"], "cl__gamma": _GAMMA, "cl__C": _C_SVM},
        {"cl__kernel": ["sigmoid"], "cl__gamma": _GAMMA, "cl__C": _C_SVM},
        {"cl__kernel": ["linear"], "cl__C": _C_SVM},
    ],
    "lr": [
        {"cl__C": _C_LIN, "cl__tol": [1e-06, 1e-03, 1e-04], "cl__dual": [False],
         "cl__fit_intercept": [True], "cl__intercept_scaling": [1], "cl__max_iter": [1000],
         "cl__penalty": ["l2"]},
    ],
    "lsvm": [
        {"cl__C": _C_LIN, "cl__tol": [1e-06, 1e-03], "cl__dual": [True], "cl__intercept_scaling": [1],
         "cl__loss": ["squared_hinge"], "cl__multi_class": ["ovr"], "cl__penalty": ["l2"]},
        {"cl__C": _C_LIN, "cl__tol": [1e-06, 1e-03], "cl__dual": [False], "cl__fit_intercept": [True],
         "cl__intercept_scaling": [1], "cl__loss": ["squared_hinge"], "cl__multi_class": ["ovr"],
         "cl__penalty": ["l1"]},
    ],
    "rf": [
        {"cl__n_estimators": [100, 200, 500, 1000], "cl__criterion": ["entropy", "gini"],
         "cl__max_features": ["sqrt"],  # published: "auto" (== "sqrt" for classifiers)
         "cl__min_samples_split": [2, 5, 10], "cl__min_samples_leaf": [1],
         "cl__class_weight": ["balanced", None]},
    ],
}

G2P = {
    "lsvm": [{"cl__C": [0.01, 0.1, 1.0, 10.0]}],
    "svm": [{"cl__C": [0.1, 1.0, 10.0], "cl__gamma": ["scale", 1e-3]}],
    "lr": [{"cl__C": [0.01, 0.1, 1.0, 10.0]}],
    "rf": [{"cl__n_estimators": [200, 500], "cl__max_depth": [None, 20]}],
}

GRIDS = {"hu2024": HU2024, "g2p": G2P}


def estimator(name: str, flavour: str = "hu2024"):
    """The base estimator of ``hyper_range.py`` (``random_state=1``, balanced classes)."""
    if flavour == "hu2024":
        return {
            "svm": lambda: SVC(random_state=1, class_weight="balanced"),
            "lr": lambda: LogisticRegression(random_state=1, class_weight="balanced"),
            "lsvm": lambda: LinearSVC(random_state=1, max_iter=10000, class_weight="balanced"),
            "rf": lambda: RandomForestClassifier(random_state=1),
        }[name]()
    if flavour == "g2p":
        return {
            "svm": lambda: SVC(random_state=1, class_weight="balanced"),
            "lr": lambda: LogisticRegression(random_state=1, class_weight="balanced", max_iter=5000),
            "lsvm": lambda: LinearSVC(random_state=1, max_iter=10000, class_weight="balanced"),
            "rf": lambda: RandomForestClassifier(random_state=1, class_weight="balanced"),
        }[name]()
    raise KeyError(f"unknown grid flavour {flavour!r}")


def grid(name: str, flavour: str = "hu2024") -> tuple[object, list[dict], list[str]]:
    """(estimator, grid, dropped keys). Keys the installed estimator does not have are dropped."""
    est = estimator(name, flavour)
    valid = set(est.get_params())
    out, dropped = [], set()
    for g in copy.deepcopy(GRIDS[flavour][name]):
        keep = {}
        for k, v in g.items():
            if k.split("__", 1)[1] in valid:
                keep[k] = v
            else:
                dropped.add(k)
        out.append(keep)
    return est, out, sorted(dropped)
