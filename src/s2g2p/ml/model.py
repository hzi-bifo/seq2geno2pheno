"""A final model for new genomes: train on every labelled genome, predict placed ones.

Model choice follows the evaluation's rule with the outer folds as the selection CV: every
candidate is grid-searched over the folds, the best by mean macro-F1 is refitted on all labelled
genomes. Only presets with ``predict: true`` qualify -- a cohort-fitted feature (``hu2024``'s
k-mer scaling, pangenome and compression) has no definition for a genome outside the cohort.

The saved bundle is a joblib pickle: load only bundles you made or trust.
"""

from __future__ import annotations

import datetime
import platform
from typing import Mapping, Sequence

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline

from .. import __version__
from . import grids
from .evaluate import assemble, feature_sets, fold_splits, inner_splits, preprocessor, select

BUNDLE_FORMAT = 1


def train(blocks: Mapping[str, pd.DataFrame], labels: pd.Series, folds: Sequence[Sequence[str]], cfg: dict,
          groups: pd.Series | None = None, n_jobs: int = 1, log=print) -> dict:
    if not cfg.get("predict"):
        raise ValueError(f"preset {cfg['name']} fits features on the whole cohort and cannot predict "
                         "unseen genomes; use a preset with predict: true (e.g. noleak)")
    ml = cfg["ml"]
    ids = list(labels.index)
    pos = {s: i for i, s in enumerate(ids)}
    y = labels.to_numpy()
    if str(ml["inner_cv"]).startswith("other_folds"):
        splits = fold_splits([np.array([pos[s] for s in f]) for f in folds])
    else:
        splits = inner_splits([len(ids)], ml["inner_cv"], y,
                              groups.loc[ids].to_numpy() if groups is not None else None)
    cands = []
    for fs_name, names in feature_sets(ml):
        A = assemble(blocks, names, ids, ml["transforms"])
        for clf in ml["classifiers"]:
            est, grid, _ = grids.grid(clf, ml["grid"])
            pipe = Pipeline([("prep", preprocessor(A.slices, ml["transforms"])), ("cl", est)])
            search = GridSearchCV(pipe, grid, scoring=ml.get("scoring", "f1_macro"), cv=splits,
                                  n_jobs=n_jobs, refit=True, error_score="raise")
            search.fit(A.X, y)
            cands.append({"feature_set": fs_name, "blocks": names, "classifier": clf,
                          "inner_mean": float(search.best_score_),
                          "inner_std": float(search.cv_results_["std_test_score"][search.best_index_]),
                          "params": search.best_params_, "pipeline": search.best_estimator_,
                          "columns": A.columns})
            log(f"{fs_name} {clf}: cv f1_macro {search.best_score_:.4f}")
    best = select(cands, ml["tie_break"])
    return {
        "format": BUNDLE_FORMAT, "s2g2p_version": __version__, "sklearn_version": sklearn.__version__,
        "python": platform.python_version(), "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "preset": cfg["name"], "preset_config": cfg, "phenotype": labels.name,
        "training_samples": ids, "feature_set": best["feature_set"], "blocks": best["blocks"],
        "columns": best["columns"], "classifier": best["classifier"], "params": best["params"],
        "cv_f1_macro": best["inner_mean"], "pipeline": best["pipeline"],
        "candidates": [{k: v for k, v in c.items() if k not in ("pipeline", "columns", "blocks")} for c in cands],
    }


def save(bundle: dict, path) -> None:
    joblib.dump(bundle, path, compress=3)


def load(path) -> dict:
    b = joblib.load(path)
    if not isinstance(b, dict) or b.get("format") != BUNDLE_FORMAT:
        raise ValueError(f"{path} is not an s2g2p model bundle of format {BUNDLE_FORMAT}")
    return b


def predict(bundle: dict, blocks: Mapping[str, pd.DataFrame], samples: Sequence[str] | None = None) -> pd.DataFrame:
    """Predictions for the samples of the given blocks. Their columns must be the training
    columns (the placement step guarantees it); any difference is an error, not a fill."""
    names = bundle["blocks"]
    if samples is None:
        samples = list(blocks[names[0]].index)
    aligned = {}
    for b in names:
        want = bundle["columns"][b]
        have = [str(c) for c in blocks[b].columns]
        if have != want:
            extra, missing = set(have) - set(want), set(want) - set(have)
            raise ValueError(f"block {b}: columns differ from training ({len(missing)} missing, "
                             f"{len(extra)} extra); place the genomes into the training pangenome")
        aligned[b] = blocks[b]
    A = assemble(aligned, names, samples, bundle["preset_config"]["ml"]["transforms"])
    pipe = bundle["pipeline"]
    out = pd.DataFrame({"sample": list(samples), "prediction": pipe.predict(A.X).astype(int)})
    if hasattr(pipe, "predict_proba"):
        try:
            out["probability_resistant"] = pipe.predict_proba(A.X)[:, 1]
        except AttributeError:  # SVC without probability=True
            pass
    if hasattr(pipe, "decision_function") and "probability_resistant" not in out:
        out["decision"] = np.asarray(pipe.decision_function(A.X)).ravel()
    out["phenotype"] = bundle["phenotype"]
    out["classifier"] = bundle["classifier"]
    return out
