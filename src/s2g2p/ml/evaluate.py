"""Nested cross-validation (Geno2Pheno's role, as Hu et al. 2024 ran it in ``main_s2p.py``).

For every outer fold and every candidate (feature set x classifier): a grid search scored by
macro-F1 on inner splits of the training folds, a refit on all training folds, and a test
prediction. Per outer fold, one candidate is then *selected* by its inner score
(``result_analysis.py``), and the selected candidates' test predictions are the result.

Row order is kept deliberately: training rows are the training folds concatenated in fold
order, each in its file order, which is the order the published code fitted on. Random forest
bootstraps depend on it.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import GridSearchCV, GroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Normalizer, StandardScaler

from ..features.compress import DedupColumns
from . import grids
from .metrics import aggregate, fold_scores

BlocksFn = Callable[[int], Mapping[str, pd.DataFrame]]


@dataclass
class Assembled:
    X: np.ndarray
    slices: dict[str, tuple[int, int]]
    columns: dict[str, list[str]]


def feature_sets(ml: dict) -> list[tuple[str, list[str]]]:
    blocks = list(ml["blocks"])
    if ml["feature_sets"] == "concat":
        return [("+".join(blocks), blocks)]
    return [(b, [b]) for b in blocks]


def assemble(blocks: Mapping[str, pd.DataFrame], names: Sequence[str], ids: Sequence[str],
             transforms: Mapping[str, str]) -> Assembled:
    """Stack the named blocks for ``ids``. ``standard_global`` is applied here, over every row
    passed in -- that is the published k-mer scaling, fitted on test genomes too."""
    parts, slices, columns, start = [], {}, {}, 0
    for b in names:
        df = blocks[b]
        missing = [s for s in ids if s not in df.index]
        if missing:
            raise KeyError(f"block {b}: {len(missing)} samples have no row, e.g. {missing[:3]}")
        x = df.loc[list(ids)].to_numpy(dtype=np.float64)
        if transforms[b] == "standard_global":
            x = StandardScaler(with_mean=True, with_std=True).fit_transform(x)
        parts.append(x)
        slices[b] = (start, start + x.shape[1])
        columns[b] = [str(c) for c in df.columns]
        start += x.shape[1]
    return Assembled(np.hstack(parts) if parts else np.zeros((len(ids), 0)), slices, columns)


def preprocessor(slices: Mapping[str, tuple[int, int]], transforms: Mapping[str, str]) -> ColumnTransformer:
    steps = []
    for b, (a, z) in slices.items():
        t = transforms[b]
        tr = {"standard": StandardScaler(with_mean=True, with_std=True), "l1": Normalizer(norm="l1"),
              "dedup": DedupColumns(), "none": "passthrough", "standard_global": "passthrough"}[t]
        steps.append((b, tr, list(range(a, z))))
    return ColumnTransformer(steps, remainder="drop", sparse_threshold=0.0)


def fold_splits(fold_pos: Sequence[np.ndarray]) -> list[tuple[np.ndarray, np.ndarray]]:
    """Leave-one-fold-out splits; training rows are the other folds concatenated in order."""
    return [(np.concatenate([fold_pos[i] for i in range(len(fold_pos)) if i != j]), np.asarray(fold_pos[j]))
            for j in range(len(fold_pos))]


def inner_splits(fold_sizes: Sequence[int], spec: str, y_train: np.ndarray,
                 groups_train: np.ndarray | None) -> list[tuple[np.ndarray, np.ndarray]]:
    """Inner splits as positions into the training rows, which are the training folds
    concatenated in order (``fold_sizes`` gives their lengths)."""
    kind, _, n = str(spec).partition(":")
    if kind == "other_folds":
        offsets = np.cumsum([0] + list(fold_sizes))
        return fold_splits([np.arange(offsets[j], offsets[j + 1]) for j in range(len(fold_sizes))])
    n = int(n or 10)
    if kind == "group_kfold" and groups_train is not None:
        return list(GroupKFold(n_splits=n).split(np.zeros(len(y_train)), y_train, groups_train))
    return list(StratifiedKFold(n_splits=n, shuffle=True, random_state=1).split(np.zeros(len(y_train)), y_train))


def select(cands: list[dict], tie_break: str) -> dict:
    """The candidate with the highest inner mean. ``first``: ties go to the earliest candidate,
    which is what result_analysis.py does in effect (it sorts by (-mean, std) and then looks the
    winning mean up with ``MEAN.index``). ``lowest_std``: ties go to the lowest inner std."""
    best = max(c["inner_mean"] for c in cands)
    tied = [c for c in cands if c["inner_mean"] == best]
    if tie_break == "lowest_std":
        return min(tied, key=lambda c: c["inner_std"])
    return tied[0]


@dataclass
class EvalResult:
    candidates: pd.DataFrame
    folds: pd.DataFrame
    predictions: pd.DataFrame
    summary: pd.DataFrame
    dropped_grid_keys: dict


def evaluate(blocks_for_fold: BlocksFn, labels: pd.Series, folds: Sequence[Sequence[str]], ml: dict,
             groups: pd.Series | None = None, n_jobs: int = 1, log=print) -> EvalResult:
    # Rows in label order, as main_s2p.py built X_all: the cohort-wide k-mer scaling then sums
    # in the published order. Folds index into it by position.
    ids = list(labels.index)
    if sorted(s for f in folds for s in f) != sorted(ids):
        raise ValueError("folds must cover exactly the labelled samples")
    y_all = labels.to_numpy()
    pos = {s: i for i, s in enumerate(ids)}
    fold_pos = [np.array([pos[s] for s in f]) for f in folds]
    g_all = groups.loc[ids].to_numpy() if groups is not None else None
    transforms = ml["transforms"]
    cand_rows, fold_rows, pred_rows, dropped = [], [], [], {}
    for k in range(len(folds)):
        blocks = blocks_for_fold(k)
        train_folds = [fold_pos[j] for j in range(len(folds)) if j != k]
        tr = np.concatenate(train_folds)
        te = fold_pos[k]
        splits = inner_splits([len(f) for f in train_folds], ml["inner_cv"], y_all[tr],
                              g_all[tr] if g_all is not None else None)
        cands = []
        for fs_name, names in feature_sets(ml):
            A = assemble(blocks, names, ids, transforms)
            for clf in ml["classifiers"]:
                est, grid, drop = grids.grid(clf, ml["grid"])
                if drop:
                    dropped[clf] = drop
                pipe = Pipeline([("prep", preprocessor(A.slices, transforms)), ("cl", est)])
                t0 = time.time()
                search = GridSearchCV(pipe, grid, scoring=ml.get("scoring", "f1_macro"), cv=splits,
                                      n_jobs=n_jobs, refit=True, error_score="raise")
                search.fit(A.X[tr], y_all[tr])
                y_pred = search.predict(A.X[te])
                bi = search.best_index_
                prep = search.best_estimator_.named_steps["prep"]
                n_used = int(prep.transform(A.X[te][:1]).shape[1])
                row = {"fold": k, "feature_set": fs_name, "classifier": clf,
                       "inner_mean": float(search.best_score_),
                       "inner_std": float(search.cv_results_["std_test_score"][bi]),
                       "params": json.dumps({p.split("__", 1)[1]: v for p, v in search.best_params_.items()},
                                            default=str, sort_keys=True),
                       "n_features": int(A.X.shape[1]), "n_features_used": n_used,
                       "seconds": round(time.time() - t0, 1)}
                row.update(fold_scores(y_all[te], y_pred))
                row["_pred"] = y_pred
                cands.append(row)
                log(f"fold {k} {fs_name} {clf}: inner {row['inner_mean']:.4f} test f1_macro {row['f1_macro']:.4f}")
        chosen = select(cands, ml["tie_break"])
        for c in cands:
            c["selected"] = c is chosen
            pred = c.pop("_pred")
            if c["selected"]:
                for s, yt, yp in zip([ids[i] for i in te], y_all[te], pred):
                    pred_rows.append({"sample": s, "fold": k, "y_true": int(yt), "y_pred": int(yp),
                                      "feature_set": c["feature_set"], "classifier": c["classifier"]})
                fold_rows.append(dict(c))
            cand_rows.append(c)
    summary = pd.DataFrame(aggregate(fold_rows))
    return EvalResult(pd.DataFrame(cand_rows), pd.DataFrame(fold_rows), pd.DataFrame(pred_rows), summary, dropped)
