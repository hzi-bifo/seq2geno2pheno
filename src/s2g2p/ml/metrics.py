"""Per-fold scores and their aggregation over folds.

Scores are those Hu et al. 2024 reported, computed from hard predictions (their AUC too). Fold
aggregates come both plain (mean, sample standard deviation) and support-weighted (weighted
mean, and the weighted variance rescaled by n/(n-1), ``extract_score.weithgted_var``). The paper
used the weighted form for homology-aware (KMA) folds and the plain form for random and
phylogeny-aware folds; both are always written, and the reader picks.
"""

from __future__ import annotations

import math
import statistics

import numpy as np
from sklearn.metrics import (accuracy_score, f1_score, matthews_corrcoef, precision_recall_fscore_support,
                             roc_auc_score)

METRICS = ["f1_macro", "accuracy", "precision_macro", "recall_macro", "f1_positive", "f1_negative",
           "precision_positive", "recall_positive", "precision_negative", "recall_negative", "mcc", "auc"]


def fold_scores(y_true, y_pred) -> dict:
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred, labels=[0, 1], zero_division=0)
    try:
        auc = roc_auc_score(y_true, y_pred)
    except ValueError:  # one class in the fold
        auc = float("nan")
    return {
        "f1_macro": f1_score(y_true, y_pred, average="macro", labels=[0, 1], zero_division=0),
        "accuracy": accuracy_score(y_true, y_pred),
        "precision_macro": float(np.mean(p)), "recall_macro": float(np.mean(r)),
        "f1_positive": f[1], "f1_negative": f[0],
        "precision_positive": p[1], "recall_positive": r[1],
        "precision_negative": p[0], "recall_negative": r[0],
        "mcc": matthews_corrcoef(y_true, y_pred), "auc": auc,
        "support": int(len(y_true)), "support_positive": int(s[1]), "support_negative": int(s[0]),
    }


def _weighted_std(values, mean, weights) -> float:
    n = len(values)
    if n < 2:
        return float("nan")
    var = np.average((np.asarray(values) - mean) ** 2, weights=weights) * n / (n - 1)
    return math.sqrt(var)


def aggregate(rows: list[dict]) -> list[dict]:
    """One row per metric: mean, std, weighted_mean, weighted_std over the folds."""
    support = [r["support"] for r in rows]
    out = []
    for m in METRICS:
        v = [float(r[m]) for r in rows]
        ok = [i for i, x in enumerate(v) if not math.isnan(x)]
        vv, ww = [v[i] for i in ok], [support[i] for i in ok]
        mean = statistics.mean(vv) if vv else float("nan")
        wmean = float(np.average(vv, weights=ww)) if vv else float("nan")
        out.append({"metric": m, "folds": len(vv), "mean": mean,
                    "std": statistics.stdev(vv) if len(vv) > 1 else float("nan"),
                    "weighted_mean": wmean, "weighted_std": _weighted_std(vv, wmean, ww)})
    return out
