"""Outer cross-validation folds.

The published evaluation used fixed fold files per (species, antibiotic) -- random, phylogeny-
aware and homology-aware (KMA) partitions -- and every preset accepts such a file unchanged.
For a user's own cohort without one, :func:`stratified_folds` gives a seeded random split.
Phylogeny- and homology-aware fold construction is not part of v3.0 (see ROADMAP.md).
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd


def stratified_folds(labels: pd.Series, n_folds: int = 10, seed: int = 0) -> list[list[str]]:
    """Seeded stratified random folds over the labelled samples."""
    from sklearn.model_selection import StratifiedKFold

    if labels.value_counts().min() < n_folds:
        raise ValueError(f"a class has fewer than {n_folds} samples; use fewer folds")
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    ids = labels.index.to_numpy()
    return [ids[test].tolist() for _, test in skf.split(np.zeros(len(ids)), labels.to_numpy())]


def restrict(folds: Sequence[Sequence[str]], samples: Sequence[str]) -> list[list[str]]:
    """Folds limited to ``samples``; a fold file may list more genomes than are labelled."""
    keep = set(samples)
    return [[s for s in f if s in keep] for f in folds]


def validate(folds: Sequence[Sequence[str]], samples: Sequence[str]) -> None:
    """Folds must be disjoint, non-empty and cover exactly ``samples``."""
    seen: dict[str, int] = {}
    for i, f in enumerate(folds):
        if not f:
            raise ValueError(f"fold {i} is empty")
        for s in f:
            if s in seen:
                raise ValueError(f"sample {s} is in folds {seen[s]} and {i}")
            seen[s] = i
    want = set(samples)
    missing = want - set(seen)
    extra = set(seen) - want
    if missing:
        raise ValueError(f"{len(missing)} labelled samples are in no fold, e.g. {sorted(missing)[:3]}")
    if extra:
        raise ValueError(f"{len(extra)} fold members are not labelled samples, e.g. {sorted(extra)[:3]}")
