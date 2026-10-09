"""Redundancy compression of binary feature tables (Seq2Geno's ``featCompress.py``).

Columns with the same pattern over the samples collapse to one representative: the first such
column in table order, renamed ``<representative>|<group size>``. Groups come out in sorted
pattern order, as pandas' ``groupby`` produced them in the original.

Where this is fitted is what separates the presets. The published pipeline compressed once
over the whole species cohort (test genomes included); :class:`DedupColumns` does the same
thing as a scikit-learn step, so inside a pipeline it is refitted on each training split.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin

# featCompress's value order: the original grouped on the strings '0' < '1' < 'NA'.
_CODES = {"0": 0, "1": 1, "NA": 2}


def _pattern_codes(df: pd.DataFrame) -> np.ndarray:
    """features x samples codes. An empty cell is an error, as featCompress's check_input made
    it (it filled NaN with '.' and rejected anything but '0', '1' and 'NA')."""
    if df.isna().to_numpy().any():
        r, c = np.argwhere(df.isna().to_numpy())[0]
        raise ValueError(f"empty cell at {df.index[r]} / {df.columns[c]}; only 0, 1 and NA are allowed")
    s = df.astype(str)
    bad = ~s.isin(list(_CODES)).to_numpy()
    if bad.any():
        r, c = np.argwhere(bad)[0]
        raise ValueError(f"only 0, 1 and NA are allowed; found {s.iat[r, c]!r} at "
                         f"{s.index[r]} / {s.columns[c]}")
    return np.vectorize(_CODES.__getitem__, otypes=[np.int8])(s.to_numpy()).T


def compress(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """featCompress on a samples x features table: (compressed table, groups)."""
    if df.shape[1] == 0:
        return df.copy(), {}
    codes = _pattern_codes(df)
    _, first, inverse, counts = np.unique(codes, axis=0, return_index=True,
                                          return_inverse=True, return_counts=True)
    inverse = inverse.reshape(-1)
    names = df.columns
    groups: dict[str, list[str]] = {}
    out_cols, out_names = [], []
    for g, (f, n) in enumerate(zip(first, counts)):
        rep = str(names[f])
        new = f"{rep}|{n}"
        members = names[np.flatnonzero(inverse == g)].tolist()
        groups[new] = [str(m) for m in members]
        out_cols.append(f)
        out_names.append(new)
    out = df.iloc[:, out_cols].copy()
    out.columns = out_names
    return out, groups


def write_groups(groups: dict[str, list[str]], path) -> None:
    """The ``*_GROUPS`` file featCompress wrote beside the compressed table."""
    with open(path, "w") as fh:
        for rep, members in groups.items():
            fh.write("{}:\t{}\n".format(rep, "\t".join(members)))


class DedupColumns(BaseEstimator, TransformerMixin):
    """featCompress fitted on the rows it is given: keeps one column per distinct pattern
    (the first carrying it), in sorted pattern order."""

    def fit(self, X, y=None):
        X = np.asarray(X)
        if X.shape[1] == 0:
            self.keep_ = np.zeros(0, dtype=int)
            return self
        _, self.keep_ = np.unique(X.T, axis=0, return_index=True)
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.keep_]
