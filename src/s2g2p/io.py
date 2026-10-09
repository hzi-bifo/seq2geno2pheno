"""Inputs and outputs shared by every stage: sample sheets, labels, folds, feature matrices.

Feature matrices are TSV, one row per sample and one column per feature, with the sample id in
the first column. Legacy Seq2Geno tables (index label ``Isolate`` / ``Isolates``) read the same.
"""

from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd

RESERVED_COLUMNS = ("sample", "assembly", "group")
_ID_RE = re.compile(r"^[A-Za-z0-9._\-]+$")
# Phenotype spellings accepted in a sample sheet. Anything else is an error, never a guess.
LABEL_VALUES = {
    "1": 1, "r": 1, "resistant": 1, "nonsusceptible": 1, "non-susceptible": 1,
    "0": 0, "s": 0, "susceptible": 0,
}


def _sep_for(path: str | os.PathLike) -> str:
    return "," if str(path).lower().endswith(".csv") else "\t"


def read_samplesheet(path: str | os.PathLike) -> pd.DataFrame:
    """The run's cohort: one row per genome.

    Columns: ``sample`` (required, unique), ``assembly`` (path to a nucleotide FASTA; relative
    paths resolve against the sheet's directory), optional ``group`` (for group-aware inner CV),
    and any number of phenotype columns. An empty phenotype cell means "in the cohort, not
    labelled for this phenotype": the genome still joins the pangenome but not that model.
    """
    df = pd.read_csv(path, sep=_sep_for(path), dtype=str, keep_default_na=False)
    if "sample" not in df.columns:
        raise ValueError(f"{path}: no 'sample' column")
    bad = [s for s in df["sample"] if not _ID_RE.match(s)]
    if bad:
        raise ValueError(f"{path}: sample ids must match {_ID_RE.pattern}; e.g. {bad[:3]}")
    dup = df["sample"][df["sample"].duplicated()].tolist()
    if dup:
        raise ValueError(f"{path}: duplicated sample ids, e.g. {dup[:3]}")
    if "assembly" in df.columns:
        base = Path(path).resolve().parent
        df["assembly"] = [str(p if os.path.isabs(p) else base / p) if p else "" for p in df["assembly"]]
    return df.set_index("sample")


def phenotype_columns(sheet: pd.DataFrame) -> list[str]:
    return [c for c in sheet.columns if c not in RESERVED_COLUMNS]


def labels_for(sheet: pd.DataFrame, phenotype: str) -> pd.Series:
    """0/1 labels of the samples labelled for ``phenotype``, in sample-sheet order."""
    if phenotype not in sheet.columns:
        raise KeyError(f"no phenotype column {phenotype!r}; have {phenotype_columns(sheet)}")
    raw = sheet[phenotype].str.strip()
    raw = raw[raw != ""]
    out = {}
    for sid, v in raw.items():
        key = v.lower()
        if key not in LABEL_VALUES:
            raise ValueError(f"{phenotype}: sample {sid} has label {v!r}; use 0/1 or S/R")
        out[sid] = LABEL_VALUES[key]
    return pd.Series(out, name=phenotype, dtype=int)


def read_labels(path: str | os.PathLike, phenotype: str | None = None) -> pd.Series:
    """Labels from a table with a ``sample`` column; ``phenotype`` names the label column
    (default: the only non-reserved column)."""
    sheet = read_samplesheet(path)
    cols = phenotype_columns(sheet)
    if phenotype is None:
        if len(cols) != 1:
            raise ValueError(f"{path}: several phenotype columns {cols}; name one")
        phenotype = cols[0]
    return labels_for(sheet, phenotype)


def read_folds(path: str | os.PathLike) -> list[list[str]]:
    """Outer folds, either the benchmark's JSON (a list of id lists) or a TSV with
    ``sample`` and ``fold`` columns (fold numbers 0..n-1)."""
    p = str(path)
    if p.endswith(".json"):
        with open(p) as fh:
            folds = json.load(fh)
        if not isinstance(folds, list) or not all(isinstance(f, list) for f in folds):
            raise ValueError(f"{p}: expected a JSON list of lists of sample ids")
        return [[str(s) for s in f] for f in folds]
    df = pd.read_csv(p, sep=_sep_for(p), dtype=str)
    if not {"sample", "fold"} <= set(df.columns):
        raise ValueError(f"{p}: needs 'sample' and 'fold' columns")
    k = df["fold"].astype(int)
    if sorted(k.unique()) != list(range(k.max() + 1)):
        raise ValueError(f"{p}: fold numbers must be 0..n-1 without gaps")
    return [df.loc[k == i, "sample"].tolist() for i in range(k.max() + 1)]


def write_folds(folds: Sequence[Sequence[str]], path: str | os.PathLike) -> None:
    with open(path, "w") as fh:
        fh.write("sample\tfold\n")
        for i, f in enumerate(folds):
            for s in f:
                fh.write(f"{s}\t{i}\n")


def read_matrix(path: str | os.PathLike, dtype=None) -> pd.DataFrame:
    """A samples x features table. Values are kept as read unless ``dtype`` is given."""
    with open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")[0]
    # The id column is always text: read as a number, '1352.11330' would become '1352.1133'.
    # It is read as an ordinary column and then made the index, because older pandas (1.2)
    # ignores a dtype given for the index_col column.
    df = pd.read_csv(path, sep="\t", dtype=str if dtype is None else {first: str},
                     keep_default_na=False, na_values=[""])
    df = df.set_index(first)
    df.index = df.index.astype(str)
    df.index.name = "sample"
    if dtype is not None:
        df = df.astype(dtype)
    return df


def write_matrix(df: pd.DataFrame, path: str | os.PathLike) -> None:
    df.to_csv(path, sep="\t", index=True, index_label="sample")


def expand_fold_path(template: str, fold: int) -> str:
    """``'fold_{fold}/gpa.tsv'`` -> ``'fold_3/gpa.tsv'``; a path without ``{fold}`` is shared by
    all folds."""
    return template.replace("{fold}", str(fold))


def write_manifest(directory: str | os.PathLike, rows: Iterable[Mapping], columns: Sequence[str]) -> Path:
    """``manifest.tsv`` beside the files it indexes, one row per file."""
    path = Path(directory) / "manifest.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            w.writerow([r.get(c, "") for c in columns])
    return path


def read_fasta(path: str | os.PathLike):
    """(header without '>', sequence) pairs; the header is kept whole."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)
