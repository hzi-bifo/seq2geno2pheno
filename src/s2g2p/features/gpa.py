"""Gene presence/absence from a Roary run (Seq2Geno's ``gpa_bin_mat`` rule).

A sample has family F when Roary's ``gene_presence_absence.csv`` names any of its genes in F's
row. Columns keep the CSV's row order, which later steps depend on: redundancy compression
keeps the first column of each pattern.
"""

from __future__ import annotations

import os
from collections import OrderedDict
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

# Roary's annotation columns ahead of the per-isolate columns (Roary 3.x).
ROARY_META_COLUMNS = (
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
    "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
    "Accessory Fragment", "Accessory Order with Fragment", "QC", "Min group size nuc",
    "Max group size nuc", "Avg group size nuc",
)


def read_roary_gpa(path: str | os.PathLike) -> pd.DataFrame:
    return pd.read_csv(path, sep=",", header=0, index_col=0, quotechar='"', low_memory=False, dtype=str)


def roary_isolates(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in ROARY_META_COLUMNS]


def gpa_from_roary(path: str | os.PathLike, samples: Sequence[str] | None = None,
                   rename: Mapping[str, str] | None = None) -> pd.DataFrame:
    """samples x families, 0/1 ints. ``samples`` defaults to every isolate column; ``rename``
    maps Roary family names to output names (Seq2Geno's ``roary_abricate.txt``)."""
    df = read_roary_gpa(path)
    cols = list(samples) if samples is not None else roary_isolates(df)
    missing = [s for s in cols if s not in df.columns]
    if missing:
        raise KeyError(f"{path}: {len(missing)} samples have no column, e.g. {missing[:3]}")
    present = df.loc[:, cols].notna().to_numpy().T.astype(np.int8)
    names = [rename.get(f, f) for f in df.index] if rename else list(df.index)
    out = pd.DataFrame(present, index=cols, columns=names)
    out.index.name = "sample"
    return out


def read_rename_map(path: str | os.PathLike) -> dict[str, str]:
    out = {}
    with open(path) as fh:
        for line in fh:
            a = line.strip().split("\t")
            if len(a) >= 2:
                out[a[0]] = a[1]
    return out


def read_clustered_proteins(path: str | os.PathLike) -> "OrderedDict[str, list[str]]":
    """Roary's ``clustered_proteins``: family -> member locus tags, in file order."""
    fam2members: OrderedDict[str, list[str]] = OrderedDict()
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            fam, first = (x.strip() for x in fields[0].split(":", 1))
            fam2members[fam] = [first] + fields[1:]
    return fam2members


def protein_to_family(fam2members: Mapping[str, Sequence[str]]) -> dict[str, str]:
    """locus tag -> family; a tag listed twice keeps its first family, as Seq2Geno's
    ``gene_clusters2multi_fasta.py`` does."""
    out: dict[str, str] = {}
    for fam, members in fam2members.items():
        for p in members:
            out.setdefault(p, fam)
    return out
