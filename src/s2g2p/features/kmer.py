"""Canonical k-mer counts per genome.

Reproduces what the published pipeline got from KMC (``kmc -k6 -m24 -fm -ci0 -cs1677215`` then
``kmc_dump``): every k-mer of every FASTA record, a k-mer and its reverse complement counted
as one, k-mers containing a non-ACGT base skipped, and no k-mer spanning two records. The
columns are the canonical vocabulary in lexicographic order, the order ``vocab_build`` of the
benchmark produced (the lexicographically smaller of each k-mer / reverse-complement pair).
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

_LUT = np.full(256, 255, dtype=np.uint8)
for _c, _v in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    _LUT[_c] = _v
# KMC's counter ceiling under -cs1677215; reached only by pathological inputs, kept for parity.
KMC_MAX_COUNT = 1677215


@lru_cache(maxsize=None)
def _canonical_map(k: int) -> tuple[np.ndarray, tuple[str, ...]]:
    """(column of every full k-mer index, canonical vocabulary)."""
    n = 4 ** k
    idx = np.arange(n, dtype=np.int64)
    rc = np.zeros(n, dtype=np.int64)
    x = idx.copy()
    for _ in range(k):  # read digits least-significant first, complement, append
        rc = rc * 4 + (3 - x % 4)
        x //= 4
    canon = np.minimum(idx, rc)
    keys = np.unique(canon)
    col = np.searchsorted(keys, canon)
    alphabet = "ACGT"
    vocab = []
    for key in keys:
        s, v = [], int(key)
        for _ in range(k):
            s.append(alphabet[v % 4])
            v //= 4
        vocab.append("".join(reversed(s)))
    return col, tuple(vocab)


def canonical_vocab(k: int) -> list[str]:
    return list(_canonical_map(k)[1])


def _records(path: str | os.PathLike):
    """Each FASTA record's sequence as bytes; headers dropped."""
    with open(path, "rb") as fh:
        chunks: list[bytes] = []
        for line in fh:
            if line.startswith(b">"):
                if chunks:
                    yield b"".join(chunks)
                chunks = []
            else:
                chunks.append(line.strip())
        if chunks:
            yield b"".join(chunks)


def count_canonical(path: str | os.PathLike, k: int = 6) -> np.ndarray:
    """Canonical k-mer counts of one FASTA file, aligned to :func:`canonical_vocab`."""
    col, vocab = _canonical_map(k)
    full = np.zeros(4 ** k, dtype=np.int64)
    for rec in _records(path):
        codes = _LUT[np.frombuffer(rec, dtype=np.uint8)]
        # split into runs of valid bases: a k-mer never spans an N or a record boundary
        bad = np.flatnonzero(codes == 255)
        start = 0
        for stop in list(bad) + [codes.size]:
            run = codes[start:stop].astype(np.int64)
            m = run.size - k + 1
            if m > 0:
                idx = np.zeros(m, dtype=np.int64)
                for j in range(k):
                    idx *= 4
                    idx += run[j:j + m]
                full += np.bincount(idx, minlength=full.size)
            start = stop + 1
    out = np.bincount(col, weights=full, minlength=len(vocab)).astype(np.int64)
    return np.minimum(out, KMC_MAX_COUNT)


def _count_one(args):
    sample, path, k = args
    return sample, count_canonical(path, k)


def kmer_matrix(assemblies: Mapping[str, str], k: int = 6, n_jobs: int = 1) -> pd.DataFrame:
    """samples x canonical k-mers, integer counts, rows in the order given."""
    vocab = canonical_vocab(k)
    jobs = [(s, p, k) for s, p in assemblies.items()]
    if n_jobs == 1:
        rows = dict(map(_count_one, jobs))
    else:
        with ProcessPoolExecutor(n_jobs) as ex:
            rows = dict(ex.map(_count_one, jobs, chunksize=4))
    mat = np.vstack([rows[s] for s in assemblies]) if jobs else np.zeros((0, len(vocab)), dtype=np.int64)
    df = pd.DataFrame(mat, index=list(assemblies), columns=vocab)
    df.index.name = "sample"
    return df


def read_kmc_dump(path: str | os.PathLike, k: int = 6) -> np.ndarray:
    """A ``kmc_dump`` text file as a count vector over :func:`canonical_vocab` (for checking
    this counter against the published one)."""
    vocab = canonical_vocab(k)
    pos = {v: i for i, v in enumerate(vocab)}
    out = np.zeros(len(vocab), dtype=np.int64)
    with open(path) as fh:
        for line in fh:
            kmer, cnt = line.split()
            out[pos[kmer]] = int(cnt)
    return out


def assemblies_from_sheet(sheet: pd.DataFrame) -> dict[str, str]:
    if "assembly" not in sheet.columns:
        raise ValueError("sample sheet has no 'assembly' column")
    missing = [s for s, p in sheet["assembly"].items() if not p or not Path(p).is_file()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} assemblies not found, e.g. {missing[:3]}")
    return dict(sheet["assembly"].items())
