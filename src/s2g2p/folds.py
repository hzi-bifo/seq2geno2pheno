"""Outer cross-validation folds.

The published evaluation (Hu et al. 2024) used three fold sets per (species, antibiotic), and
every preset accepts such a file unchanged. For a cohort without one, ``s2g2p folds`` makes them
with the same three methods:

``random``     seeded stratified k-fold over the labelled genomes (:func:`stratified_folds`).
               The published random folds were made by the Geno2Pheno server, also stratified,
               with an ordering and seed that are not recorded; they cannot be regenerated.
``phylogeny``  the leaf order of a neighbour-joining tree of Roary's core-gene alignment,
               cut into k equal blocks (:func:`tree_folds`). The published phylogeny-aware
               folds were made by the Geno2Pheno server from such a tree; its code was never
               released. The supplement says it "generates folds of equal size ... based on the
               clade structure", and the published fold sizes are exactly ``np.array_split``'s
               and not stratified. On E. faecium vancomycin, the published folds are as
               tree-shaped on our own tree as equal leaf-order blocks are (docs/folds.md).
``kma``        Aytan-Aktug's homology-aware folds as Hu et al. refined them: kma_clustering of
               the genomes, then clusters allocated to folds (:func:`cluster_split`, a port of
               AMR_benchmarking src/cv_folds/cluster2folds.py, MIT).
"""

from __future__ import annotations

import collections
import copy
import os
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

METHODS = ("random", "phylogeny", "kma")
# Published fold allocation: cluster2folds.prepare_folders(cv, Random_State=42, ..., 'new')
KMA_SEED = 42
# merge_scaffolds.extract_info(..., 16): Ns between contigs, longer than the k-mer size
KMA_SPACER = 16


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


# ------------------------------------------------------------------------------------ phylogeny
def newick_leaves(text: str) -> list[str]:
    """Leaf names of a Newick tree, in the order they are written.

    A label right after ``(`` or ``,`` is a leaf; one after ``)`` names an internal node (or
    carries a support value) and is skipped. Branch lengths, [comments] and 'quoted' labels are
    handled.
    """
    leaves: list[str] = []
    i, n = 0, len(text)
    prev = "("
    while i < n:
        c = text[i]
        if c in "(),;":
            prev = c
            i += 1
        elif c == "[":
            i = text.index("]", i) + 1
        elif c == ":":
            i += 1
            while i < n and text[i] not in "(),;[":
                i += 1
        elif c.isspace():
            i += 1
        else:
            if c == "'":
                j = i + 1
                label = []
                while True:
                    j2 = text.index("'", j)
                    label.append(text[j:j2])
                    if j2 + 1 < n and text[j2 + 1] == "'":   # '' is an escaped quote
                        label.append("'")
                        j = j2 + 2
                    else:
                        i = j2 + 1
                        break
                name = "".join(label)
            else:
                j = i
                while j < n and text[j] not in "(),;:[" and not text[j].isspace():
                    j += 1
                name = text[i:j]
                i = j
            if prev in "(,":
                leaves.append(name)
            prev = "label"
    if len(set(leaves)) != len(leaves):
        dup = [s for s, c in collections.Counter(leaves).items() if c > 1]
        raise ValueError(f"tree has duplicated leaf names, e.g. {dup[:3]}")
    return leaves


def tree_folds(leaf_order: Sequence[str], samples: Iterable[str], n_folds: int = 10) -> list[list[str]]:
    """The labelled genomes in tree leaf order, cut into ``n_folds`` contiguous blocks of equal
    size (the first ``n % n_folds`` blocks one larger, as ``np.array_split``)."""
    want = set(samples)
    missing = want - set(leaf_order)
    if missing:
        raise ValueError(f"{len(missing)} labelled genomes are not leaves of the tree, e.g. {sorted(missing)[:3]}")
    order = [s for s in leaf_order if s in want]
    if len(order) < n_folds:
        raise ValueError(f"{len(order)} genomes cannot make {n_folds} folds")
    return [b.tolist() for b in np.array_split(np.array(order, dtype=object), n_folds)]


# ------------------------------------------------------------------------------------------ kma
def write_kma_input(assemblies: Sequence[tuple[str, str]], out: str | os.PathLike,
                    spacer: int = KMA_SPACER) -> None:
    """One record per genome, its contigs joined by ``spacer`` Ns: Aytan-Aktug's
    merge_scaffolds.extract_info, byte for byte (an N run also precedes the first contig, and a
    record has no final newline, so the next header follows on the same line -- kma_clustering
    reads it as a new record)."""
    with open(out, "w") as new_file:
        for sample, path in assemblies:
            new_file.write(">%s" % sample)
            new_file.write("\n")
            with open(path) as fh:
                for each in fh:
                    if ">" not in each:
                        new_file.write(each.replace("\n", ""))
                    else:
                        new_file.write("N" * spacer)


def read_kma_clusters(path: str | os.PathLike) -> dict[str, list[str]]:
    """kma_clustering's log as {cluster id: genomes}, in the order the log names them
    (cluster2folds.prepare_cluster, single-species branch). Ids are 0-based strings."""
    dict_cluster: dict[str, list[str]] = collections.defaultdict(list)
    with open(path) as fh:
        for each_cluster in fh:
            splitted = each_cluster.split()
            if splitted != []:
                if splitted[0].isdigit():
                    dict_cluster[str(int(splitted[0]) - 1)].append(splitted[3])
                if splitted[0] == "Similar":
                    splitted_2 = each_cluster.split(":")
                    dict_cluster[str(int(splitted_2[1].split()[0]) - 1)].append(splitted[6])
    return dict_cluster


def cluster_split(dict_cluster: Mapping[str, list[str]], random_state: int, cv: int) -> list[list[str]]:
    """Clusters into ``cv`` folds: shuffled, cut into equal runs of clusters, then clusters moved
    from the largest fold into any fold below 20% of the expected size.

    cluster2folds.cluster_split ("khu: modified"), unchanged but for its prints.
    """
    from sklearn import utils

    all_data_splits = []
    all_samples = []
    for i in dict_cluster:
        for each in dict_cluster[i]:
            all_samples.append(each)

    shuffled = list(utils.shuffle(list(dict_cluster.keys()), random_state=random_state))
    r = int(len(shuffled) / cv)
    a = 0
    b = r
    for i in range(cv):
        all_data_splits.append(shuffled[a:b])
        if i != cv - 2:
            a = b
            b = b + r
        else:
            a = b
            b = len(shuffled)

    totals = []
    for folder_cluster in all_data_splits:
        tem = []
        for e in folder_cluster:
            tem.append(len(dict_cluster[str(e)]))
        totals.append(sum(tem))
    for i in range(len(all_data_splits)):
        extracted = all_data_splits[i]
        sum_sub = totals[i]
        b = 0
        while (sum_sub < 0.2 * (len(all_samples) / float(cv))) and b < 100:
            b += 1
            m_from = np.argmax(totals)
            extra = list(utils.shuffle(all_data_splits[m_from], random_state=random_state))[0]
            a = 0
            while len(dict_cluster[extra]) >= 1.0 * (len(all_samples) / float(cv)) and a < 5:
                m_from = np.argmax(totals)
                extra = list(utils.shuffle(all_data_splits[m_from], random_state=a))[0]
                a += 1
            a = 0
            totals_sub = copy.deepcopy(totals)
            while len(dict_cluster[extra]) >= 1.0 * (len(all_samples) / float(cv)) and a < 5:
                totals_sub[m_from] = 0
                m_from = np.argmax(totals_sub)
                extra = list(utils.shuffle(all_data_splits[m_from], random_state=a))[0]
                a += 1
            # make sure the fold giving out the cluster still has enough samples left
            tem_Nsamples = []
            for e in sorted(list(set(all_data_splits[m_from]) - set([extra]))):
                tem_Nsamples.append(len(dict_cluster[str(e)]))
            sum_from = sum(tem_Nsamples)
            extracted = extracted + [extra]
            tem_Nsamples = []
            for e in extracted:
                tem_Nsamples.append(len(dict_cluster[str(e)]))
            sum_sub = sum(tem_Nsamples)
            if len(dict_cluster[extra]) < 1.0 * (len(all_samples) / float(cv)) and sum_from > sum_sub:
                totals[i] = sum_sub
                totals[m_from] = sum_from
                all_data_splits[i] = extracted
                all_data_splits[m_from] = sorted(list(set(all_data_splits[m_from]) - set([extra])))
    return all_data_splits


def kma_folds(dict_cluster: Mapping[str, list[str]], samples: Sequence[str], n_folds: int = 10,
              seed: int = KMA_SEED) -> list[list[str]]:
    """Genomes per fold from kma clusters (cluster2folds.prepare_folders, f_version 'new'):
    members are listed cluster by cluster, and only labelled genomes are kept.

    A labelled genome that kma_clustering did not report ("Too short") would be in no fold (the
    published code dropped it silently); each such genome is added as a cluster of its own,
    after the reported ones, and :func:`unclustered` lists them for the manifest."""
    dict_cluster = collections.OrderedDict((k, list(v)) for k, v in dict_cluster.items())
    seen = {s for v in dict_cluster.values() for s in v}
    for s in samples:
        if s not in seen:
            dict_cluster[str(len(dict_cluster))] = [s]
    names = set(samples)
    splits = cluster_split(dict_cluster, seed, n_folds)
    return [[e for cl in splits[k] for e in dict_cluster[cl] if e in names] for k in range(n_folds)]


def unclustered(dict_cluster: Mapping[str, list[str]], samples: Sequence[str]) -> list[str]:
    seen = {s for v in dict_cluster.values() for s in v}
    return [s for s in samples if s not in seen]
