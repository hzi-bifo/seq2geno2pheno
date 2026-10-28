import collections
import copy
import json

import numpy as np
import pandas as pd
import pytest
from sklearn import utils

from s2g2p import cli, folds as F


# ------------------------------------------------- published code (AMR_benchmarking, MIT), near-verbatim
def cluster_split_reference(dict_cluster, Random_State, cv):
    """AMR_benchmarking src/cv_folds/cluster2folds.py cluster_split ("khu: modified"), prints removed."""
    all_data_splits = []
    all_samples = []
    for i in dict_cluster:
        for each in dict_cluster[i]:
            all_samples.append(each)
    shuffled = list(utils.shuffle(list(dict_cluster.keys()), random_state=Random_State))
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
            elements = dict_cluster[str(e)]
            tem.append(len(elements))
        totals.append(sum(tem))
    for i in range(len(all_data_splits)):
        extracted = all_data_splits[i]
        sum_sub = totals[i]
        b = 0
        while (sum_sub < 0.2 * (len(all_samples) / float(cv))) and b < 100:
            b += 1
            m_from = np.argmax(totals)
            extra = list(utils.shuffle(all_data_splits[m_from], random_state=Random_State))[0]
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
            tem_Nsamples = []
            for e in sorted(list(set(all_data_splits[m_from]) - set([extra]))):
                elements = dict_cluster[str(e)]
                tem_Nsamples.append(len(elements))
            sum_from = sum(tem_Nsamples)
            extracted = extracted + [extra]
            tem_Nsamples = []
            for e in extracted:
                elements = dict_cluster[str(e)]
                tem_Nsamples.append(len(elements))
            sum_sub = sum(tem_Nsamples)
            if len(dict_cluster[extra]) < 1.0 * (len(all_samples) / float(cv)) and sum_from > sum_sub:
                totals[i] = sum_sub
                totals[m_from] = sum_from
                all_data_splits[i] = extracted
                all_data_splits[m_from] = sorted(list(set(all_data_splits[m_from]) - set([extra])))
    return all_data_splits


def prepare_cluster_reference(cluster_summary):
    """cluster2folds.prepare_cluster, single-species branch, on the file's lines."""
    dict_cluster = collections.defaultdict(list)
    for each_cluster in cluster_summary:
        splitted = each_cluster.split()
        if splitted != []:
            if splitted[0].isdigit():
                dict_cluster[str(int(splitted[0]) - 1)].append(splitted[3])
            if splitted[0] == "Similar":
                splitted = each_cluster.split()
                splitted_2 = each_cluster.split(":")
                dict_cluster[str(int(splitted_2[1].split()[0]) - 1)].append(splitted[6])
    return dict_cluster


def prepare_folders_reference(cv, Random_State, names, cluster_summary):
    """cluster2folds.prepare_folders(..., 'new'), single species, then main_s2p's index -> name."""
    dict_cluster = prepare_cluster_reference(cluster_summary)
    all_data_splits = cluster_split_reference(dict_cluster, Random_State, cv)
    folders_sample = []
    for out_cv in range(cv):
        folders_sample_sub = []
        for cl_ID in all_data_splits[out_cv]:
            for element in dict_cluster[cl_ID]:
                if element in names:
                    folders_sample_sub.append(names.index(element))
        folders_sample.append(folders_sample_sub)
    return [[names[s] for s in f] for f in folders_sample]


def merge_scaffolds_reference(input_path, files, output, spacer):
    """Aytan-Aktug data_preparation/merge_scaffolds.extract_info (ids given as a list)."""
    new_file = open(output, "w")
    for i in files:
        open_scaf = open("%s/%s.fna" % (input_path, i), "r")
        scaf = open_scaf.readlines()
        new_file.write(">%s" % i)
        new_file.write("\n")
        for each in scaf:
            if ">" not in each:
                each_new = each.replace("\n", "")
                new_file.write(each_new)
            else:
                new_file.write(str("N" * spacer))
    new_file.close()


# ------------------------------------------------------------------------------------- helpers
def kma_log(sizes, rng, names):
    """A kma_clustering log: one '# Added' line per cluster, then 'Similar to q:' lines for the
    other members, interleaved as the program prints them, with noise lines in between."""
    lines = ["# Reading inputfile: \t--", "# Preparing compressed DB."]
    k = 0
    pending = []
    for c, size in enumerate(sizes, start=1):
        lines.append(f"{c}\t# Added:\t{names[k]}")
        k += 1
        for _ in range(size - 1):
            pending.append((c, names[k]))
            k += 1
        while pending and rng.random() < 0.5:
            q, s = pending.pop(rng.integers(len(pending)))
            lines.append(f"Similar to q:{q} t:{q}\t# Skipped:\t{s}")
    for q, s in pending:
        lines.append(f"Similar to q:{q} t:{q}\t# Skipped:\t{s}")
    lines.append("# Templates key-value pairs:\t1234.")
    return [l + "\n" for l in lines]


@pytest.mark.parametrize("seed", range(8))
def test_cluster_split_matches_published(seed):
    rng = np.random.default_rng(seed)
    n_clusters = int(rng.integers(12, 80))
    sizes = rng.geometric(0.15, size=n_clusters)  # many small clusters, a few big ones
    if seed % 2:
        sizes[0] = sizes.sum()  # one cluster larger than a fold
    d = collections.OrderedDict((str(i), [f"g{i}_{j}" for j in range(s)]) for i, s in enumerate(sizes))
    for cv in (5, 10):
        assert F.cluster_split(d, 42, cv) == cluster_split_reference(d, 42, cv)


@pytest.mark.parametrize("seed", range(5))
def test_kma_folds_match_published(tmp_path, seed):
    rng = np.random.default_rng(seed)
    sizes = rng.geometric(0.2, size=int(rng.integers(15, 60)))
    names = [f"1352.{10000 + i}" for i in range(int(sizes.sum()))]
    log = kma_log(sizes, rng, names)
    (tmp_path / "c.txt").write_text("".join(log))
    labelled = [s for s in names if rng.random() < 0.9]  # the log may name more genomes
    clusters = F.read_kma_clusters(tmp_path / "c.txt")
    assert dict(clusters) == dict(prepare_cluster_reference(log))
    assert list(clusters) == list(prepare_cluster_reference(log))
    got = F.kma_folds(clusters, labelled, 10, 42)
    want = prepare_folders_reference(10, 42, labelled, log)
    assert got == want
    F.validate(got, labelled)


def test_kma_unclustered_genome_gets_its_own_cluster(tmp_path):
    log = ["1\t# Added:\ta\n", "Similar to q:1 t:1\t# Skipped:\tb\n", "2\t# Added:\tc\n",
           "Too short\t# Skipped:\td\n", "3\t# Added:\te\n", "4\t# Added:\tf\n"]
    (tmp_path / "c.txt").write_text("".join(log))
    clusters = F.read_kma_clusters(tmp_path / "c.txt")
    samples = list("abcdef")
    assert F.unclustered(clusters, samples) == ["d"]
    folds = F.kma_folds(clusters, samples, 2)
    F.validate(folds, samples)
    where = {s: k for k, f in enumerate(folds) for s in f}
    assert where["a"] == where["b"]  # a cluster stays in one fold


def test_kma_input_matches_merge_scaffolds(tmp_path):
    (tmp_path / "asm").mkdir()
    ids = ["s1", "s2", "s3"]
    for i, s in enumerate(ids):
        contigs = "".join(f">{s}_c{c} len=9\nACGTN\nTTGCA\n" for c in range(i + 1))
        (tmp_path / "asm" / f"{s}.fna").write_text(contigs)
    merge_scaffolds_reference(tmp_path / "asm", ids, tmp_path / "ref.fna", 16)
    F.write_kma_input([(s, tmp_path / "asm" / f"{s}.fna") for s in ids], tmp_path / "ours.fna")
    assert (tmp_path / "ours.fna").read_bytes() == (tmp_path / "ref.fna").read_bytes()


# ------------------------------------------------------------------------------------ phylogeny
def test_newick_leaves():
    t = "((a:0.1,'b c':0.2)95:0.3,(d,e)[&comment]:1e-3,'f''g':2)root;"
    assert F.newick_leaves(t) == ["a", "b c", "d", "e", "f'g"]
    ape = "(1352.10474:0.0012,(1352.10148:3.1e-05,1352.10167:-2e-06):0.004,1352.11265:0.1);\n"
    assert F.newick_leaves(ape) == ["1352.10474", "1352.10148", "1352.10167", "1352.11265"]
    with pytest.raises(ValueError, match="duplicated"):
        F.newick_leaves("(a,(b,a));")


def test_tree_folds_are_equal_leaf_order_blocks():
    order = [f"s{i}" for i in range(300)]
    labelled = order[3:280]  # 277, like E. faecium vancomycin
    folds = F.tree_folds(order[::-1], labelled, 10)
    assert [len(f) for f in folds] == [28] * 7 + [27] * 3
    assert sum(folds, []) == [s for s in order[::-1] if s in set(labelled)]
    with pytest.raises(ValueError, match="not leaves"):
        F.tree_folds(order[:100], labelled, 10)


# ------------------------------------------------------------------------------------------ cli
def sheet(tmp_path, n=40):
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"sample": [f"g{i}" for i in range(n)],
                       "assembly": [f"asm/g{i}.fna" for i in range(n)],
                       "drugA": rng.choice(["R", "S"], size=n), "drugB": rng.choice(["1", "0", ""], size=n)})
    df.loc[:1, "drugA"] = ["R", "S"]
    df.loc[:3, "drugB"] = ["1", "1", "0", "0"]
    (tmp_path / "asm").mkdir()
    for s in df["sample"]:
        (tmp_path / "asm" / f"{s}.fna").write_text(f">{s}_1\nACGT\n>{s}_2\nGGCC\n")
    df.to_csv(tmp_path / "samples.tsv", sep="\t", index=False)
    return df


@pytest.mark.parametrize("method", ["random", "phylogeny", "kma"])
def test_cli_folds(tmp_path, method):
    df = sheet(tmp_path)
    args = ["folds", "--samplesheet", str(tmp_path / "samples.tsv"), "--method", method,
            "--n-folds", "3", "--out-dir", str(tmp_path / "folds")]
    if method == "phylogeny":
        (tmp_path / "t.nwk").write_text("(" + ",".join(f"{s}:1" for s in df["sample"][::-1]) + ");")
        args += ["--tree", str(tmp_path / "t.nwk")]
    if method == "kma":
        for ph in ("drugA", "drugB"):
            ids = [s for s, v in zip(df["sample"], df[ph]) if v != ""]
            log = [f"{i + 1}\t# Added:\t{s}\n" for i, s in enumerate(ids)]
            (tmp_path / f"{ph}.clusters.txt").write_text("".join(log))
        args += ["--kma-clusters", str(tmp_path / "{phenotype}.clusters.txt")]
    cli.main(args)
    man = pd.read_csv(tmp_path / "folds" / "manifest.tsv", sep="\t")
    assert list(man["file"]) == [f"drugA_{method}_cv.json", f"drugB_{method}_cv.json"]
    for ph in ("drugA", "drugB"):
        got = json.loads((tmp_path / "folds" / f"{ph}_{method}_cv.json").read_text())
        labelled = [s for s, v in zip(df["sample"], df[ph]) if v != ""]
        F.validate(got, labelled)
        if method == "phylogeny":
            assert sum(got, []) == [s for s in df["sample"][::-1] if s in set(labelled)]
    # the files are what --folds reads
    cli.main(["split", "--samplesheet", str(tmp_path / "samples.tsv"),
              "--folds", str(tmp_path / "folds" / f"{{phenotype}}_{method}_cv.json"), "--out-dir", str(tmp_path / "sp")])


def test_cli_kma_input(tmp_path):
    df = sheet(tmp_path, n=12)
    cli.main(["kma-input", "--samplesheet", str(tmp_path / "samples.tsv"), "--phenotype", "drugB",
              "--out", str(tmp_path / "m.fna")])
    ids = [s for s, v in zip(df["sample"], df["drugB"]) if v != ""]
    merge_scaffolds_reference(tmp_path / "asm", ids, tmp_path / "ref.fna", 16)
    assert (tmp_path / "m.fna").read_bytes() == (tmp_path / "ref.fna").read_bytes()
