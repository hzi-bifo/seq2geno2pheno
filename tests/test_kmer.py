import itertools
import random

import numpy as np

from s2g2p.features.kmer import canonical_vocab, count_canonical, kmer_matrix

COMP = str.maketrans("ACGT", "TGCA")


def ref_vocab(k):
    """The benchmark's vocab_build(canonical=True): lexicographic, smaller of each pair."""
    vocab = ["".join(x) for x in itertools.product("ACGT", repeat=k)]
    keep = []
    for i, v in enumerate(vocab):
        rc = v.translate(COMP)[::-1]
        if rc not in vocab[:i]:
            keep.append(v)
    return keep


def ref_count(records, k):
    vocab = ref_vocab(k)
    pos = {v: i for i, v in enumerate(vocab)}
    out = np.zeros(len(vocab), dtype=np.int64)
    for seq in records:
        seq = seq.upper()
        for i in range(len(seq) - k + 1):
            w = seq[i:i + k]
            if set(w) <= set("ACGT"):
                rc = w.translate(COMP)[::-1]
                out[pos[min(w, rc)]] += 1
    return out


def test_vocab_matches_vocab_build():
    for k in (1, 2, 3, 4, 6):
        assert canonical_vocab(k) == ref_vocab(k)
    assert len(canonical_vocab(6)) == 2080


def test_counts_match_brute_force(tmp_path):
    rng = random.Random(7)
    records = []
    for _ in range(5):
        s = "".join(rng.choice("ACGTacgtN") if rng.random() < 0.1 else rng.choice("ACGT")
                    for _ in range(rng.randint(0, 400)))
        records.append(s)
    f = tmp_path / "g.fna"
    with open(f, "w") as fh:
        for i, s in enumerate(records):
            fh.write(f">c{i} desc\n")
            for j in range(0, len(s), 60):  # wrapped lines must not break k-mers
                fh.write(s[j:j + 60] + "\n")
    for k in (2, 3, 6):
        assert np.array_equal(count_canonical(f, k), ref_count(records, k))


def test_no_kmer_spans_records(tmp_path):
    f = tmp_path / "g.fna"
    f.write_text(">a\nAAA\n>b\nAAA\n")
    v = dict(zip(canonical_vocab(4), count_canonical(f, 4)))
    assert v["AAAA"] == 0


def test_matrix_rows_follow_input(tmp_path):
    a, b = tmp_path / "a.fa", tmp_path / "b.fa"
    a.write_text(">x\nACGTACGT\n")
    b.write_text(">x\nTTTTTT\n")
    df = kmer_matrix({"b": str(b), "a": str(a)}, k=3)
    assert list(df.index) == ["b", "a"]
    assert df.loc["b", "AAA"] == 4  # TTT folds onto AAA
