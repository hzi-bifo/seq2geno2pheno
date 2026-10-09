import numpy as np
import pandas as pd
import pytest

from s2g2p.features.compress import DedupColumns, compress


def feat_compress_reference(df):
    """Seq2Geno denovo/featCompress.py, minus snakemake I/O."""
    t_df = df.transpose()
    grouped = t_df.groupby(t_df.columns.values.tolist())
    new_feat, groups = {}, {}
    for name, group in grouped:
        repre = str(group.index.tolist()[0])
        new_repre = "|".join([repre, str(len(group))])
        new_feat[new_repre] = t_df.loc[repre, :]
        groups[new_repre] = group.index.tolist()
    return pd.DataFrame(data=new_feat), groups


@pytest.mark.parametrize("seed", range(5))
def test_matches_featcompress(seed):
    rng = np.random.default_rng(seed)
    n, m = 12, 60
    base = rng.integers(0, 2, size=(n, 8))
    cols = base[:, rng.integers(0, 8, size=m)]  # many duplicated patterns
    vals = np.where(cols == 1, "1", "0").astype(object)
    vals[rng.random(vals.shape) < 0.05] = "NA"
    df = pd.DataFrame(vals, index=[f"s{i}" for i in range(n)], columns=[f"f{j}" for j in range(m)])
    ref, ref_groups = feat_compress_reference(df)
    out, groups = compress(df)
    assert list(out.columns) == list(ref.columns)
    assert (out.to_numpy() == ref.to_numpy()).all()
    assert groups == ref_groups


def test_empty_cell_is_an_error():
    df = pd.DataFrame({"a": ["1", None], "b": ["0", "1"]}, index=["x", "y"])
    with pytest.raises(ValueError):
        compress(df)


def test_dedup_transformer_keeps_compress_columns():
    rng = np.random.default_rng(3)
    x = rng.integers(0, 2, size=(20, 4))[:, rng.integers(0, 4, size=30)]
    df = pd.DataFrame(x.astype(str), columns=[f"f{j}" for j in range(30)])
    out, _ = compress(df)
    d = DedupColumns().fit(x)
    assert [f"f{j}" for j in d.keep_] == [c.split("|")[0] for c in out.columns]
    assert np.array_equal(d.transform(x), out.to_numpy().astype(int))
