import pandas as pd

from s2g2p import io


def test_numeric_looking_ids_survive(tmp_path):
    df = pd.DataFrame({"f": [1, 2]}, index=["1352.11330", "1352.1133"])
    io.write_matrix(df, tmp_path / "m.tsv")
    for dtype in (None, float, int):
        back = io.read_matrix(tmp_path / "m.tsv", dtype=dtype)
        assert list(back.index) == ["1352.11330", "1352.1133"]


def test_folds_round_trip(tmp_path):
    folds = [["a", "b"], ["c"]]
    io.write_folds(folds, tmp_path / "f.tsv")
    assert io.read_folds(tmp_path / "f.tsv") == folds
