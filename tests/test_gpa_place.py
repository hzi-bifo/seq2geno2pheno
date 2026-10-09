import pandas as pd

from s2g2p.features.gpa import gpa_from_roary, read_clustered_proteins
from s2g2p.place import best_hits, place_gpa, training_sequence_families

from test_indel import roary_csv


def test_gpa_from_roary(tmp_path):
    roary_csv(tmp_path / "g.csv", [("famA", 2, {"s1": "s1_1", "s2": "s2_1"}), ("famB", 1, {"s2": "s2_2"})])
    df = gpa_from_roary(tmp_path / "g.csv")
    assert list(df.index) == ["s1", "s2"]
    assert list(df.columns) == ["famA", "famB"]
    assert df.to_numpy().tolist() == [[1, 0], [1, 1]]
    df2 = gpa_from_roary(tmp_path / "g.csv", samples=["s2"], rename={"famA": "x,y,famA"})
    assert list(df2.columns) == ["x,y,famA", "famB"]


def test_best_hits_uses_roary_criterion(tmp_path):
    cols = "q s pident length mm go qs qe ss se evalue bitscore slen qlen".split()
    rows = [
        ["q0", "s0", 99, 100, 0, 0, 1, 100, 1, 100, 0, 200, 100, 100],   # ok
        ["q0", "s1", 99, 100, 0, 0, 1, 100, 1, 100, 0, 300, 100, 100],   # ok, higher bitscore
        ["q1", "s0", 95, 100, 0, 0, 1, 100, 1, 100, 0, 500, 100, 100],   # identity not > 95
        ["q2", "s0", 99, 95, 0, 0, 1, 95, 1, 95, 0, 500, 100, 95],        # covers 95% of subject only
    ]
    f = tmp_path / "h.tsv"
    f.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
    assert best_hits(f) == {"q0": 1}


def test_place_exact_matches_need_no_blast(tmp_path):
    (tmp_path / "cp").write_text("famA: t1_1\nfamB: t1_2\n")
    (tmp_path / "t1.faa").write_text(">t1_1 x\nMKV\n>t1_2 y\nMAA\n")
    (tmp_path / "q1.faa").write_text(">q1_1 x\nMAA\n")
    seq2fam, stats = training_sequence_families(tmp_path / "cp", {"t1": tmp_path / "t1.faa"})
    assert seq2fam == {"MKV": "famA", "MAA": "famB"}
    table, assign, st = place_gpa(seq2fam, {"q1": tmp_path / "q1.faa"}, ["famA", "famB"],
                                  blastp="/nonexistent", makeblastdb="/nonexistent")
    assert table.loc["q1"].tolist() == [0, 1]
    assert st["query_blasted"] == 0
    assert assign.to_dict("records") == [{"sample": "q1", "locus": "q1_1", "family": "famB"}]
    assert list(read_clustered_proteins(tmp_path / "cp")) == ["famA", "famB"]
