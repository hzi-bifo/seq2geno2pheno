import pandas as pd

from s2g2p import cli, io

from test_evaluate import synthetic, small_grids  # noqa: F401  (fixture)


def write_inputs(tmp_path, labels, folds, blocks):
    sheet = pd.DataFrame({"sample": labels.index, "drug": labels.map({0: "S", 1: "R"}).to_numpy()})
    sheet.loc[len(sheet)] = ["unlabelled1", ""]  # in the cohort, not in the model
    sheet.to_csv(tmp_path / "samples.tsv", sep="\t", index=False)
    io.write_folds(folds, tmp_path / "folds.tsv")
    for name, df in blocks.items():
        for k in range(len(folds)):
            (tmp_path / f"f{k}").mkdir(exist_ok=True)
            io.write_matrix(df, tmp_path / f"f{k}" / f"{name}.tsv")
    io.write_matrix(blocks["kmer"], tmp_path / "kmer.tsv")


def test_evaluate_train_predict(tmp_path, small_grids):
    labels, folds, blocks = synthetic(3)
    write_inputs(tmp_path, labels, folds, blocks)
    common = ["--samplesheet", str(tmp_path / "samples.tsv"), "--folds", str(tmp_path / "folds.tsv"),
              "--kmer", str(tmp_path / "kmer.tsv"), "--gpa", str(tmp_path / "f{fold}" / "gpa.tsv"),
              "--indel", str(tmp_path / "f{fold}" / "indel.tsv")]
    cli.main(["evaluate", "--preset", "noleak", "--out-dir", str(tmp_path / "out"), *common])
    out = tmp_path / "out" / "drug"
    man = pd.read_csv(out / "manifest.tsv", sep="\t")
    assert set(man["file"]) >= {"summary.tsv", "predictions.tsv", "folds.tsv", "candidates.tsv", "run.json"}
    pred = pd.read_csv(out / "predictions.tsv", sep="\t")
    assert len(pred) == len(labels) and "unlabelled1" not in set(pred["sample"])

    flat = [a.replace("f{fold}", "f0") for a in common]
    cli.main(["train", "--preset", "noleak", "--phenotype", "drug", "--out", str(tmp_path / "m.joblib"), *flat])
    cli.main(["predict", "--model", str(tmp_path / "m.joblib"), "--out", str(tmp_path / "p.tsv"),
              "--kmer", str(tmp_path / "kmer.tsv"), "--gpa", str(tmp_path / "f0" / "gpa.tsv"),
              "--indel", str(tmp_path / "f0" / "indel.tsv")])
    assert len(pd.read_csv(tmp_path / "p.tsv", sep="\t")) == len(labels)


def test_presets_listing(capsys):
    cli.main(["presets"])
    out = capsys.readouterr().out
    for name in ("hu2024", "noleak", "fast", "g2p"):
        assert name in out
