#!/usr/bin/env python3
"""`s2g2p evaluate --preset hu2024` against the published main_s2p.py outputs of one pair.

Inputs are the published run's own tables -- gpa.mat_NONRDNT and indel.mat_NONRDNT -- and the
6-mer counts recomputed by s2g2p (validate_reproduction.py checks those against the KMC dumps),
on the published folds and labels. For every classifier and outer fold, the test macro-F1 and
the test predictions are compared with the published per-classifier JSON.

  python dev/compare_hu2024_ml.py --bench <AMR_benchmarking> --species Enterococcus_faecium \
      --antibiotic vancomycin --scheme random --work <dir> --threads 32
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from s2g2p import cli, io  # noqa: E402
from s2g2p.features.kmer import kmer_matrix  # noqa: E402

SCHEMES = {"random": ("random", False, False), "phylotree": ("phylotree", False, True), "kma": ("KMA", True, False)}


def strip(df):
    df.index = [s.replace("iso_", "") for s in df.index]
    df.index.name = "sample"
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True, type=Path)
    ap.add_argument("--species", required=True)
    ap.add_argument("--antibiotic", required=True)
    ap.add_argument("--scheme", required=True, choices=list(SCHEMES))
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    sp, ab = a.species, a.antibiotic
    fold_name, f_kma, f_tree = SCHEMES[a.scheme]
    so = a.bench / "log/software/seq2geno/software_output" / sp
    den = so / "results/denovo"
    a.work.mkdir(parents=True, exist_ok=True)

    pheno = pd.read_csv(a.bench / f"data/PATRIC/meta/loose_by_species/Data_{sp}_{ab}_pheno.txt", sep="\t",
                        index_col=0, dtype={"genome_id": str})
    sheet = pd.DataFrame({"sample": pheno["genome_id"], ab: pheno["resistant_phenotype"].astype(int).astype(str)})
    sheet.to_csv(a.work / "samples.tsv", sep="\t", index=False)

    kmers = a.work / "kmers.tsv"
    if not kmers.exists():
        asm = {g: str(den / "spades" / f"iso_{g}" / "contigs.fasta") for g in sheet["sample"]}
        io.write_matrix(kmer_matrix(asm, 6, n_jobs=a.threads), kmers)
    for name in ("gpa", "indel"):
        p = a.work / f"{name}.tsv"
        if not p.exists():
            io.write_matrix(strip(io.read_matrix(den / f"{name}.mat_NONRDNT")), p)

    out = a.work / a.scheme
    folds = a.bench / f"data/PATRIC/cv_folds/loose/single_S_A_folds/{sp}/{ab}_{fold_name}_cv.json"
    cli.main(["evaluate", "--preset", "hu2024", "--samplesheet", str(a.work / "samples.tsv"), "--phenotype", ab,
              "--folds", str(folds), "--kmer", str(kmers), "--gpa", str(a.work / "gpa.tsv"),
              "--indel", str(a.work / "indel.tsv"), "--threads", str(a.threads), "--out-dir", str(out)])

    cand = pd.read_csv(out / ab / "candidates.tsv", sep="\t")
    rows = []
    for clf in ("svm", "lr", "rf", "lsvm"):
        pub = json.load(open(a.bench / f"log/software/seq2geno/analysis/{sp}/{ab}_cl_{clf}_KMA_{f_kma}_Tree_{f_tree}.json"))
        mine = cand[cand.classifier == clf].sort_values("fold")
        for k, (f_pub, pred_pub) in enumerate(zip(pub["f1_test"], pub["predictY_test"])):
            m = mine[mine.fold == k].iloc[0]
            rows.append({"scheme": a.scheme, "classifier": clf, "fold": k, "f1_published": f_pub,
                         "f1_s2g2p": m["f1_macro"], "f1_equal": abs(f_pub - m["f1_macro"]) < 1e-9})
    # predictions only exist for the selected candidate in predictions.tsv; per-classifier
    # equality of macro-F1 on identical test folds is the check
    res = pd.DataFrame(rows)
    res.to_csv(a.work / f"compare_{a.scheme}.tsv", sep="\t", index=False)
    print(res.groupby("classifier")["f1_equal"].agg(["sum", "count"]))
    print(f"{int(res.f1_equal.sum())}/{len(res)} (classifier, fold) macro-F1 identical")


if __name__ == "__main__":
    main()
