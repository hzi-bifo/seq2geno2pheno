# Validation against the published pipeline

The reference is the Hu et al. 2024 reproduction in the amr_benchmark tree: Seq2Geno run as
published on PATRIC assemblies, then `main_s2p.py` and `result_analysis.py`. Every check compares
v3's code with that run's own files. No check compares v3 with itself.

## 1. Feature ports (`dev/validate_reproduction.py`)

Each published intermediate is recomputed from the same run's upstream files.

| Species | Genomes | Checks identical |
|---|---:|---:|
| *Enterococcus faecium* | 277 | 7/7 |
| *Campylobacter jejuni* | 395 | 7/7 |
| *Neisseria gonorrhoeae* | 680 | 7/7 |
| *Acinetobacter baumannii* | 1,144 | 7/7 |
| *Salmonella enterica* | 1,922 | 7/7 |
| *Staphylococcus aureus* | 3,325 | 7/7 |
| *Streptococcus pneumoniae* | 4,757 | 7/7 |

The seven checks are:

1. `gpa.mat`;
2. featCompress of `gpa.mat` and its groups;
3. featCompress of `indel.mat` and its groups;
4. the set of aligned families;
5. every `vcf2indel` call;
6. `indel.mat`;
7. canonical 6-mer counts vs the KMC dumps (all 277 genomes for *E. faecium*, 100 sampled
   elsewhere).

Per-species detail is in `<species>.tsv` here.

*S. aureus* is the species that found the published indel rule's quirks (`docs/design.md`).
On *E. faecium* alone, an incomplete port passed.

## 2. Models (`dev/compare_hu2024_ml.py`)

The inputs are the published `gpa.mat_NONRDNT` and `indel.mat_NONRDNT`, s2g2p's 6-mers
(identical to KMC, section 1), and the published folds and labels. `s2g2p evaluate --preset
hu2024` is then compared with the published per-classifier JSON and with the reproduction's
final score (`result/hu2024-repro/reproduced_f1_macro.tsv`).

*E. faecium* vancomycin:

| Folds | Identical, scikit-learn 1.6.1 | Identical, scikit-learn 0.24.1 | Final macro-F1, s2g2p | Final macro-F1, published code |
|---|---:|---:|---:|---:|
| random | 37/40 | **40/40** | 0.9737736929124945 | 0.9737736929124945 |
| phylotree | 35/40 | **40/40** | 0.9892442396313365 | 0.9892442396313365 |
| KMA (support-weighted) | 39/40 | **40/40** | 0.947718517674583 | 0.947718517674583 |

The scikit-learn 0.24.1 column ran in the published model environment (`amr_env`: Python 3.7,
numpy 1.20.1, pandas 1.2.2; `containers/hu2024ml/`). There, **every classifier on every fold
matches, 120 of 120**.

**The final score is identical in all three fold schemes, to every printed digit.**

The per-classifier differences (3, 5 and 1 of 40) come from scikit-learn 0.24.1 in the
published `amr_env` vs 1.6.1 here. Logistic regression and random forest both changed
internally between those versions. Linear SVM is selected in all 30 folds and matches
everywhere, so none of the differences reaches the reported score.

That is a property of this pair, not a guarantee. On a pair where random forest or logistic
regression wins or nearly wins, a version-induced shift in its inner score can change the
selection. In the published environment that risk is gone (120/120 above), which is why the
`hu2024` preset runs its models there.

Per-fold detail: `dev-runs/hu2024-ml/*/compare_<scheme>.tsv` (not committed; regenerate with
the script).
