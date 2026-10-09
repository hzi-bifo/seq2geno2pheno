<!--
SPDX-FileCopyrightText: 2021 Ehsaneddin Asgari and Tzu-Hao Kuo
SPDX-FileCopyrightText: 2026 Hadi Foroughmand Araabi

SPDX-License-Identifier: GPL-3.0-or-later
-->

# Seq2Geno2Pheno

Predict a bacterial phenotype, such as antimicrobial resistance, from genome assemblies, and
measure how well that prediction works.

- **Seq2Geno** turns assemblies into genotype tables: gene presence/absence (Prokka + Roary),
  indels in gene families, and canonical 6-mer counts.
- **Geno2Pheno** fits and evaluates classifiers on those tables: SVM, logistic regression,
  random forest and linear SVM, with nested cross-validation.

> **v3 is in development and has not been released.** The 2021 release (v2.00001, the version
> cited by Hu et al. 2024) is preserved unchanged at the tag `v2.00001` and on the branch
> `legacy-v2`.

## What changed in v3

| | v2 (2021) | v3 |
|---|---|---|
| Install | Seq2Geno's conda envs + a Geno2Pheno **web server** | one workflow (`nextflow run`) with pinned containers; `pip install` for the modelling part |
| Modelling | upload to the Geno2Pheno server | local, `s2g2p evaluate` |
| Input | read pairs (assemblies via a later adaptation) | a sample sheet of assemblies |
| Test-set information | pangenome, indels, compression and k-mer scaling fitted on all genomes | **`noleak` (default)**: every feature fitted on training genomes only |
| New genomes | not supported | `--mode train`, then `--mode predict` |
| Published configurations | — | named presets that reproduce them |

## Presets

A preset fixes every method choice. One YAML file per preset (`src/s2g2p/presets/`) is read by
both the workflow and the Python package. `s2g2p presets` lists them.

| Preset | Status | What it is |
|---|---|---|
| `noleak` | **recommended** | `hu2024`'s tools, flags, grids and selection rule, but the pangenome, alignments and indel calls are built per outer fold from training genomes. Held-out genomes are placed into them one at a time; k-mer scaling and compression are refitted inside every training split. The only preset that can train and predict. |
| `hu2024` | reproduction | Seq2Geno2Pheno exactly as Hu et al. 2024 ran it. 6-mers + Roary GPA + indels, compressed and concatenated; the published grids; per outer fold, the classifier with the best inner macro-F1. Fitted on all genomes, test genomes included. |
| `fast` | experimental | `noleak` without Roary's `-e --mafft` core alignment, which no feature reads (about two thirds of Roary's runtime). |
| `g2p` | experimental | Seq2Geno tables modelled under Geno2Pheno's protocol: separate feature sets, L1-normalised k-mers, group-aware inner CV. The grid is a stand-in until Geno2Pheno's own grid is recovered. |

## Quick start

```bash
# sample sheet: sample, assembly, [group], one column per phenotype (R/S or 1/0; empty = unlabelled)
nextflow run hzi-bifo/seq2geno2pheno -profile singularity \
    --input samples.tsv --preset noleak --outdir results

# reproduce the published configuration on fixed folds ('{phenotype}' is replaced per column)
nextflow run hzi-bifo/seq2geno2pheno -profile singularity,slurm \
    --input samples.tsv --preset hu2024 --folds 'folds/{phenotype}_random_cv.json'

# a final model, then predictions for new genomes
nextflow run hzi-bifo/seq2geno2pheno --input train.tsv --mode train   --outdir model_run
nextflow run hzi-bifo/seq2geno2pheno --input new.tsv   --mode predict --model_dir model_run/model
```

The modelling half also runs on its own, on tables you already have:

```bash
pip install .            # Python >= 3.10
s2g2p evaluate --preset noleak --samplesheet samples.tsv --folds folds.tsv \
    --kmer kmers.tsv --gpa 'fold_{fold}/gpa.tsv' --indel 'fold_{fold}/indel.tsv' --out-dir results
```

Each phenotype gets `results/<phenotype>/`:

| File | Contents |
|---|---|
| `summary.tsv` | metric × {mean, std, weighted_mean, weighted_std} over folds |
| `folds.tsv` | the selected candidate per outer fold |
| `candidates.tsv` | every feature-set × classifier candidate per fold |
| `predictions.tsv` | held-out predictions |
| `run.json` | preset, versions, inputs |
| `manifest.tsv` | index of the files above |

Hu et al. reported the support-weighted mean for homology-aware (KMA) folds and the plain mean
otherwise; both are always written.

## Layout

```
main.nf, modules/, conf/   Nextflow workflow (Seq2Geno: Prokka -> Roary -> GPA / indels / k-mers)
src/s2g2p/                 Python package and CLI (features, placement, nested CV, train/predict)
src/s2g2p/presets/         the presets
containers/                image definitions
docs/                      design notes, the published method, legacy material
tests/                     unit tests; each port is checked against the published code
dev/                       development environment and test runner
```

## Development

```bash
conda env create -f dev/environment-dev.yml && conda activate s2g2p-dev
PYTHONPATH=src pytest -q
nextflow run . -profile test -stub-run          # wiring only, no tools
```

See [ROADMAP.md](ROADMAP.md) for what v3.0 still needs, and [docs/design.md](docs/design.md)
for why each preset is built the way it is.

## Citation

v3 has not been released or published. Until it is, cite the version you used. For the
published configuration:

- Hu K. *et al.* Assessing computational predictions of antimicrobial resistance phenotypes from
  microbial genomes. *Briefings in Bioinformatics* 25(3):bbae206 (2024).
- Asgari E., Kuo T.-H., Bremges A., Robertson G., Weimann A., McHardy A. C. Seq2Geno2Pheno,
  v2.00001 (2021).

## License

GPL-3.0-or-later, as v2. The ported Seq2Geno code (`featCompress.py`, `vcf2indel.py`,
`generate_indel_features.py`, `gene_clusters2multi_fasta.py`, `filter.awk`) is © 2021
Tzu-Hao Kuo and Aaron Weimann, GPL-3.0-or-later.
