<!--
SPDX-FileCopyrightText: 2021 Ehsaneddin Asgari and Tzu-Hao Kuo
SPDX-FileCopyrightText: 2026 Hadi Foroughmand Araabi

SPDX-License-Identifier: LicenseRef-ASL-1.0
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
| Install | Seq2Geno's conda envs + a Geno2Pheno **web server** | `envs/install.sh`: pinned conda environments, no containers, no root; `pip install` for the modelling part |
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
| `hu2024` | **original** (published method) | Seq2Geno2Pheno exactly as Hu et al. 2024 ran it. 6-mers + Roary GPA + indels, compressed and concatenated; the published grids; per outer fold, the classifier with the best inner macro-F1. Fitted on all genomes, test genomes included. |
| `fast` | experimental | `noleak` without Roary's `-e --mafft` core alignment, which no feature reads (about two thirds of Roary's runtime). |
| `g2p` | experimental | Seq2Geno tables modelled under Geno2Pheno's protocol: separate feature sets, L1-normalised k-mers, group-aware inner CV. The grid is a stand-in until Geno2Pheno's own grid is recovered. |

## Install

Linux x86-64, [conda](https://conda-forge.org/download/) (or mamba/micromamba) and
[Nextflow](https://www.nextflow.io) ≥ 24.04:

```bash
conda create -n nextflow -c conda-forge -c bioconda nextflow openjdk=17 && conda activate nextflow
git clone https://github.com/hzi-bifo/seq2geno2pheno.git && cd seq2geno2pheno
bash envs/install.sh          # once: one locked environment per tool, under envs/installed/
```

Every tool is installed at the version the published method used, from lock files (`envs/`).
Container images are an alternative (`containers/`, `-profile apptainer|singularity|docker`).

## Quick start

**[The tutorial](docs/tutorial.md) runs everything below on 30 public genomes.**

```bash
# sample sheet: sample, assembly, [group], one column per phenotype (R/S or 1/0; empty = unlabelled)

# the original (published) method
nextflow run . -profile conda --input samples.tsv --preset hu2024 --outdir results-hu2024

# the recommended alternative: no test-set information in any feature
nextflow run . -profile conda --input samples.tsv --preset noleak --outdir results-noleak

# fixed folds ('{phenotype}' is replaced per column); on a cluster add slurm to the profile
nextflow run . -profile conda,slurm --input samples.tsv --preset hu2024 \
    --folds 'folds/{phenotype}_random_cv.json' --outdir results

# a final model, then predictions for new genomes
nextflow run . -profile conda --input train.tsv --mode train   --outdir model_run
nextflow run . -profile conda --input new.tsv   --mode predict --model_dir model_run/model
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
envs/                      lock files and install.sh: the tools, without containers
containers/                optional image definitions (same lock files)
tutorial/                  the tutorial's sample sheets and download script
docs/                      tutorial, design notes, validation, release policy, legacy material
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

v3 has not been released or published. Until it is, cite the version you used
(`CITATION.cff`; GitHub's "Cite this repository" reads it). For the published configuration:

- Hu K. *et al.* Assessing computational predictions of antimicrobial resistance phenotypes from
  microbial genomes. *Briefings in Bioinformatics* 25(3):bbae206 (2024).
- Asgari E., Kuo T.-H., Bremges A., Robertson G., Weimann A., McHardy A. C. Seq2Geno2Pheno,
  v2.00001 (2021).

## License

From v3, the [Academic Software License v1.0](LICENSE) (ASL): free for academic, non-commercial
use, with GPL-style reciprocity; commercial use needs a separate licence from the licensor
([NOTICE](NOTICE)). v2.00001 and earlier stay under GPL-3.0-or-later. NOTICE also lists the code
ported from Seq2Geno and from Hu et al. 2024, and the tools `envs/install.sh` downloads, which
keep their own licences.
