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

> **v3 is a pre-release:** usable and tested (see [Performance](#performance)), with no version
> tag yet ([ROADMAP.md](ROADMAP.md)). The 2021 release (v2.00001, the version cited by Hu et al.
> 2024) is preserved unchanged at the tag `v2.00001` and on the branch `legacy-v2`.

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
[Nextflow](https://www.nextflow.io) ≥ 24.04 (tested: 24.04.4 and 26.10.0):

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

## Performance

### It reproduces the published method

The reference is Seq2Geno2Pheno as Hu *et al.* 2024 ran it, rerun from their published code
([docs/validation/](docs/validation/README.md)).

- **Features: 49 of 49 checks identical** on 7 species, 277 to 4,757 genomes. Each check is one
  of k-mer counts, gene presence/absence, both compressed tables and their groups, the aligned
  families, every indel call and the indel table.
- **Models: the published score, to every printed digit.** On the published features for
  *E. faecium* vancomycin, `--preset hu2024` gives:

  | Folds | Macro-F1, published code | Macro-F1, `--preset hu2024` | Per-classifier scores identical |
  |---|---:|---:|---:|
  | random | 0.973774 | 0.973774 | 40/40 |
  | phylotree | 0.989244 | 0.989244 | 40/40 |
  | KMA (support-weighted) | 0.947719 | 0.947719 | 40/40 |

- **The whole pipeline, from assemblies:** the 277 *E. faecium* genomes, published random folds.
  Macro-F1 is 0.981 ± 0.041, against the published 0.974. The difference comes from Roary, which
  is not deterministic: 97% of gene families have identical members in the two pangenomes, and
  2 of the 277 held-out predictions differ.

### Running time

These times are for the whole pipeline on 277 *E. faecium* genomes, `--preset hu2024`, SLURM,
with up to 16 CPUs per task:

| Step | Wall time |
|---|---:|
| Prokka | 1-2 min per genome, genomes in parallel |
| Roary (pangenome) | 9 h 31 min |
| Indels | 1 h 50 min |
| k-mers, gene presence, compression | under 1 min |
| Nested cross-validation, 10 outer folds | 1 h |

**Roary dominates, and it grows faster than the cohort.** In the Hu *et al.* reproduction, the
same Roary took 4.6 h on 395 genomes (*C. jejuni*) and 39 h on 4,757 (*S. pneumoniae*), at
10-32 CPUs.

**`noleak` builds one pangenome per outer fold**, from that fold's training genomes, so it does
about *k* times `hu2024`'s Roary work. With folds in parallel, the wall time grows less. `fast`
skips Roary's core-gene alignment, which no feature reads. That alignment is about two thirds of
Roary's time on large cohorts, and almost all of it on small ones.

### The presets on the tutorial data

[The tutorial](docs/tutorial.md) uses 30 *E. faecium* genomes (15 resistant, 15 susceptible to
vancomycin) and 5 outer folds. These numbers are from two clean-room runs on SLURM. Each started
from an empty directory: a fresh clone, a new conda environment, `envs/install.sh`, the genome
download, then every command in the tutorial (`dev/cleanroom.sbatch`).

| Run | Macro-F1, run 1 | Macro-F1, run 2 | Wall time, run 1 / run 2 |
|---|---:|---:|---:|
| `--preset hu2024` (published method) | 0.865 ± 0.140 | 0.830 ± 0.118 | 1 h 43 min / 1 h 40 min |
| `--preset hu2024`, local executor, 8 CPUs | 0.829 ± 0.210 | 0.796 ± 0.072 | 2 h 22 min / 2 h 40 min |
| `--preset noleak` (recommended) | 0.764 ± 0.089 | 0.830 ± 0.118 | 2 h 44 min / 2 h 59 min |
| `--preset fast` | 0.790 ± 0.150 | 0.830 ± 0.118 | 1 h 21 min / 1 h 38 min |
| `--preset g2p` (experimental) | 0.762 ± 0.194 | 0.730 ± 0.194 | 2.5 min / 2.4 min |
| `noleak` on user-supplied folds | 0.597 ± 0.317 | 0.644 ± 0.266 | 2 h 51 min / 3 h 25 min |
| `--mode train`, then `--mode predict` on 5 new genomes | 5/5 correct | 5/5 correct | 1 h 21 min + 6 min / 1 h 42 min + 6 min |

The wall times include Prokka. `g2p` reuses `hu2024`'s features.

`envs/install.sh` took 16-17 min, and the 35 genomes took under a minute to download. All 35 are
byte-identical to the copies Hu *et al.* used.

**At 30 genomes, these runs do not rank the presets.** A test fold holds six genomes, so one
wrong prediction moves its macro-F1 by about 0.15. The runs show this:

- The same code on the same data gave `hu2024` 0.865 in one run and 0.830 in the other.
- In run 2, three presets landed on exactly the same score.

The comparison between the presets that matters, at the scale of Hu *et al.*, is in
[ROADMAP.md](ROADMAP.md).

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
