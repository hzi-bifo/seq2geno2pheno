# Tutorial: vancomycin resistance in *Enterococcus faecium*

This walks through a complete analysis on 30 public genomes. Most of the time goes to Roary: one
pangenome for `hu2024`, one per fold for `noleak`, one for the final model. Measured on a single
8-core node, `hu2024` takes about 2.5 hours; on a SLURM cluster `noleak` takes about 3 hours, with its
five folds in parallel ([What you should get](#what-you-should-get) lists every step). You will:

1. install the pipeline (no containers, no root);
2. evaluate the published method (`hu2024`) and the recommended alternative (`noleak`);
3. read the results;
4. train a final model and predict five genomes it has never seen;
5. run the same thing on a SLURM cluster.

The genomes and labels are a subset of the PATRIC data of Hu *et al.* 2024: 15 resistant, 15
susceptible, plus 5 more for the prediction step.

## 1. Install

You need Linux x86-64, [conda](https://conda-forge.org/download/) (or mamba / micromamba),
internet access and about 6 GB of disk.

```bash
# Nextflow (once)
conda create -n nextflow -c conda-forge -c bioconda nextflow openjdk=17
conda activate nextflow

# the pipeline
git clone https://github.com/hzi-bifo/seq2geno2pheno.git
cd seq2geno2pheno
bash envs/install.sh            # ~15-20 min, once
```

Tested with Nextflow 24.04.4 and 26.10.0 (the version `conda create` installs today).

`envs/install.sh` creates one conda environment per tool under `envs/installed/`, each from a
lock file, so the versions are exactly those the published method used (Prokka 1.14.6, Seq2Geno's
Roary 3.8.2 fork with BLAST+ 2.4.0, mafft 7.310, msa2vcf 1.0, scikit-learn 0.24.1 for the
published models). It ends with a self-check. Running it again only redoes what changed.

To install somewhere else, use `bash envs/install.sh --prefix /path/to/envs` and pass
`--env_dir /path/to/envs` to every run below. On a cluster, that path must be visible from the
compute nodes.

## 2. Get the genomes

```bash
bash tutorial/fetch_genomes.sh
```

This downloads 35 assemblies from [BV-BRC](https://www.bv-brc.org) into `tutorial/genomes/` and
checks them against the copies Hu *et al.* used.

The sample sheet is `tutorial/samples.tsv`:

```
sample      assembly                 vancomycin
1352.10167  genomes/1352.10167.fna   0
1352.11265  genomes/1352.11265.fna   1
...
```

- `sample`: any unique name.
- `assembly`: a FASTA file, absolute or relative to the sheet.
- Every other column is a phenotype: `1`/`0` or `R`/`S`, empty for unknown.
- An optional `group` column keeps related genomes in the same inner fold (`g2p` preset).

## 3. Evaluate the published method

```bash
nextflow run . -profile conda \
    --input tutorial/samples.tsv --preset hu2024 --n_folds 5 \
    --outdir tutorial-results/hu2024
```

`hu2024` is Seq2Geno2Pheno exactly as Hu *et al.* 2024 ran it:

1. Prokka annotation;
2. one Roary pangenome over all genomes;
3. gene presence/absence, indels and 6-mers, compressed and concatenated;
4. SVM, logistic regression, random forest and linear SVM with the published grids, and per
   outer fold the one with the best inner-CV macro-F1.

The models run in the published scikit-learn 0.24.1 environment.

## 4. Evaluate the recommended alternative

```bash
nextflow run . -profile conda \
    --input tutorial/samples.tsv --preset noleak --n_folds 5 \
    --outdir tutorial-results/noleak -resume
```

`hu2024` builds the pangenome, the indel alignments, the compression and the k-mer scaling from
all genomes, test genomes included. No labels leak, but the test genomes shape the features.

`noleak` builds all of them from each fold's training genomes only. Each held-out genome is then
placed into that fold's pangenome as a new genome would be. This costs one Roary run per fold.

`-resume` reuses the Prokka annotations from step 3. It resumes the most recent run launched from
the same directory, so run the steps in order.

Other presets:

| Preset | Use it for |
|---|---|
| `hu2024` | the published method, to compare with published numbers |
| `noleak` | **an honest estimate of performance on new genomes; also the only preset that can train a model for prediction** |
| `fast` | `noleak` without Roary's core-gene alignment (`-e --mafft`), which no feature reads. On these 24-genome folds Roary takes 2-3 min instead of 73-103 min; on large species the saving is about 3x |
| `g2p` | experimental: Geno2Pheno's modelling protocol (separate feature sets, group-aware CV) |

`s2g2p presets` (or `cat src/s2g2p/presets/<name>.yml`) shows every choice a preset makes.

To use your own folds instead of seeded stratified ones, pass `--folds folds.json` (a JSON list
of lists of sample names, or a TSV of `sample`, `fold`). With several phenotype columns, use
`--folds 'folds/{phenotype}.json'`.

## 5. Read the results

Each run writes `evaluation/results/<phenotype>/`:

| File | What it holds |
|---|---|
| `summary.tsv` | macro-F1, accuracy, precision, recall, MCC, AUC: mean and std over folds |
| `folds.tsv` | per outer fold, the chosen classifier, its hyperparameters and inner-CV score |
| `candidates.tsv` | every classifier tried on every fold |
| `predictions.tsv` | the held-out prediction for every genome (`y_true`, `y_pred`) |
| `run.json` | the version, the preset and the inputs, enough to rerun |
| `manifest.tsv` | an index of these files |

```bash
column -t tutorial-results/hu2024/evaluation/results/vancomycin/summary.tsv | head -3
column -t tutorial-results/noleak/evaluation/results/vancomycin/summary.tsv | head -3
```

With 30 genomes and 5 folds, each test fold has six genomes, so a single error moves a fold's
macro-F1 by about 0.15. Treat the differences between presets here as illustration. Hu *et al.*
used all 277 labelled genomes, and that run reproduces their published macro-F1
(`docs/validation/`).

The feature tables are in `features/` (`hu2024`). For `noleak`, each fold's tables are in that
fold's work directory, since they differ per fold. `pipeline_info/` has Nextflow's timeline and
resource report.

### What you should get

These are the numbers from two clean runs of this tutorial, on 2026-10-13 and 2026-10-14. Each
started from a fresh clone and a new conda environment, ran `envs/install.sh`, then every command
on this page. Nothing was reused from an earlier installation (`dev/cleanroom.sbatch`). Macro-F1
is the mean ± standard deviation over the 5 outer folds.

| Run | Command | Macro-F1, run 1 | Macro-F1, run 2 | Wall time, run 1 / run 2 |
|---|---|---:|---:|---:|
| Published method | `--preset hu2024` | 0.865 ± 0.140 | 0.830 ± 0.118 | 1 h 43 min / 1 h 40 min |
| Same, local executor, 8 CPUs | `-profile conda` only | 0.829 ± 0.210 | 0.796 ± 0.072 | 2 h 22 min / 2 h 40 min |
| Recommended | `--preset noleak` | 0.764 ± 0.089 | 0.830 ± 0.118 | 2 h 44 min / 2 h 59 min |
| `noleak` without the core alignment | `--preset fast` | 0.790 ± 0.150 | 0.830 ± 0.118 | 1 h 21 min / 1 h 38 min |
| Geno2Pheno protocol (experimental) | `--preset g2p` | 0.762 ± 0.194 | 0.730 ± 0.194 | 2.5 min / 2.4 min (reuses `hu2024`'s features) |
| `noleak` on your own folds | `--folds folds.tsv` | 0.597 ± 0.317 | 0.644 ± 0.266 | 2 h 51 min / 3 h 25 min |

Times are for SLURM task jobs queued at low priority, except the local run, and include
Prokka. Each step on its own:

| Step | Measured |
|---|---|
| `envs/install.sh` | 16-18 min |
| `tutorial/fetch_genomes.sh` | 1 min; all 35 genomes byte-identical to the benchmark's copies |
| Prokka | 2.5-3 min per genome, genomes in parallel |
| Roary, 30 genomes (`hu2024`) | 56 min |
| Roary, 24 genomes per fold (`noleak`, folds in parallel) | 73-103 min |
| Roary, 24 genomes per fold (`fast`) | 2-3 min |
| Indels | 16-33 min (`hu2024`), 25-77 min per fold (`noleak`) |
| Placing the held-out genomes (`noleak`) | 2-35 min per fold |
| Models | 2-5 min |

How to read these numbers:

- **They are noisy, and the differences between rows are not findings.** A test fold here has six
  genomes, so one wrong prediction moves that fold's macro-F1 by about 0.15.
- **Every `hu2024` result here comes from the same code on the same data.** The results differ
  because Roary is not deterministic. Between the SLURM and local runs of run 1:
  - Roary made 7,084 gene families in one and 7,083 in the other.
  - That changed 4 of about 1,200 compressed gene-presence columns.
  - The change flipped the classifier chosen in one fold, and 3 of the 30 held-out predictions
    differ.

  Run 2 shows the same effect differently: three presets landed on exactly the same score. Your
  own run will land near these numbers, not on them.
- **The `--folds` row uses a different split**, five folds alternating within each class, so its
  folds hold different genomes. Its spread of ± 0.27-0.32 is the same small-sample effect.
- **On the full 277 genomes** of Hu *et al.*, `hu2024` reproduces the published macro-F1
  (`docs/validation/`).

## 6. Train a model and predict new genomes

```bash
nextflow run . -profile conda --mode train \
    --input tutorial/samples.tsv --preset noleak --outdir tutorial-results/model -resume

nextflow run . -profile conda --mode predict \
    --input tutorial/new_genomes.tsv --model_dir tutorial-results/model/model \
    --outdir tutorial-results/predict
```

Training builds the pangenome from all 30 genomes. It selects the classifier by cross-validation
exactly as in `evaluate`, then refits it on all of them.

`model/` is self-contained: the pangenome, the alignments and one `<phenotype>.joblib` per
phenotype. Copy it anywhere.

Prediction annotates the new genomes and places them into that pangenome, without changing it.
It writes `predictions/vancomycin.tsv`, with the columns `sample`, `prediction` (1 = resistant),
`decision` (the classifier's score), `phenotype` and `classifier`. Compare with the labels BV-BRC
has for them:

```bash
join -t $'\t' <(tail -n +2 tutorial-results/predict/predictions/vancomycin.tsv | sort) \
              <(tail -n +2 tutorial/new_genomes_truth.tsv | sort)
```

In the clean run, all five predictions matched. Training took 81 min and prediction 6 min, each
including annotation and placement:

| Sample | Prediction | Decision | Classifier | BV-BRC label |
|---|---:|---:|---|---:|
| 1352.10148 | 0 | -0.95 | lsvm | 0 |
| 1352.10474 | 0 | -0.51 | lsvm | 0 |
| 1352.11262 | 1 | 0.39 | lsvm | 1 |
| 1352.11266 | 1 | 1.64 | lsvm | 1 |
| 1352.11327 | 1 | 0.33 | lsvm | 1 |

Only `noleak` and `fast` can train: their features are defined without the test genomes, so a
new genome is handled exactly as in evaluation.

## 7. On a SLURM cluster

Add `slurm` to the profile. Every step then becomes its own job:

```bash
nextflow run . -profile conda,slurm \
    --input tutorial/samples.tsv --preset noleak --outdir tutorial-results/noleak \
    --max_cpus 16 --max_memory 64.GB --max_time 24.h
```

- **Launch Nextflow itself as a small job** (2 CPUs, 8 GB, as long as the whole run), not on the
  login node. `dev/run-e2e.sbatch` is an example.
- **Installing** can also run as a job: `nextflow run . -profile conda,slurm --mode install`.
  The node needs internet access.
- **Partition, account and other `sbatch` options:** put them in a config file and pass it with
  `-c`, for example `process.queue = 'cpu'` and
  `process.clusterOptions = '--account=myproj'`. `conf/hzi.config` is the HZI clusters' version.

| Step | 277 *E. faecium* genomes, measured |
|---|---|
| Prokka | ~1.5 min per genome on 8 CPUs, genomes in parallel |
| Roary | 9.5 h on 16 CPUs (once for `hu2024`; once per fold for `noleak`, on 90% of the genomes) |
| Indels | 1.8 h on 16 CPUs |
| Models | 1.8 h on 16 CPUs (`hu2024`) |

## 8. Only the modelling

If you already have genotype tables (genome × feature TSVs), you do not need the workflow:

```bash
pip install git+https://github.com/hzi-bifo/seq2geno2pheno.git
s2g2p evaluate --preset hu2024 --samplesheet samples.tsv --phenotype vancomycin \
    --folds folds.json --kmer kmers.tsv --gpa gpa.tsv --indel indel.tsv --out-dir results
```

`s2g2p --help` lists the other commands: feature building, compression, folds, train and
predict.

## Troubleshooting

- **`-profile conda: no environments under --env_dir`**: run `bash envs/install.sh` (or pass
  the `--env_dir` you installed to).
- **`BLAST+ 2.4.0 needs system libraries this host lacks`**: BLAST+ 2.4.0 is NCBI's binary and
  needs the system's `libnsl.so.1`. Install it with `dnf install libnsl` (RHEL, Alma, Rocky) or
  your distribution's package that ships it.
- **A step fails**: look at `.command.err` in the work directory Nextflow names. After fixing
  the cause, rerun with `-resume`; finished steps are not repeated.
