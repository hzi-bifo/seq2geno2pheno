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
| Folds | random or tree-based, on the server | `--mode folds`: random, phylogeny-aware or homology-aware (KMA), as Hu et al. 2024 |
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

What `envs/install.sh` puts under `envs/installed/` (or `--prefix DIR`, then `--env_dir DIR` in
every run):

| Environment | For | Contents |
|---|---|---|
| `prokka`, `roary-perl`, `roary-deps`, `mafft`, `msa2vcf` + Roary, BLAST+ 2.4.0, GNU parallel | features (step 5) | Seq2Geno's tools at the published versions |
| `s2g2p`, `hu2024ml` | models | this package; scikit-learn 0.24.1 for `hu2024` |
| `folds` | `--mode folds` with `phylogeny` or `kma` (step 4) | R 4.2.3, ape 5.6-2, phangorn 2.10.0 (`nj_tree`); `kma_clustering`, which the script compiles from `envs/vendor/kma_clustering.c` with the environment's own C compiler and zlib |

- **Re-running it is cheap and safe:** an environment whose lock file has not changed is kept,
  so on an existing install it only adds what is new (e.g. `folds`). A workflow run that needs a
  tool the install lacks stops at the start and says so.
- **It needs** conda, mamba or micromamba, internet access (conda-forge, bioconda, github.com,
  ftp.ncbi.nlm.nih.gov, ftp.gnu.org, CPAN), Linux x86-64, and `libnsl.so.1` for BLAST+ 2.4.0
  (RHEL/Alma/Rocky: `dnf install libnsl`). It checks every tool's version before it finishes.
- **On a cluster**, run it where the compute nodes see the result, as a job if the login node is
  not for compute: `nextflow run . -profile conda,slurm --mode install`.

## Step by step

[The tutorial](docs/tutorial.md) runs every step on 30 public genomes. In short:

```bash
nextflow run . -profile conda --input samples.tsv --mode folds --fold_method kma \
    --outdir run -work-dir /shared/scratch/run-work                         # step 4 (optional)
nextflow run . -profile conda --input samples.tsv --preset noleak \
    --folds 'run/folds/{phenotype}_kma_cv.json' --outdir run -work-dir /shared/scratch/run-work
```

### 1. Sequence files

One nucleotide FASTA per genome: an assembly (contigs or a complete genome), **uncompressed**,
since Prokka does not read `.gz`. The file name does not matter; the sample sheet names the
genome. Keep the files on storage every compute node can read: tasks link to them rather than
copy them.

### 2. Sample sheet

A TSV (or `.csv`) with one row per genome:

```
sample      assembly                 vancomycin  ciprofloxacin
1352.10167  genomes/1352.10167.fna   0           R
1352.11265  genomes/1352.11265.fna   1
```

- **`sample`**: a unique id of `A-Z a-z 0-9 . _ -`. It becomes Prokka's locus tag and the row
  name in every table.
- **`assembly`**: an absolute path, or one relative to the sheet's own folder. Relative paths
  with the sheet next to the genomes make the sheet work on any machine that mounts them.
- **Every other column is a phenotype**: `1`/`0`, `R`/`S`, `resistant`/`susceptible`. An empty
  cell means not labelled for that phenotype: the genome is still annotated and in the
  pangenome, but not in that phenotype's model or folds.
- **`group`** (optional) keeps related genomes in one inner fold (`g2p` preset).
- **Row order** only matters for `kma` folds, because kma_clustering is greedy: to regenerate
  published KMA folds, list the genomes in the order of the published phenotype file.

### 3. Where the files go

- **`-work-dir`** (default: `./work` where you launch) holds every task's files. Prokka's
  annotations, every Roary pangenome (one per fold for `noleak`) and the alignments exist only
  here; `--outdir` gets copies of the tables and results (step 6). For 277 genomes, `hu2024`'s
  work directory was 6.8 GB against 164 MB in `--outdir`, and `noleak` builds *k* times as many
  pangenomes. Choose it explicitly, on shared storage (not a node's `/tmp`), and keep it until
  the run is finished and checked: `-resume` reuses it, and without it Prokka and Roary run
  again. The workflow warns when no `-work-dir` is given.
- **`--outdir`** holds the results. Two runs can share it (e.g. `--mode folds`, then
  `evaluate`); two *concurrent* runs need separate `-work-dir`s and Nextflow caches
  (`NXF_CACHE_DIR`, as `dev/run-e2e.sbatch` sets).

### 4. Folds

Every evaluation splits each phenotype's labelled genomes into *k* outer folds (`--n_folds`,
default 10); each fold is held out once. Three ways to get them:

- **Your own**: `--folds folds.json`, a JSON list of lists of sample ids, or a TSV of `sample`
  and `fold` (0..k-1). With several phenotypes, `--folds 'folds/{phenotype}.json'`. Genomes in
  the file but not labelled are ignored; a labelled genome missing from it is an error. **To
  compare with Hu *et al.* 2024, use their published fold files** ([Mendeley
  data](https://data.mendeley.com/datasets/6vc2msmsxj/2), or the
  [AMR_benchmarking](https://github.com/hzi-bifo/AMR_benchmarking) repository): none of the
  three can be regenerated byte for byte (below).
- **None**: seeded stratified random folds (`--seed`, default 0), made inside the run.
- **`--mode folds --fold_method random|phylogeny|kma`**: writes the folds to files first, with
  the three methods of Hu *et al.*, for a later `--folds`:

| Method | What it does | Fold sizes | Versus the published folds |
|---|---|---|---|
| `random` | stratified *k*-fold, seeded (`--seed`) | equal, same R/S ratio | Geno2Pheno made the published ones, stratified, with an unrecorded order and seed |
| `phylogeny` | Prokka and Roary (`-e --mafft`) over the whole sheet; Hamming distances of the core-gene alignment; phangorn neighbour-joining tree (`nj_tree`, as Hu *et al.*); the labelled genomes in leaf order, cut into *k* contiguous blocks. `--tree file.newick` uses an existing tree instead. | equal, not stratified | Geno2Pheno split Hu *et al.*'s trees with code that was never released. This rebuilds it from what it is documented and observed to do ([docs/folds.md](docs/folds.md)); Roary is not deterministic, so the tree, and the folds, differ slightly |
| `kma` | per phenotype, its labelled genomes with contigs joined by 16 Ns; `kma_clustering -k 16 -ht/-hq <identity> -NI`; clusters allocated to folds by Hu *et al.*'s refinement of Aytan-Aktug's method (seed 42) | uneven: a cluster is never split | the published procedure, ported and tested against its code ([docs/folds.md](docs/folds.md)) |

```bash
nextflow run . -profile conda --mode folds --fold_method phylogeny --input samples.tsv --outdir run
nextflow run . -profile conda --mode folds --fold_method kma --kma_identity 0.98 \
    --input ngonorrhoeae.tsv --outdir run
```

- **`--kma_identity`** is kma_clustering's template and query identity (`-ht`, `-hq`), 0.9 by
  default. Hu *et al.* used 0.98 for *N. gonorrhoeae*, and for *M. tuberculosis* 0.99
  (amikacin, ethambutol, ethiomide, ethionamide, kanamycin, ofloxacin) or 0.98 (the other
  drugs); run those as separate `--phenotype`s.
- `phylogeny` and `kma` need `-profile conda` for now: no container image carries R/phangorn or
  kma_clustering yet.
- Output: `<outdir>/folds/<phenotype>_<method>_cv.json`, a `manifest.tsv` (genomes, fold
  sizes and resistant genomes per fold), plus `folds/tree/nj_tree.newick` or
  `folds/kma/<phenotype>.clusters.txt`. The same step without the workflow:
  `s2g2p folds --method ... --tree ... | --kma-clusters ... --out-dir folds`.

### 5. Evaluate

```bash
nextflow run . -profile conda --input samples.tsv --preset noleak --folds 'run/folds/{phenotype}_kma_cv.json' \
    --outdir run -work-dir /shared/scratch/run-work
```

For each phenotype and each outer fold:

1. **Features**: Prokka annotates every genome once; Roary builds the pangenome; gene
   presence/absence, indels and 6-mer counts are tabulated. `noleak` builds the pangenome from
   that fold's training genomes and places the held-out genomes into it; `hu2024` builds one
   pangenome over all genomes, as published.
2. **Model selection, on the training genomes only**: each of SVM, logistic regression, random
   forest and linear SVM is tuned over the published grid by an inner cross-validation over the
   other folds; the candidate with the best inner macro-F1 is chosen.
3. **Test**: the chosen model predicts the held-out fold.

The modelling half also runs on its own, on tables you already have:

```bash
pip install .            # Python >= 3.10
s2g2p evaluate --preset noleak --samplesheet samples.tsv --folds folds.tsv \
    --kmer kmers.tsv --gpa 'fold_{fold}/gpa.tsv' --indel 'fold_{fold}/indel.tsv' --out-dir results
```

### 6. Results

```
<outdir>/
  evaluation/results/<phenotype>/   summary.tsv, folds.tsv, candidates.tsv, predictions.tsv, run.json, manifest.tsv
  features/                         kmers.tsv; hu2024 also gpa_NONRDNT.tsv, indel_NONRDNT.tsv
  splits/<phenotype>/               noleak: each fold's train/test genome lists, manifest.tsv
  folds/                            --mode folds (step 4)
  model/, predictions/              --mode train / predict (step 7)
  annotations/                      Prokka's output, only with --publish_annotations
  pipeline_info/                    Nextflow's timeline, report and trace
```

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

### 7. A final model, and predictions for new genomes

```bash
nextflow run . -profile conda --input train.tsv --mode train   --outdir model_run
nextflow run . -profile conda --input new.tsv   --mode predict --model_dir model_run/model
```

`train` fits one model per phenotype on all labelled genomes (`noleak`-style presets only) and
saves it with its pangenome; `predict` places new genomes into that pangenome and writes
`predictions/<phenotype>.tsv`.

### 8. Running on SLURM

Add `slurm` to the profile and Nextflow's SLURM executor submits every task (each Prokka genome,
each Roary pangenome, each fold's evaluation) as its own job:

```bash
nextflow run . -profile conda,slurm --input samples.tsv --preset noleak --outdir results \
    --max_cpus 16 --max_memory 64.GB --max_time 24.h
```

- **Resources** come from the process labels in `conf/base.config` (`process_single`, `_low`,
  `_medium`, `_high`; Roary and the indel step are `_high`: 16 CPUs, 64 GB, doubled on the one
  retry). Every request is capped by `--max_cpus`, `--max_memory` and `--max_time` (defaults
  32, 128 GB, 48 h). Tasks killed by the scheduler (exit 137, 140, 143, ...) are retried once.
- **Partition, account and other `sbatch` options** go in a config file passed with `-c`, e.g.
  `process.queue = 'cpu'`, `process.clusterOptions = '--account=myproj'`,
  `executor.queueSize = 10`.
- **Run Nextflow itself as a small job** (2 CPUs, 8 GB, as long as the whole run), not on the
  login node; `dev/run-e2e.sbatch` is an example. If that job hits the walltime limit, which a
  large species can (Roary alone took 39 h on 4,757 *S. pneumoniae*), resubmit it with
  `-resume` and finished tasks are reused.
- **Installing** can run as a job too: `nextflow run . -profile conda,slurm --mode install`
  (the node needs internet access).

**HZI clusters (bioinf, metis):** add the `hzi` profile, `-profile hzi,slurm,conda`
(`conf/hzi.config`). It charges every job to `bifo`, keeps at most 10 jobs queued
(`--slurm_queue_size`), and passes `--slurm_extra` (default `--nice=30000`, so background runs
queue behind interactive work) to every task. The partition differs per cluster:
`--slurm_partition cpu` on bioinf, `cpugpu` on metis. The conda tools are built for AVX2 and die
with SIGILL (exit 132) on older nodes, so add the placement to `--slurm_extra`:
`'--nice=30000 --exclude=bioinf009'` on bioinf, `'--nice=30000 --constraint=avx2'` on metis.

```bash
# bioinf: dev/run-e2e.sbatch launches the driver job (its task placement is bioinf's exclude)
env INPUT=samples.tsv PRESET=noleak OUT=results-noleak \
    sbatch -o log/%x_%j.out dev/run-e2e.sbatch
```

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
src/s2g2p/R/               nj_tree.R (phylogeny-aware folds)
envs/                      lock files and install.sh: the tools, without containers
envs/vendor/                kma_clustering.c, as Hu et al. 2024 shipped it (Apache-2.0)
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
