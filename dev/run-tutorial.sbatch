#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-tutorial
#SBATCH --account=bifo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --nice=30000
#
# docs/tutorial.md, step by step, as a user would run it -- but with SLURM task jobs (section 7)
# instead of the local executor. Output: tutorial-results/ (gitignored); work: tutorial-work/.
#   sbatch -o log/%x_%j.out dev/run-tutorial.sbatch
# NICE (default 30000) is passed to every task job; ENV_DIR defaults to envs/installed.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate s2g2p-dev
export NXF_HOME="$PWD/.nextflow-home" NXF_OPTS="-Xms512m -Xmx4g" NXF_CACHE_DIR="$PWD/tutorial-results/.nextflow"
common=(-profile hzi,slurm,conda -work-dir tutorial-work --env_dir "${ENV_DIR:-$PWD/envs/installed}"
        --slurm_extra "--nice=${NICE:-30000} --exclude=bioinf009" --max_cpus 8 --max_memory 32.GB --max_time 12.h)

echo "== 2. genomes"; bash tutorial/fetch_genomes.sh
echo "== 3. hu2024"
nextflow run . "${common[@]}" --input tutorial/samples.tsv --preset hu2024 --n_folds 5 \
    --outdir tutorial-results/hu2024
echo "== 4. noleak"
nextflow run . "${common[@]}" --input tutorial/samples.tsv --preset noleak --n_folds 5 \
    --outdir tutorial-results/noleak -resume
echo "== 6. train"
nextflow run . "${common[@]}" --mode train --input tutorial/samples.tsv --preset noleak \
    --outdir tutorial-results/model -resume
echo "== 6. predict"
nextflow run . "${common[@]}" --mode predict --input tutorial/new_genomes.tsv \
    --model_dir tutorial-results/model/model --outdir tutorial-results/predict

echo "== results"
for p in hu2024 noleak; do echo "$p"; head -3 "tutorial-results/$p/evaluation/results/vancomycin/summary.tsv"; done
join -t $'\t' <(tail -n +2 tutorial-results/predict/predictions/vancomycin.tsv | sort) \
              <(tail -n +2 tutorial/new_genomes_truth.tsv | sort)
