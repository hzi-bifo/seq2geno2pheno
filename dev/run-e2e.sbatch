#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-e2e
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --nice=30000
#
# Nextflow driver for an end-to-end run on the HZI clusters; tasks go to SLURM as their own jobs.
# Launch from the repository root under /vol (see conf/hzi.config), e.g.
#   env INPUT=dev-runs/e2e/efaecium_30.tsv PRESET=noleak OUT=dev-runs/e2e/smoke-noleak \
#     sbatch -p cpu -A b_algbio -o log/%x_%j.out dev/run-e2e.sbatch [extra nextflow args]
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate s2g2p-dev
: "${INPUT:?}" "${PRESET:?}" "${OUT:?}"
# One Nextflow cache per output dir: with a shared one, -resume attaches to whichever run
# was launched last from this directory, and a second concurrent run dies on its lock.
export NXF_HOME="$PWD/.nextflow-home" NXF_OPTS="-Xms512m -Xmx4g" NXF_CACHE_DIR="$OUT/.nextflow"
ACCOUNT="${ACCOUNT:-${SLURM_JOB_ACCOUNT:-b_algbio}}"
PARTITION="${PARTITION:-${SLURM_JOB_PARTITION:-cpu}}"
nextflow run . -profile hzi,slurm,apptainer -resume \
  --input "$INPUT" --preset "$PRESET" --outdir "$OUT" \
  --slurm_account "$ACCOUNT" --slurm_partition "$PARTITION" \
  --slurm_extra "--nice=30000 --exclude=bioinf009" \
  -work-dir "$OUT/work" "$@"
