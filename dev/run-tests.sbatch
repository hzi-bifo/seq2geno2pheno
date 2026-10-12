#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-tests
#SBATCH --account=bifo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --nice=30000
#
# Unit tests and a Nextflow stub run, in the s2g2p-dev env (dev/environment-dev.yml).
# HZI clusters: submit from the repo root with your partition/account, e.g.
#   sbatch -o log/%x_%j.out dev/run-tests.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate s2g2p-dev
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m pytest -q "$@"
if [[ "${STUB:-1}" == 1 ]]; then
  NXF_OFFLINE=true bash dev/stub-tests.sh
fi
