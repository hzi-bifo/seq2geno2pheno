#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-tests
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --nice=30000
#
# Unit tests and a Nextflow stub run, in the s2g2p-dev env (dev/environment-dev.yml).
# HZI clusters: submit from the repo root with your partition/account, e.g.
#   sbatch -p cpu -A b_algbio -o log/%x_%j.out dev/run-tests.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate s2g2p-dev
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m pytest -q "$@"
if [[ "${STUB:-1}" == 1 ]]; then
  export NXF_HOME="$PWD/.nextflow-home" NXF_OFFLINE=true
  for preset in noleak hu2024; do
    nextflow -q run . -profile test -stub-run --preset "$preset" --outdir "results-test/$preset" -work-dir "work-test"
    echo "stub run $preset: OK"
  done
  nextflow -q run . -profile test -stub-run --mode train --outdir results-test/train -work-dir work-test
  echo "stub run train: OK"
  nextflow -q run . -profile test -stub-run --mode predict --model_dir results-test/train/model \
      --outdir results-test/predict -work-dir work-test
  echo "stub run predict: OK"
fi
