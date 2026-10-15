#!/usr/bin/env bash
# Nextflow stub runs: the wiring of every mode and of -profile conda, with no tool installed
# (each process runs its stub block). Needs nextflow on PATH; run from the repository root.
# Used by dev/run-tests.sbatch and by CI (.github/workflows/ci.yml).
set -euo pipefail
export NXF_HOME="${NXF_HOME:-$PWD/.nextflow-home}"
for preset in noleak hu2024; do
  nextflow -q run . -profile test -stub-run --preset "$preset" --outdir "results-test/$preset" -work-dir "work-test"
  echo "stub run $preset: OK"
done
nextflow -q run . -profile test -stub-run --mode train --outdir results-test/train -work-dir work-test
echo "stub run train: OK"
nextflow -q run . -profile test -stub-run --mode predict --model_dir results-test/train/model \
    --outdir results-test/predict -work-dir work-test
echo "stub run predict: OK"
nextflow -q run . -profile test -stub-run --mode install --outdir results-test/install -work-dir work-test
echo "stub run install: OK"
# -profile conda: every task gets the installed wrappers on PATH; EVALUATE gets the preset's
# model environment first (hu2024: bin-hu2024ml)
for preset in hu2024 noleak; do
  out=results-test/conda-$preset
  nextflow -q run . -profile test,conda -stub-run --preset "$preset" --env_dir /nonexistent/envs \
      --outdir "$out" -work-dir "work-test"
  # task names carry the workflow prefix (run_pipeline:EVALUATE)
  run=$(awk -F'\t' 'NR>1 && $4 ~ /(^|:)EVALUATE$/ {print $2; exit}' "$out/pipeline_info/trace.tsv")
  [[ -n $run ]] || { echo "conda profile: no EVALUATE task in $out/pipeline_info/trace.tsv"; exit 1; }
  script=$(ls -d work-test/"$run"*/.command.run)
  if [[ $preset == hu2024 ]]; then want=/nonexistent/envs/bin-hu2024ml:; else want="export PATH=/nonexistent/envs/bin:"; fi
  grep -q "$want" "$script" || { echo "conda profile: EVALUATE ($preset) lacks $want"; exit 1; }
  echo "stub run conda $preset: OK"
done
