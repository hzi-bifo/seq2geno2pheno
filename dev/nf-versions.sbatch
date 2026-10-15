#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-nfversions
#SBATCH --account=bifo
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#
# The workflow must compile and wire under every Nextflow the README accepts (>= 24.04) and under
# the newest, whose strict syntax is the default from 26.x. For each version: `nextflow lint`
# (strict-syntax errors; 25.04+ only) and the stub run of -profile test.
#   env ENVROOT=/shared/dir VERSIONS="24.04.4 26.10.0" sbatch -o log/%x_%j.out dev/nf-versions.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
: "${ENVROOT:?a directory for the per-version Nextflow environments}"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=safe.directory GIT_CONFIG_VALUE_0="$PWD"
rc=0
for v in ${VERSIONS:-24.04.4 26.10.0}; do
    e="$ENVROOT/nextflow-$v"
    [ -x "$e/bin/nextflow" ] || conda create -y -q -p "$e" -c conda-forge -c bioconda "nextflow=$v" openjdk=17 > /dev/null || { echo "$v: env failed"; rc=1; continue; }
    conda activate "$e"
    export NXF_HOME="$ENVROOT/.nextflow-home-$v"
    echo "== Nextflow $v"
    if nextflow lint -h > /dev/null 2>&1; then
        nextflow lint -o concise main.nf modules conf nextflow.config 2>&1 | grep -E "^(Error|Warn)|files? had" ; [ "${PIPESTATUS[0]}" = 0 ] || { echo "$v: lint FAILED"; rc=1; }
    else
        echo "(no lint in $v)"
    fi
    # the three shapes of a run: cohort (hu2024), per fold (noleak), and a final model
    for args in "--preset hu2024" "--preset noleak" "--preset noleak --mode train"; do
        tag="$v ${args// /_}"; log="stub-$v-${args// /_}.log"
        rm -rf "work-stub-$v"
        # shellcheck disable=SC2086
        nextflow run . -ansi-log false -profile test -stub-run -work-dir "work-stub-$v" \
            --outdir "results-stub-$v" $args > "$log" 2>&1 \
            && echo "$tag: stub run OK" \
            || { echo "$tag: stub run FAILED"; tail -20 "$log"; rc=1; }
        rm -rf "work-stub-$v" "results-stub-$v" .nextflow
    done
    conda deactivate
done
exit $rc
