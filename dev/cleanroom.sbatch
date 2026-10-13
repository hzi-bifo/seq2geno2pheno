#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-cleanroom
#SBATCH --account=bifo
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --nice=30000
#
# Release check: install and run Seq2Geno2Pheno the way a new user would, in an empty directory,
# from a fresh clone, with a new conda environment and nothing from this checkout's envs/.
# Runs docs/tutorial.md in every setting the tutorial documents and collects what each produced.
#
#   ROOT=/path/to/empty/dir  where everything goes (must not exist yet; visible to compute nodes)
#   SRC=<git url or path>    what to clone (default: this checkout)    REF=<commit> (default HEAD)
#   PARTITION, SLURM_EXTRA   how the task jobs are submitted (default: bioinf's cpu)
#
# metis (no MaxJobs cap; AVX2 nodes only; long QOS for >2 days):
#   env ROOT=/net/hpcprojects/hforough/tmp/s2g2p-cleanroom-$(date +%F) PARTITION=cpugpu \
#       SLURM_EXTRA='--nice=30000 --constraint=avx2' \
#       sbatch -p cpugpu --qos=long --constraint=avx2 --export=ALL -o log/%x_%j.out dev/cleanroom.sbatch
#
# Every run is attempted even if an earlier one failed; ROOT/report/runs.tsv records each exit
# status, and ROOT/report/ holds each run's summary, predictions and trace (manifest.tsv).
set -uo pipefail
here="${SLURM_SUBMIT_DIR:-$PWD}"
: "${ROOT:?set ROOT to a new, empty directory}"
SRC="${SRC:-$here}"; REF="${REF:-$(git -C "$here" rev-parse HEAD)}"
PARTITION="${PARTITION:-cpu}"; SLURM_EXTRA="${SLURM_EXTRA:---nice=30000 --exclude=bioinf009}"
[ -e "$ROOT" ] && { echo "ROOT $ROOT exists; the check needs an empty start" >&2; exit 1; }
mkdir -p "$ROOT/report" && cd "$ROOT" || exit 1
log() { printf '\n== %s  %s\n' "$(date +%T)" "$*"; }
printf 'run\texit\tseconds\toutdir\n' > report/runs.tsv
step() {  # name outdir cmd...: run, time and record a step without stopping the check
    local name=$1 out=$2; shift 2; log "$name"; local t0=$SECONDS
    "$@"; local rc=$?
    printf '%s\t%s\t%s\t%s\n' "$name" "$rc" "$((SECONDS - t0))" "$out" >> report/runs.tsv
    echo "[$name] exit $rc after $((SECONDS - t0)) s"; return $rc
}
{ echo "host $(hostname)  cluster ${SLURM_CLUSTER_NAME:-?}  job ${SLURM_JOB_ID:-local}  $(date -Is)"
  echo "source $SRC @ $REF"; } | tee report/context.txt

# -- 1. install, exactly as README / docs/tutorial.md section 1 ------------------------------------
source "$HOME/miniforge3/etc/profile.d/conda.sh"
step clone - git clone -q "$SRC" seq2geno2pheno || exit 1
git -C seq2geno2pheno checkout -q "$REF" || exit 1
step nextflow-env - conda create -y -q -p "$ROOT/nextflow-env" -c conda-forge -c bioconda nextflow openjdk=17 || exit 1
conda activate "$ROOT/nextflow-env"
cd seq2geno2pheno
export NXF_HOME="$ROOT/.nextflow-home" NXF_OPTS="-Xms512m -Xmx4g"
{ nextflow -version 2>&1 | grep -i version; java -version 2>&1 | head -1; } | tee -a ../report/context.txt
step install - bash envs/install.sh || exit 1
cp envs/installed/manifest.tsv ../report/install_manifest.tsv

# -- 2. genomes -------------------------------------------------------------------------------------
step genomes - bash tutorial/fetch_genomes.sh || exit 1
md5sum -c --quiet tutorial/genomes.md5 > ../report/genomes_md5.txt 2>&1 && echo "all identical" >> ../report/genomes_md5.txt

# -- 3. the settings ---------------------------------------------------------------------------------
# folds supplied by the user (docs/tutorial.md: "--folds folds.json ... or a TSV of sample, fold"),
# here 5 folds alternating within each class
awk -F'\t' 'NR>1 {print $3 "\t" $1}' tutorial/samples.tsv | sort | \
    awk -F'\t' 'BEGIN {print "sample\tfold"} {print $2 "\t" (n[$1]++ % 5)}' > ../folds.tsv
slurm=(-profile hzi,slurm,conda --slurm_partition "$PARTITION" --slurm_extra "$SLURM_EXTRA"
       --max_cpus 8 --max_memory 32.GB --max_time 12.h)
nf() {  # outdir args...: one Nextflow run from the clone, with its default cache and work/, run
    # one after another as a user would -- which is what lets -resume reuse earlier Prokka output
    local out=$1; shift
    nextflow run . -ansi-log false --outdir "$ROOT/results/$out" "$@"
}
step stub-test - nextflow run . -ansi-log false -profile test -stub-run -work-dir work-stub
for preset in hu2024 noleak fast g2p; do
    step "evaluate-$preset" "results/$preset" nf "$preset" "${slurm[@]}" -resume \
        --input tutorial/samples.tsv --preset "$preset" --n_folds 5
done
step evaluate-noleak-folds results/noleak-folds nf noleak-folds "${slurm[@]}" -resume \
    --input tutorial/samples.tsv --preset noleak --folds "$ROOT/folds.tsv"
step train results/model nf model "${slurm[@]}" -resume --mode train --input tutorial/samples.tsv --preset noleak
step predict results/predict nf predict "${slurm[@]}" --mode predict --input tutorial/new_genomes.tsv \
    --model_dir "$ROOT/results/model/model"
# a single machine, no scheduler (docs/tutorial.md sections 3-6 as written): the local executor
# inside this 8-CPU allocation, with its own work directory so nothing is reused
step evaluate-hu2024-local results/hu2024-local \
    nextflow run . -ansi-log false -profile conda -work-dir work-local --outdir "$ROOT/results/hu2024-local" \
    --input tutorial/samples.tsv --preset hu2024 --n_folds 5 --max_cpus 8 --max_memory 30.GB

# -- 4. collect ------------------------------------------------------------------------------------------
cd "$ROOT"
printf 'file\trun\tkind\trows\n' > report/manifest.tsv
for d in results/*/; do
    r=$(basename "$d")
    for f in "$d"evaluation/results/*/summary.tsv "$d"evaluation/results/*/folds.tsv \
             "$d"evaluation/results/*/predictions.tsv "$d"predictions/*.tsv "$d"pipeline_info/trace.tsv; do
        [ -s "$f" ] || continue
        k=$(basename "$f" .tsv); dst="report/${r}__${k}.tsv"; cp "$f" "$dst"
        printf '%s\t%s\t%s\t%s\n' "$(basename "$dst")" "$r" "$k" "$(($(wc -l < "$f") - 1))" >> report/manifest.tsv
    done
done
[ -s results/predict/predictions/vancomycin.tsv ] && \
    join -t $'\t' <(tail -n +2 results/predict/predictions/vancomycin.tsv | sort) \
                  <(tail -n +2 seq2geno2pheno/tutorial/new_genomes_truth.tsv | sort) > report/predict_vs_truth.tsv
log "done"; column -t report/runs.tsv
