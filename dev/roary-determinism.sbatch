#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-roarydet
#SBATCH --account=bifo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-3
#
# Is Roary's run-to-run variation Perl's hash randomisation? Runs Roary (the tools image) on the
# same GFFs four times: tasks 0-1 with PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0, tasks 2-3 with the
# default random seed. Then compare clustered_proteins of 0 vs 1 and 2 vs 3 as partitions
# (dev/compare_roary_stages.py).
#   sbatch -o log/%x_%A_%a.out --exclude=bioinf009 dev/roary-determinism.sbatch dev-runs/roary-determinism
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
out=$1
t=${SLURM_ARRAY_TASK_ID:-0}
run=$out/run$t
rm -rf "$run" && mkdir -p "$run/gff"
# copies, not links: the link targets are not mounted inside the container
while read -r g; do cp "$g" "$run/gff/"; done < "$out/gffs.txt"
seed=()
if [ "$t" -lt 2 ]; then seed=(--env PERL_HASH_SEED=0 --env PERL_PERTURB_KEYS=0); fi
cd "$run"
apptainer exec --env TMPDIR=/tmp --env SHELL=/bin/bash "${seed[@]}" \
    "$SLURM_SUBMIT_DIR/containers/images/s2g2p-tools.sif" \
    bash -c 'echo "PERL_HASH_SEED=${PERL_HASH_SEED:-random}"; roary -f roary -v gff/*.gff -p 8 -e --mafft -g 700000 -z' \
    > roary.log 2>&1
md5sum roary/clustered_proteins roary/_clustered.clstr roary/_inflated_mcl_groups roary/_uninflated_mcl_groups
