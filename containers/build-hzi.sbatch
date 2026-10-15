#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-images
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --nice=30000
#
# Build the Apptainer images into containers/images/ on an HZI compute node (fakeroot works
# there, not on the login node; without a subuid mapping Apptainer uses a root-mapped
# namespace, and dnf only works there with --ignore-fakeroot-command). From the repository root:
#   sbatch -p cpu -A bifo -o log/%x_%j.out containers/build-hzi.sbatch [s2g2p|prokka|tools ...]
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p containers/images
# The build sandbox must be node-local: the nodes' $TMPDIR is on NFS, where rpm's
# unpack-then-rename fails ("cpio: open failed").
scratch=$(mktemp -d /tmp/s2g2p-build.XXXXXX)
trap 'rm -rf "$scratch"' EXIT
export APPTAINER_TMPDIR="$scratch"
declare -A DEF=([s2g2p]=containers/s2g2p/s2g2p.def [prokka]=containers/prokka/prokka.def [tools]=containers/tools/tools.def [hu2024ml]=containers/hu2024ml/hu2024ml.def)
declare -A SIF=([s2g2p]=s2g2p.sif [prokka]=prokka-1.14.6.sif [tools]=s2g2p-tools.sif [hu2024ml]=s2g2p-hu2024ml.sif)
for img in "${@:-s2g2p prokka tools hu2024ml}"; do
  for i in $img; do
    echo "=== $i $(date -Is)"
    # one cache per image: two builds unpacking the same base layer into one NFS cache collide
    export APPTAINER_CACHEDIR="$PWD/containers/images/.cache-$i"
    apptainer build --fakeroot --ignore-fakeroot-command --force "containers/images/${SIF[$i]}.tmp" "${DEF[$i]}"
    mv "containers/images/${SIF[$i]}.tmp" "containers/images/${SIF[$i]}"
    apptainer inspect --deffile "containers/images/${SIF[$i]}" | sha256sum | sed "s/-$/${SIF[$i]} deffile/"
  done
done
