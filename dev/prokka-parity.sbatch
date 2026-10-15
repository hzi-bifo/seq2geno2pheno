#!/usr/bin/env bash
#SBATCH --job-name=s2g2p-prokka-parity
#SBATCH --account=bifo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#
# Does the conda Prokka (envs/prokka.lock.txt) annotate exactly as the reproduction's environment
# did? Re-annotates genomes with the workflow's PROKKA command and compares every output that
# downstream steps read (.faa .ffn .gff, plus .tsv .tbl) byte for byte with a finished run's.
#   sbatch -o log/%x_%j.out dev/prokka-parity.sbatch <sample sheet> <finished run outdir> <N> [env_dir]
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
sheet=$1 run=$2 n=$3 envs=${4:-$PWD/envs/installed}
out=$run/compare/prokka-parity
rm -rf "$out" && mkdir -p "$out"
trace=$run/pipeline_info/trace.tsv
printf 'sample\tfile\tidentical\n' > "$out/result.tsv"
# the first n samples (awk, not head: under pipefail, head's early exit fails the pipe)
picked=$(tail -n +2 "$sheet" | sort | awk -v n="$n" 'NR <= n')
while IFS=$'\t' read -r s fasta _; do
    h=$(awk -F'\t' -v s="PROKKA ($s)" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["name"]==s{print $c["hash"]}' "$trace")
    ref=$(ls -d "$run"/work/"$h"*)/$s
    (cd "$out" && "$envs/bin/prokka" --locustag "$s" --prefix "$s" --force --cpus 8 \
        --metagenome --compliant --outdir "$s" "$fasta" > /dev/null 2>&1)
    for ext in faa ffn gff tsv tbl; do
        if cmp -s "$ref/$s.$ext" "$out/$s/$s.$ext"; then r=yes; else r=no; fi
        printf '%s\t%s\t%s\n' "$s" "$ext" "$r" >> "$out/result.tsv"
    done
done <<< "$picked"
awk -F'\t' 'NR>1{n++; if($3=="yes") y++} END{print y" of "n" files identical"}' "$out/result.tsv"
grep -P '\tno$' "$out/result.tsv" || true
