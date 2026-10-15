# Images (optional)

The supported way to install is without containers: `envs/install.sh` and `-profile conda`
(README, docs/tutorial.md). The images here are built from the same lock files in `envs/`, for
sites that prefer containers. Images built before the lock files moved to `envs/` used the
reproduction's Prokka list (`docs/legacy/prokka_env.reproduction.lock.txt`).

| Image | Definition | Used by | Contents |
|---|---|---|---|
| `s2g2p.sif` | `s2g2p/s2g2p.def` (`s2g2p/Dockerfile` for Docker) | KMERS, GPA, COMPRESS, SPLIT, EVALUATE | s2g2p and its Python dependencies |
| `prokka-1.14.6.sif` | `prokka/prokka.def` | PROKKA | `envs/prokka.lock.txt`: the tool versions of Seq2Geno's `prokka_env` as the reproduction ran it |
| `s2g2p-tools.sif` | `tools/tools.def` | ROARY, INDELS, PLACE, TRAIN, PREDICT | Seq2Geno's Roary fork (hzi-bifo/Roary @ 7402d1b) with its own BLAST 2.4.0; mafft 7.310; msa2vcf 1.0; s2g2p |
| `s2g2p-hu2024ml.sif` | `hu2024ml/hu2024ml.def` | EVALUATE with `hu2024` | the published model environment (scikit-learn 0.24.1) and s2g2p |

## Why exact package lists

Every conda environment in the legacy images is created from `conda list --explicit --md5`
of the environment Seq2Geno built and ran in the Hu et al. 2024 reproduction
(`envs/*.lock.txt`). That reproduction's tables are what `dev/validate_reproduction.py` checks the ports
against, so these are the versions known to produce them. Re-solving the 2016-era pins instead
gives whatever the channels hold today.

Prokka is pinned with its whole stack because gene names become Roary family names, and family
names decide which indel columns exist under the legacy rules (`docs/design.md`).

The exception is Prokka: 60 packages of the reproduction's `prokka_env` (system libraries,
Perl, Java) came from Anaconda's own channel. `envs/prokka.lock.txt` keeps every Bioconda tool
at the same version and takes the rest from conda-forge (`envs/specs/prokka.txt`,
`dev/lock-envs.sbatch`).

## Building on the HZI clusters

```bash
sbatch -p cpu -A bifo -o log/%x_%j.out containers/build-hzi.sbatch            # all four
sbatch -p cpu -A bifo -o log/%x_%j.out containers/build-hzi.sbatch tools      # one
```

- Builds run on compute nodes. There is no subuid mapping, so Apptainer uses a root-mapped
  namespace, and `dnf` only works there with `--ignore-fakeroot-command`.
- Each image gets its own cache directory. Two builds unpacking the same base layer into one
  NFS cache fail with "stale file handle".
- Images land in `containers/images/` (gitignored), which the `hzi` profile reads.

## Fidelity checks

Roary's BLAST is the step that silently failed in the runs behind the published numbers. Before
an image is used for a reported result:

1. The image's `%test` must pass: BLAST 2.4.0 and mafft 7.310 start, msa2vcf runs, and Roary
   reports 3.8. The build also fails if any Roary dependency or any of its 13 Perl modules is
   missing.
2. Roary on the *E. faecium* cohort must reproduce the reproduction's families as partitions,
   within the ~1% run-to-run noise.
