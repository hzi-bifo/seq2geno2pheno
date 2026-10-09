# Images

| Image | Definition | Used by | Contents |
|---|---|---|---|
| `s2g2p.sif` | `s2g2p/s2g2p.def` (`s2g2p/Dockerfile` for Docker) | KMERS, GPA, COMPRESS, SPLIT, EVALUATE | s2g2p and its Python dependencies |
| `prokka-1.14.6.sif` | `prokka/prokka.def` | PROKKA | the exact package list of Seq2Geno's `prokka_env` as the reproduction ran it |
| `s2g2p-tools.sif` | `tools/tools.def` | ROARY, INDELS, PLACE, TRAIN, PREDICT | Seq2Geno's Roary fork (hzi-bifo/Roary @ 7402d1b) with its own BLAST 2.4.0; mafft 7.310; msa2vcf 1.0; s2g2p |

## Why exact package lists

Every conda environment in the legacy images is created from `conda list --explicit --md5`
of the environment Seq2Geno built and ran in the Hu et al. 2024 reproduction (`*.explicit.txt`
here). That reproduction's tables are what `dev/validate_reproduction.py` checks the ports
against, so these are the versions known to produce them. Re-solving the 2016-era pins instead
gives whatever the channels hold today.

Prokka is pinned with its whole stack because gene names become Roary family names, and family
names decide which indel columns exist under the legacy rules (`docs/design.md`).

60 packages of `prokka_env` come from the `anaconda` channel, as in the published
environment. Check Anaconda's terms of service before redistributing that image.

## Building on the HZI clusters

```bash
sbatch -p cpu -A b_algbio -o log/%x_%j.out containers/build-hzi.sbatch            # all three
sbatch -p cpu -A b_algbio -o log/%x_%j.out containers/build-hzi.sbatch tools      # one
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
