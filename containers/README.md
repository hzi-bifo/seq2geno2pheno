# Images

| Image | Used by | Contents | State |
|---|---|---|---|
| `seq2geno2pheno` (`s2g2p/Dockerfile`) | KMERS, GPA, COMPRESS, SPLIT, EVALUATE | s2g2p and its Python dependencies | builds; not published |
| `seq2geno2pheno-tools` (`tools/Dockerfile`) | ROARY, INDELS, PLACE, TRAIN, PREDICT | Seq2Geno's Roary fork with its own BLAST 2.4.0, mafft 7.310, msa2vcf 1.0, s2g2p | **draft** |
| Prokka 1.14.6 | PROKKA | BioContainers | tag to verify |

The `hzi` profile expects Apptainer images under `containers/images/` (gitignored). Build them
in a SLURM job, never on a login node.

Fidelity of the `tools` image is a correctness question, not a convenience: Roary's BLAST is
the step that silently failed in the runs behind the published numbers. Before an image is
used for a reported result, check that:

1. `blastp -version` reports 2.4.0 and starts;
2. Roary on the *E. faecium* cohort reproduces the reproduction's families as partitions, within
   the ~1% run-to-run noise.
