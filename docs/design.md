# Design notes

## Why v3

The 2021 release could not be run as published by 2026:

- Seq2Geno's conda channels (`defaults`, `conda-forge/label/broken`) are unreachable or
  deprecated.
- Its per-rule environments carry absolute paths.
- Roary's bundled BLAST 2.4.0 fails on `libnsl.so.1` on current systems, and Snakemake
  completes anyway. The pangenome is then built without BLAST and no error is reported.
- Modelling depended on a Geno2Pheno web server that is no longer reachable.

Hu et al. 2024 ran Seq2Geno on assemblies, followed by their own local scikit-learn script
(`main_s2p.py`). That script needed patches to run (KMC was only `echo`ed, among others).

v3 keeps the method and replaces the plumbing.

## Split between workflow and package

- **Workflow (Nextflow)** runs the external tools, each in a pinned container: Prokka, Roary,
  mafft, msa2vcf and BLAST.
- **Package (`s2g2p`)** is everything that is Python: reading tool outputs, the ported Seq2Geno
  scripts, placement, and the models.
  - It installs with pip, so the modelling half runs without the workflow.
  - It is unit-tested against the published code.

## What "the published method" is, precisely

From the code, not the paper (amr_benchmark `REPRODUCE_HU2024.md`):

1. **Annotation.** Prokka 1.14.6: `--metagenome --compliant --locustag <id> --prefix <id>`.
2. **Pangenome.** Roary 3.8.2 from Seq2Geno's fork, hzi-bifo/Roary @ `7402d1b`. The fork adds
   a coverage filter to the BLAST edges: identity > 95% and coverage > 95% of both sequences.
   Run as `-e --mafft -g 700000 -z`.
3. **GPA.** Every Roary family, 1 if the genome has a member.
4. **Indels.** Gene families present in more than 2 isolates are aligned and called as follows:
   1. mafft aligns the family's nucleotide sequences;
   2. `msa2vcf` turns the alignment into variant calls;
   3. `vcf2indel` sets a genome to 0 if any copy carries an indel of 8 bp or more relative to
      the majority. A genome without the family is also 0.

   A family becomes a column only if some genome is 0. Three quirks decide *which* families,
   all reproduced under `indel_legacy_rules` and confirmed against the reproduction's own
   `indel.mat` on seven species:
   - **The `head -n` subset.** Seq2Geno aligns the first N lines of `clustered_proteins`, N
     being the number of families in more than 2 isolates, not those N families. A rare family
     that sorts early is aligned; a common one that sorts late is not.
   - **Columns come from the `.Rtab`.** The rule declares the core list as an input but reads
     every family from `gene_presence_absence.Rtab`. A family gets a column when an indel file
     of exactly its name exists. So a rare but aligned family can get one: *S. aureus*
     `group_650` (2 isolates) does.
   - **Non-word characters.** File names replace every non-word character with `_`.
     - Colliding families share one alignment: *S. aureus* `blaR1-2` and `blaR1_2` are aligned
       together.
     - A family whose name has such a character never finds a file of its own name.
5. **Compression.** featCompress keeps one column per distinct pattern, in sorted pattern
   order, named `<first>|<n>`.
6. **6-mers.** KMC canonical counts per record. They are StandardScaler-transformed over every
   genome of the (species, antibiotic) pair.
7. **Model.** For each of `svm lr rf lsvm`:
   - a grid search scored by macro-F1, over the inner splits "each of the 9 other outer folds
     held out once";
   - a refit on the 9 folds;
   - a prediction for the test fold.
8. **Selection.** Per outer fold, the classifier with the highest inner mean. The code sorts by
   (−mean, std) but then looks the winning mean up with `MEAN.index`, so ties go to the first
   classifier and the std tiebreak never fires. `tie_break: first` reproduces this.
9. **Aggregation.** Mean and std over folds. KMA folds use the support-weighted versions.

## What leaks, and what `noleak` does

No step sees test labels. Four unsupervised steps are fitted on all genomes, test genomes
included: the Roary families, the indel alignments and majority calls, the compression, and the
6-mer scaling.

`noleak` fixes each of them:

- **Pangenome and indels.** Built per outer fold, from the cohort minus that fold. Genomes
  unlabelled for this phenotype stay in, matching the published per-species scope.
- **Held-out genomes** are placed one at a time (`s2g2p.place`):
  - each protein joins its best training hit under Roary's own edge criterion;
  - its copies of an indel family are `mafft --add`-ed to the training alignment and called
    with that genome as the only non-training sequence.
- **Scaling and compression** are pipeline steps, so they are refitted on every inner and outer
  training split.

One thing stays shared: the inner folds are scored against families built from all 9 training
folds. This affects model selection only, never the outer estimate. Removing it would cost 9
more pangenomes per outer fold.

The same placement is what makes `--mode predict` possible. A final model's pangenome is built
from all labelled genomes, and a new genome is placed into it.

## Noise floor

Roary is not deterministic. Two runs on identical *E. faecium* input share 98.9% of family
member sets, and group names are renumbered. Compare pangenomes as partitions, and treat any
preset difference smaller than that run-to-run difference as noise.
