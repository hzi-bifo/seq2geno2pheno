# Roadmap to v3.0.0

Nothing is published. Status as of 2026-10-09.

## Done (in the tree)

- Repository: v2 history kept, `v2.00001` tag and `legacy-v2` branch; v3 on `main`; the push
  URL is disabled.
- `s2g2p` package:
  - features: k-mers, GPA, compression, indels;
  - placement of held-out genomes;
  - nested CV with the published grids and selection rule;
  - train/predict;
  - presets and CLI.
- Unit tests compare each port with a near-verbatim copy of the published code: `featCompress`,
  `vcf2indel`, `generate_indel_features`, the k-mer vocabulary, and the `main_s2p` +
  `result_analysis` procedure.
- Nextflow workflow covering the evaluate, train and predict modes, plus a stub test profile.

## Before a release candidate

1. **Equivalence on real data.**
   - **Feature ports: done.** `dev/validate_reproduction.py` on 7 species, 277 to 4,757
     genomes, 49/49 checks identical (`docs/validation/`):
     - k-mers vs the KMC dumps;
     - `gpa.mat`;
     - both featCompress tables and their groups;
     - the aligned families;
     - every `vcf2indel` call;
     - `indel.mat`.
   - **To do:** `s2g2p evaluate --preset hu2024` vs the published per-classifier, per-fold scores
     (`dev/compare_hu2024_ml.py`, *E. faecium* vancomycin, 3 schemes, submitted).
   - **To do:** `--preset noleak` vs the `hu2024-variants/noleak` run.
2. **Images.**
   - `seq2geno2pheno-tools`: hzi-bifo/Roary @ `7402d1b` with its `install_dependencies.sh`
     (BLAST 2.4.0, CD-HIT 4.6.6, MCL 14-137, PRANK, bedtools 2.26), the `perl5_22_env` pins,
     mafft 7.310, jvarkit-msa2vcf 1.0, and s2g2p.
   - `seq2geno2pheno`: s2g2p only.
   - Verify the Prokka 1.14.6 BioContainers tag.
   - Build as Apptainer images on the cluster (SLURM job), then Docker for CI.
3. **End-to-end run** on a small real cohort (one species, a few hundred genomes) for every
   preset, on SLURM.
4. **The `g2p` grid.** Recover Geno2Pheno's server-side grids (ask E. Asgari) or drop the
   preset from v3.0.
5. **Authors, name and licence.** Confirm with T.-H. Kuo, E. Asgari and A. C. McHardy. Then add
   `CITATION.cff` and a Zenodo DOI.

## After 3.0

- Phylogeny-aware and homology-aware fold construction for users' own cohorts (today: fixed
  fold files, or stratified random).
- A `modern` preset (Bakta + Panaroo), clearly marked as a non-reproduction.
- nf-core-style `nextflow_schema.json`, `nf-test` tests, GitHub Actions CI, a MultiQC-style
  report.
- skops instead of pickle for model bundles.
- Record convergence warnings per candidate. The protocol rerun found 18% of folds selected a
  non-converged fit.
- Reads input (Seq2Geno's SNP and expression branches): not tested in this project; out of scope
  until it is.
