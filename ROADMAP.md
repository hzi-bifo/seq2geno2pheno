# Roadmap to v3.0.0

Nothing is published. Status as of 2026-10-12.

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
   - **Models: done for one pair.** `s2g2p evaluate --preset hu2024` reproduces the published
     final macro-F1 exactly on *E. faecium* vancomycin under random, phylotree and KMA folds.
     35-39 of 40 per-classifier scores match; the rest differ through scikit-learn 1.6.1 vs
     0.24.1 (`docs/validation/README.md`).
   - **To do:** a `hu2024` run environment pinned to scikit-learn 0.24.1, so selection cannot
     flip on pairs where RF or LR compete. Then repeat the model check on a species where RF
     wins.
   - **To do:** `--preset noleak` vs the `hu2024-variants/noleak` run.
2. **Installation without containers** (`envs/install.sh`, `-profile conda`, `--mode install`).
   - Lock files in `envs/`: the reproduction's exact lists for Roary's Perl, mafft, msa2vcf and
     the models; Prokka re-locked on conda-forge + bioconda at the same tool versions (the
     reproduction's list took 60 system packages from Anaconda's channel); Roary's programs
     from Bioconda at the fork's versions; BLAST+ 2.4.0, GNU parallel and the Roary fork as
     checksummed downloads.
   - **To do:** generate the new locks and checksums (`dev/lock-envs.sbatch`); install on a
     compute node; Prokka output identical to the container's on the same genomes; the
     tutorial end to end with `-profile conda` (`docs/tutorial.md`), and its pangenome vs the
     container run as partitions.
3. **Images (optional): built on the HZI clusters** (`containers/build-hzi.sbatch`).
   - The legacy envs come from the reproduction's exact package lists.
   - Roary's Perl modules are at the versions the reproduction loaded.
   - Each image's `%test` checks BLAST 2.4.0, mafft 7.310, msa2vcf, Roary and `ps`.
   - A 3-genome Roary run inside the image produced BLAST results.
   - The definitions now read the lock files in `envs/`; **rebuild** before using an image
     again.
4. **End-to-end runs.**
   - **`hu2024`, 30 *E. faecium* genomes: succeeded.** Roary 32 min, indels 9 min; the k-mer
     table equals KMC for all 30 genomes.
   - **`noleak`, same 30: succeeded.**
     - One pangenome per fold from 20 training genomes; Roary 18-21 min per fold.
     - Placement put 10 held-out genomes per fold: 93-97% of their proteins found a family on
       average. One genome dropped to 70%, plausible with a 20-genome reference.
     - 2 of 30 predictions differ from `hu2024`.
   - **`hu2024`, all 277 genomes, published random folds: done** (containers; Roary 9.5 h,
     indels 1.8 h on 16 CPUs). Macro-F1 0.98122 vs the published 0.97377. 4,642 vs 4,637 GPA
     and 631 vs 629 indel columns. Of 277 held-out predictions, 2 differ from the published
     features' (both resistant genomes the published features missed); linear SVM is selected
     in every fold in both. **To do:** the partition comparison with the reproduction's Roary
     (`dev/compare_pangenome.py`, job queued) and the re-evaluation in scikit-learn 0.24.1
     (resume queued).
   - Problems these runs found and fixed:
     - images lacked `ps`;
     - `$SHELL` was the host's zsh;
     - assemblies were not staged into `KMERS`;
     - JVM SIGBUS on a full `/tmp`;
     - Roary's `bin/` was not on `PATH`;
     - concurrent runs shared one Nextflow cache (now `NXF_CACHE_DIR` per output).
5. **The `g2p` grid.** Recover Geno2Pheno's server-side grids (ask E. Asgari) or drop the
   preset from v3.0.
6. **Licence: Academic Software License v1.0** (ASL, as ACE.jl; decided 2026-10-12; v2 was
   GPL-3.0-or-later). **Blockers:**
   - The ported Seq2Geno code (`src/s2g2p/features/`) is © T.-H. Kuo and A. Weimann under
     GPL-3.0-or-later, and the ASL itself says it cannot be combined with GPL code for
     distribution. Kuo and Weimann (and HZI, if it holds the rights) must agree in writing to
     its relicensing under the ASL; record that in `NOTICE`.
   - Confirm the licensor and contact in `NOTICE` (for commercial licences, possibly HZI's
     technology transfer rather than the group).
   - ASL is not an OSI licence: check Bioconda's acceptance before planning a recipe.
7. **Authors:** the v2 list, unchanged (`pyproject.toml`, `CITATION.cff`). Zenodo DOI at release.
8. **Versioning:** `docs/RELEASING.md`; one version in four files, checked by
   `tests/test_version.py`.
9. **CI:** `.github/workflows/ci.yml` (pytest on Python 3.10-3.13, Nextflow stub runs) is
   written but has never run (nothing is pushed).

## After 3.0

- Phylogeny-aware and homology-aware fold construction for users' own cohorts (today: fixed
  fold files, or stratified random).
- A `modern` preset (Bakta + Panaroo), clearly marked as a non-reproduction.
- nf-core-style `nextflow_schema.json`, `nf-test` tests, a MultiQC-style report.
- A Bioconda recipe for the `s2g2p` package, and PyPI.
- skops instead of pickle for model bundles.
- Record convergence warnings per candidate. The protocol rerun found 18% of folds selected a
  non-converged fit.
- Reads input (Seq2Geno's SNP and expression branches): not tested in this project; out of scope
  until it is.
