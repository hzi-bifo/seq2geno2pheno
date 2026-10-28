# Changelog

All notable changes; versions follow [docs/RELEASING.md](docs/RELEASING.md) (Semantic
Versioning, judged by what happens to results).

## 3.0.0 (unreleased)

Rewrite. v2.00001 is kept at the tag `v2.00001` and the branch `legacy-v2`.

- Licence changed from GPL-3.0-or-later to the Academic Software License v1.0 (ASL;
  non-commercial academic use). See LICENSE and NOTICE.
- Installed without containers: `envs/install.sh` (or `--mode install`) creates one pinned conda
  environment per tool, and `-profile conda` runs them. Container images remain optional.
- One workflow (Nextflow) from assemblies to evaluated models, each tool at a pinned version.
- Modelling runs locally (`s2g2p`), replacing the Geno2Pheno web server.
- Presets: `noleak` (default; no test-set information in any feature), `hu2024` (the published
  configuration), `fast`, `g2p` (experimental).
- Train a final model and predict new genomes (`--mode train`, `--mode predict`).
- Outer folds for any cohort (`--mode folds`, `s2g2p folds`): stratified random,
  phylogeny-aware (NJ tree leaf order in equal blocks) and homology-aware (KMA, ported from
  Hu et al. 2024). See docs/folds.md.
- A warning when no `-work-dir` is given; the README walks through a run step by step.
- Seq2Geno's Python 2 helpers ported to Python 3 and tested against the originals.
- Git submodules `seq2geno` and `Geno2PhenoClient` removed from `main`.

## 2.00001 (2021-07-07)

See `docs/legacy/CHANGELOG-v2`.
