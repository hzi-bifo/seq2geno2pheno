# Changelog

## 3.0.0 (unreleased)

Rewrite. v2.00001 is kept at the tag `v2.00001` and the branch `legacy-v2`.

- One workflow (Nextflow) from assemblies to evaluated models, each tool in a pinned container.
- Modelling runs locally (`s2g2p`), replacing the Geno2Pheno web server.
- Presets: `noleak` (default; no test-set information in any feature), `hu2024` (the published
  configuration), `fast`, `g2p` (experimental).
- Train a final model and predict new genomes (`--mode train`, `--mode predict`).
- Seq2Geno's Python 2 helpers ported to Python 3 and tested against the originals.
- Git submodules `seq2geno` and `Geno2PhenoClient` removed from `main`.

## 2.00001 (2021-07-07)

See `docs/legacy/CHANGELOG-v2`.
