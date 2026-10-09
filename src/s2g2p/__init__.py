"""Seq2Geno2Pheno v3.

Seq2Geno turns bacterial assemblies into genotype tables (gene presence/absence, indels,
k-mers); Geno2Pheno fits and evaluates phenotype predictors on them. v3 runs both locally, with
named presets that reproduce the published configurations and a default that fits every
feature on training genomes only.
"""

__version__ = "3.0.0.dev0"
