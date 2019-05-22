#!/bin/bash
#$ -V
#$ -l h_core=20
#$ -l h_vmem=200G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/test/v2
source activate snakemake_env

#seq2geno  -f options.yaml 
#snakemake -p --use-conda \
#  --conda-prefix=/net/metagenomics/data/from_moni/old.tzuhao/seq2geno2pheno/env \
#  --configfile=options.sgp.yaml \
#  --snakefile=/net/metagenomics/data/from_moni/old.tzuhao/seq2geno2pheno/main.smk
/net/metagenomics/data/from_moni/old.tzuhao/seq2geno2pheno/main.py -f options.sgp.yaml
