# Seq2Geno2Pheno
---
### Contents
- [What is Seq2Geno2Pheno?](#what_is)
- [Recommended workflow](#workflow)
- [License](#license)
- [Download](#download)
- [Contact](#contact)
---

### <a name='what_is'></a>What is Seq2Geno2Pheno?

The study of microbial genotype-phenotype relationships plays a key role in many biological contexts, including environmental health, agriculture, antimicrobial resistance (AMR), and infectious diseases diagnostics for precision-medicine. The complex molecular mechanisms underlying these phenotypes motivate developing domain-specific data analysis frameworks. Machine-learning approaches together with domain-knowledge can efficiently translate the ever-increasing amount of sequence data into biological insight. However, the complexity of existing computational workflows and the diverse choices of algorithms and hyperparameters pose a challenge for attaining reproducible and robust results.

Seq2Geno2Pheno provides a comprehensive and configurable workflow for inferring reproducible and robust genotype-phenotype associations. The software includes the Seq2Geno (desktop-app) and Geno2Pheno (web-server) packages running jointly or independently using user-defined configuration files. Seq2Geno computes genomic features from microbial sequence data, from which Geno2Pheno then performs predictive-modeling, evaluates those predictors, and suggests biomarkers.

###  <a name='workflow'></a>Recommended workflow
1. Use [seq2geno](seq2geno/) to compute the features 
2. Use the generator of genyml file (seq2geno/submission_tool/create_genyml.py) to make the input file for Geno2Pheno. The generated file can be further manually updated.
3. Use the [data validator](Geno2PhenoClient/) for Geno2Pheno server to validate and compress data into a zip file
4. Submit the zip file to [Geno2Pheno server](http://genopheno.bifo.helmholtz-hzi.de)

###  <a name='download'></a>Download

Please clone this repositroy __recursively__ like:
```
git clone --recurse-submodules https://github.com/hzi-bifo/seq2geno2pheno.git
cd seq2geno2pheno
git submodule update --init --recursive
```

### <a name='license'></a>License

Please refer to LICENSE

### <a name='contact'></a>Contact

(Note: Please remember to state how your problem can be reproduced and, if accessible, what solutions have been tried)

method 1. Open an issue in this repository

method 2. Send email to 

- Tzu-Hao Kuo (Tzu-Hao.Kuo@helmholtz-hzi.de): questions about Seq2Geno

- Ehsaneddin Asgari (asgari@berkeley.edu): questions about Geno2Pheno


