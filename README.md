<!--
SPDX-FileCopyrightText: 2021 Ehsaneddin Asgari and Tzu-Hao Kuo

SPDX-License-Identifier: GPL-3.0-or-later
-->

# Seq2Geno2Pheno
---
### Contents
- [What is Seq2Geno2Pheno?](#what_is)
- [Recommended workflow](#workflow)
- [License](#license)
- [Download](#download)
- [Contact](#contact)
- [Citation](#citation)
---

### <a name='what_is'></a>What is Seq2Geno2Pheno?

The study of microbial genotype-phenotype relationships plays a key role in many biological contexts, including environmental health, agriculture, antimicrobial resistance (AMR), and infectious diseases diagnostics for precision-medicine. The complex molecular mechanisms underlying these phenotypes motivate developing domain-specific data analysis frameworks. Machine-learning approaches together with domain-knowledge can efficiently translate the ever-increasing amount of sequence data into biological insight. However, the complexity of existing computational workflows and the diverse choices of algorithms and hyperparameters pose a challenge for attaining reproducible and robust results.

Seq2Geno2Pheno provides a comprehensive and configurable workflow for inferring reproducible and robust genotype-phenotype associations. The software includes the Seq2Geno (desktop-app) and Geno2Pheno (web-server) packages running jointly or independently using user-defined configuration files. Seq2Geno computes genomic features from microbial sequence data, from which Geno2Pheno then performs predictive-modeling, evaluates those predictors, and suggests biomarkers.

###  <a name='workflow'></a>Recommended utility 
- Include `--to_gp` when using Seq2Geno to compute the features (more details in the sub-repository), which will automatically
   submit the results to the Geno2Pheno server
- Use the client tool of Geno2Pheno to validate and pack pre-computed data and
  then manually submit the packed materials to the [Geno2Pheno server](http://genopheno.bifo.helmholtz-hzi.de)

###  <a name='download'></a>Download

Please clone this repositroy __recursively__ like:
```
git clone --recurse-submodules https://github.com/hzi-bifo/seq2geno2pheno.git
cd seq2geno2pheno
```

### <a name='license'></a>License

Please refer to LICENSE. All figures in this repository is protected under
[Creative Commons Attribution 3.0 Unported License](https://creativecommons.org/licenses/by/3.0/deed.en_US)

### <a name='contact'></a>Contact

(Note: Please remember to state how your problem can be reproduced and, if accessible, what solutions have been tried)

method 1. Open an issue in this repository

method 2. Send an email to 

- Tzu-Hao Kuo (Tzu-Hao.Kuo@helmholtz-hzi.de): questions about Seq2Geno

- Dr. Ehsaneddin Asgari (asgari@berkeley.edu): questions about Geno2Pheno

### <a name="citation"></a>Citation
We will be publishing the paper for the joint work of Seq2Geno and Geno2Pheno.
Before that, please use 

```
Asgari, E., Kuo, T.-H., Bremges, A., Robertson, G., Weimann, A. & McHardy, A. C. (2021). Seq2Geno2Pheno [A Computational Workflow for Phenotype Predictive-Modeling and Biomarker Detection from Microbial Sequence Data]
```
or 
```
@software{seq2geno2pheno2021,
  author = {Ehsaneddin Asgari, Tzu-Hao Kuo, Andreas Bremges, Gary Robertson, Aaron Weimann, Alice C. McHardy},
  title = {Seq2Geno2Pheno: A Computational Workflow for Phenotype Predictive-Modeling and Biomarker Detection from Microbial Sequence Data},
  version = {v2.00001},
  date = {2021-07-07},
}
```
