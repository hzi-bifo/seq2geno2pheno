# Seq2Geno2Pheno
---
### Contents
- [What is Seq2Geno2Pheno?](#what_is)
- [What can Seq2Geno2Pheno do?](#what_can)
- [Get started](#get_started)
	- [Prerequisites](#prerequisites)
	- [Installation](#installation)
- [Usage](#usage)
- [Example dataset](#example)
- [Input yaml file](#input)
- [License](#license)
- [Contact](#contact)
---
### <a name='what_is'></a>What is Seq2Geno2Pheno?
This is a computational framework for genomic studies of bacterial population. Compared to conventional methods, such as shell scripts, this package emphasizes data reproducibility and computational environment management. 

### <a name='what_can'></a>What can Seq2Geno2Pheno do?
Seq2Geno2Pheno includes two main stages. The first stage Seq2Geno aims to compute genomic features with sequencing reads. Briefly speaking, it covers (1) variant calling, (2) expression analysis, (3) _de novo_ assembly and gene content identification, and (4) phylogeny inference. These are followed by optional functionals of (5) differential expression analysis and (6) ancestral reconstruction. 

The outputs from Seq2Geno are formatted for the subsequent stage: Geno2Pheno. This stage mainly trains phenotype predictors with the genomic data. Furthermore, it reports lists of genomic factors that are potentially linked to the target phenotype using feature selection techniques. For standalone running of Geno2Pheno please refer to the documentation in the Geno2Pheno directory.

![SEQ2GENO2PHENO-Copy of Overview diagram](https://user-images.githubusercontent.com/8551117/59326225-24d7cd00-8ce6-11e9-86ef-137856f1714c.png)


### <a name='get_started'></a>Get started
- <a name='prerequisites'></a>Requirements

    - [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) (tested version: 4.6.14)
    - file [.condarc](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html) that includes these channels and is detectable by your conda
		- conda-forge/label/broken
		- bioconda
		- conda-forge
		- defaults
    - [python](https://www.python.org/downloads/) (tested verson: 3.6)
    - [Linux](https://www.cyberciti.biz/faq/find-linux-distribution-name-version-number/) (tested version: Debian GNU/Linux 8.8 jessie)
    - [git](https://git-scm.com/downloads) (tested version: 2.18)

- <a name='installation'></a>Installation

    1. Download Seq2Geno2Pheno:

	```
	git clone --recurse-submodules https://github.com/hzi-bifo/seq2geno2pheno.git
	cd seq2geno2pheno/
	git submodule update --init --recursive
	```

	The option `--recurse-submodules` helps to download the submodules that are located at another repository (i.e. Seq2Geno and Geno2Pheno). The flag is available only in git version >2.13, and users of earlier git versions may need to find the substitute.  

    2. Install the core environment:

	Please refer to `install/INSTALL.md`. 	

	The core environment means the set of computational tools that are required to initiate Seq2Geno2Pheno. The environment, like those in the other parts of this package, is stated in a file that is recognizable by Conda, which makes the installation of everything achievable in one command. 	

    3. Finally, set the variable `PATH` and `SGP_HOME`. You may either manually set them up using:

	```
	export SGP_HOME=/where/the/seq2geno2pheno/is/at
	export PATH=$SGP_HOME:$PATH
	```

	or insert the above line in the `.profile` or `.bashrc` of your home directory.

	Until here, the core environment of Seq2Geno2Pheno should have been installed. However, it still lacks the process-specific environments, which are computational environments required by only part of the workflows. These environments will be automatically installed when for the first time called by a scheduled process. Therefore, we encourage every user to run Seq2Geno2Pheno with a small dataset, such as our example dataset, to install these environments and facilitate future usages.

### <a name='usage'></a>Usage

```
usage: sgp [-h] [-v] -f CONFIG_F

Seq2Geno2Pheno: the pipline tool for genomic features computation and phenotype predictor training

optional arguments:
  -h, --help
    show this help message and exit
  -v
    show program's version number and exit
  -f CONFIG_F
    the yaml file where the arguments are listed
```

Example:
```
   source activate sgp_env
   sgp -f options.yml
   source deactivate
```

The only input file (i.e. `options.yml`) is a yaml file where all options are described. We also offer an example dataset where reads and commands are included (see below).

### <a name='example'></a>Example dataset

We offer an example dataset `example_sgp_dataset.tar.gz`. Please decompress it by 

```
tar xzvf example_sgp_dataset.tar.gz
cd example_sgp_dataset
```
, and then refer to the tutorial `USAGE.md`. 

### <a name='input'></a>Input yaml file

(Note: in __general__ and __prediction__ session, the file or directory paths must be __absolute path__)

(Note: Geno2Pheno requires at least ~ 100 samples, with a reasonable balance in the positive and negative labels in order to perform properly) 

1. config_f: the input filename itself

2. features: target computational processes

| option | action | values ([default])|
| --- | --- | --- |
| dryrun | display the processes and exit | [Y]/N |
| s | identify variants (SNPs) | Y/[N] |
| d | create _de novo_ assemblies and analyze gene contents | Y/[N] |
| e | count expression levels | Y/[N] |
| p | infer the phylogeny | Y/[N] |
| de | conduct analysis of differential expression | Y/[N] |
| ar | ancestrally reconstruct the expression levels | Y/[N] |

3. general: materials and resources

    - cores: available number of CPUs 

    - old_config: do not overwrite the previous config files of seq2geno. (Please make sure the old config file really are accessible.)

    - sgp_output_d: the working directory

    The intermediate and final files will be created under the folder. 

    - dna_reads: The list of DNA-seq data 

    It should be a two-column list, where the first column includes all samples and the second column lists the __paired-end reads files__. The two reads file are separated by a comma. The first line is the first sample.
    ```
    sample01	/paired/end/reads/sample01_1.fastq.gz,/paired/end/reads/sample01_2.fastq.gz
    sample02	/paired/end/reads/sample02_1.fastq.gz,/paired/end/reads/sample02_2.fastq.gz
    sample03	/paired/end/reads/sample03_1.fastq.gz,/paired/end/reads/sample03_2.fastq.gz
    ```

    - rna_reads: The list of RNA-seq data

    It should be a two-column list, where the first column includes all samples and the second column lists the __short reads files__. The first line is the first sample.
    ```
    sample01	/transcription/reads/sample01.rna.fastq.gz
    sample02	/transcription/reads/sample02.rna.fastq.gz
    sample03	/transcription/reads/sample03.rna.fastq.gz
    ```

    - phe_table: The phenotype table

    The table is tab-separated. For n samples with m phenotypes, the table is (n+1)-by-(m+1) as shown below. The first column should list the sample names. The header line is supposed to includes names of phenotypes. Missing values are acceptable.
    ```
    strains	virulence
    sample01	high
    sample02	mediate
    sample03	low
    ```

    - ref_fa, ref_gff, ref_gbk:	The reference data

    The fasta, gff, and genbank files of a reference genome. They should have the same sequence accessions. 

    - adaptor: The adaptor file (optional)

    The fasta file of adaptors of DNA-seq data. Pruning DNA-seq reads will be activated if the file is properly specified.

4. prediction: parameters about machine learning classification and feature selection 

    - pred: a name for the prediction -- this parameter is useful when different classifications schemes are applied on the same data

    - models: the classification method 
		- choices (multiple selection is possible): 
			- SVM (for support vector machine)
			- RF (for random forests)
			- LR (for logistic regression]

    - part: the method to partition samples 
		- choices (each prediction block can have 1 partitioning method): 
			- tree (for phylogenetic based splitting in CV and test set creation)
			- rand (for random splitting)])

    - fold_n: the number of fold for cross-validation (integer; commonly 5 or 10)

    - optimize: the target metric to optimize 
		- choices:
			- 'accuracy': 'accuracy',
			- 'scores_p_1': 'precision of the positive class'
			- 'scores_r_1': 'recall of the positive class'
			- 'scores_p_0': 'precision of the negative class'
			- 'scores_r_0': 'recall of the negative class'
			- 'scores_f1_1': 'f1 of the positive class'
			- 'scores_f1_0': 'f1 of the negative class'
			- 'precision_micro': 'precision micro'
			- 'precision_macro': 'precision macro'
			- 'recall_macro': 'recall macro'
			- 'recall_micro': 'recall micro'
			- 'f1_macro': 'f1_macro' 
			- 'f1_micro': 'f1_micro'

    - test_ratio: the percentage of samples  (float number between 0.0 and 1.0; commonly 0.1)
    The proportion of samples will be isolated from the whole set for independent testing. These samples will not be used to train the predictor. 

    - kmer: the k-mer size for encoding genome sequences

    - classes_f: classification labels (filename)
    The file specifies classification labels based on the phenotypes. It is a two-column tab-separated file, where the first column is the phenotypes and the second column includes their prediction labels.  Label 1 is used as the positive class and 0 as the negative class. For the multiclass classification the optimizing score needs to be chosen from scores that are not specific to the positive nor negative class.
    ```
    high	1
    mediate	1
    low	0
    ```

    - ovrd: whether to overwrite the existing geno2pheno outputs (choices: 1:yes, 0:no)
    

### <a name='license'></a>License
Apache 2.0 (please refer to LICENSE)

### <a name='contact'></a>Contact

(Note: Please remember to state how your problem can be reproduced and, if accessible, what solutions have been tried)

method 1. Open an issue in this repository

method 2. Send email to 

- Tzu-Hao Kuo (Tzu-Hao.Kuo@helmholtz-hzi.de): questions about Seq2Geno

- Ehsaneddin Asgari (asgari@berkeley.edu): questions about Geno2Pheno


