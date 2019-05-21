
# GENO2PHENO PIPELINE

## Summary
Geno2Pheno is the second part of Seq2Geno2Pheno bioinformatics software, which is developped for phenotype prediction and characterization using the genotype information. When the Genotype tables already exist it Geno2Pheno can be used as a stand-alone package.
Geno2Pheno has two main functionalities: (i) predictive models of the phenotype (ii) marker detection based on the learned predictive models. The features and functionalities of Geno2Pheno 1.0.0 are detailed in the next sections.

Geno2Pheno works as a data analaysis pipeline that can be customized by the users using a markup format serving as the input of Geno2Pheno, called GeneML. Here are the information will be provided to the Geno2Pheno using GeneML:

![Geno2Pheno](https://user-images.githubusercontent.com/8551117/55457427-243a0080-55ea-11e9-9a55-ce057a2b4b7e.png)

You may see examples of GeneML in the GeneML directory. Seq2Geno automatically generates a GeneML file which can then be modified by the user.

### Running the pipeline

Running the pipeline is pretty straightforward. The only information needed are the GeneML file, number of cores to be used, and if the user wants to override the existing files. Here is an example:

```
python3 geno2pheno.py --genoparse GeneML/pseudomonas.genml --override 0 --cores 30
```


## GeneML and Geno2Pheno functionalities

To make modification of the Geno2Pheno's parameters easier for the users we designed a specific markup input called GeneML.
Here the functionalities together with their related GeneML tags are detailed.

### Input/Output and preprocessings

The main input to the Geno2Pheno are the genotype tables, phylogenetic tree, phenotype table, and output directory.

<ul>
<li> Specifying the output directory and the project name</li>
```
<project output="/net/sgi/metagenomics/projects/pseudo_genomics/results/geno2pheno_package/Pseudogenomics/" name="pseudogenomics_infection">
```

<li>Performing normalizations on top of Genotype tables</li>
Sever

<li>Automatic combining of genotypes</li>
<li>Supporting mutliple phenotypes at once</li>
<li>Generation of the table of results for different classifiers and features</li>
</ul>



### Classification of Phenotype


Geno2Pheno supports:
<ul>
<li>Machine learning classifiers: SVM, Logistic Regression, Random Forests, and Neural Networks.</li>
<li>Randomized and phylogenetic-based K-fold cross validation for hyperparameter tuning for a specified target metric</li>
<li>Validation on a separate held-out dataset with specified size</li>
<li>Allowing user to perform modify phentype mapping</li>
<li>Generation of graph of results</li>
<li>Generation of detailed pickle file for each case</li>
</ul>


### Marker detection

Geno2Pheno selects top Genotypes based on avergae reciprocal rankings of top features identified by the several predictive models and outputs
a final list of important features for further investigations. These markers can be basis for a clinical assays.



# Installation


## (Step 1) Miniconda installation

The first step is to install the latest version of conda on your home directory.

```
cd ~
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Then run the installation and follow the steps:
```
bash ./Miniconda3-latest-Linux-x86_64.sh
```

Then you need to add the conda to your path. Please replace the path with the path you are provided in the installation:
```
export PATH="/mounts/YOUR_MINI_CONDA_PATH/bin:$PATH"
```

Then you need to add conda channels:

```
conda config --add channels conda-forge
```


## (Step 2) Installation of dependencies in the virtual environment

The next step would be installation of the dependencies:

```
conda create --name geno2pheno --file installation/requirements.txt
```

Then you need to activate the created virtual environment:

```
source activate geno2pheno
```



