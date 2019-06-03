#!/bin/bash

conda env create --file=sgp_env.yml
source activate sgp_env
cd ../seq2geno/denovo/lib/Roary/
./install_dependencies.sh
source deactivate 
