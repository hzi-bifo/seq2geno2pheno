#!/usr/bin/env python3

def run_pipeline(config_f):
    import snakemake
    import os
    sgp_home=os.path.dirname(os.path.realpath(__file__))
    snakemake.snakemake(
            snakefile= os.path.join(sgp_home, 'main.smk'),
            use_conda=True,
            conda_prefix= os.path.join(sgp_home, 'env'),
            restart_times= 3, 
            configfile=config_f,
            workdir= os.path.dirname(config_f),
            printshellcmds= True,
            force_incomplete= True,
            notemp=True
            )

def main():
    import argparse
    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='Seq2Geno2Pheno: the pipline tool '
                'for genomic features computation and phenotype predictor '
                'training')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')
    parser.add_argument('-f', dest= 'config_f', required= True, 
        help= 'the yaml file where the arguments are listed')
    parser.add_argument('-d', dest= 'dryrun', required= True, 
        help= 'the yaml file where the arguments are listed')
    args= parser.parse_args()
    run_pipeline(args.config_f)

if __name__== '__main__':
    main()
