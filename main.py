#!/usr/bin/env python3

def run_pipeline(config_f):
    import snakemake
    import os
    sgp_home=os.path.dirname(os.path.realpath(__file__))
    try:
        outcome= snakemake.snakemake(
            snakefile= os.path.join(sgp_home, 'main.smk'),
            use_conda=True,
            conda_prefix= os.path.join(sgp_home, 'env'),
            configfile=config_f,
            workdir= os.path.dirname(config_f),
            printshellcmds= True,
            force_incomplete= True,
            notemp=True
            )
        if not outcome:
            raise RuntimeError
    except RuntimeError as e:
        print('\n'*2)
        print('\033[91m' + "Seq2Geno2Pheno was stopped."
            "If this is unexpected, checking files or results under the "
            "output directory may help you to find the cause." + '\033[0m')
        print('\n'*2)
    else:
        print('\n'*2)
        print('\033[91m'
              'Seq2Geno2Pheno complete!\n'
              'Please consider citing Seq2Geno2Pheno if it helped '
              '\033[0m')
        print('\n'*2)

    
def main():
    import argparse
    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='Seq2Geno2Pheno: the pipline tool '
                'for computing genomic features and training phenotype '
                'classifiors')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')
    parser.add_argument('-f', dest= 'config_f', required= True, 
        help= 'the yaml file where the arguments are listed')
    args= parser.parse_args()
    run_pipeline(args.config_f)


if __name__== '__main__':
    main()
