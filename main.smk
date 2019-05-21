
__license__ = "Apache 2"
__version__ = "1.0.0"
__url__='https://github.com/hzi-bifo/seq2geno2pheno.git'

'''
As an interface to link seq2geno, genml creator, and geno2pheno
'''

import os

#sgp_home= os.path.dirname(os.path.realpath(__file__))
sgp_home= '/net/sgi/metagenomics/data/from_moni/old.tzuhao/seq2geno2pheno'
sys.path.append(os.path.join(sgp_home, 'lib'))

seq2geno_home= os.path.join(sgp_home, 'seq2geno')
sys.path.append(seq2geno_home)
os.environ['SEQ2GENO_HOME']=seq2geno_home 
geno2pheno_home= os.path.join(sgp_home, 'geno2pheno')
sys.path.append(geno2pheno_home)
os.environ['PATH']= ':'.join([
    os.path.join(seq2geno_home, 'main'), 
    geno2pheno_home,
    os.environ['PATH']])
    

## accept arguments
import UserOptionsSGP
#args= UserOptionsSGP.main()
dsply_args= (True if config['d'] == 'Y' else False)
opt_f= config['f']
args= UserOptionsSGP.read_options(opt_f, dsply_args)

rule all:
    input:
#        sg_out= os.path.join(args.wd,'RESULTS'),
#        genml= os.path.join(args.sgp_output_d, 'info', 'seq2geno.gml')
#        gp_out= directory(args.gp_output_d)

rule create_options_yaml:
    '''
    Create the input yaml file for seq2geno
    ''' 
    input:
        opt_yml= opt_f
    output:
        target_yml= os.path.join(args.sgp_output_d, 'info', 'seq2geno_inputs.yml')
    run:
        import yaml
        yml_fh=open(input.opt_yml, 'r')
        opt_dict= {}
        with open(input.opt_yml, 'r') as yml_fh:
            opt_dict= yaml.safe_load(yml_fh)
        available_functions= [
            'snps', 'expr', 'denovo', 'phylo', 'de', 'ar', 'dryrun']
        for k in available_functions: 
            if not k in opt_dict['functions']:
                opt_dict['functions'][k]='N'

        target_opt_dict = {x: opt_dict[x] for x in ['functions', 'files']}
        target_opt_dict['files']['wd']= args.wd
        del target_opt_dict['files']['sgp_output_d']
        
        with open(output.target_yml, 'w') as out_fh:
            yaml.dump(target_opt_dict, out_fh, default_flow_style=False)

rule seq2geno:
    input: 
        target_yml= os.path.join(args.sgp_output_d, 'info', 'seq2geno_inputs.yml')
    output: 
        sg_key_output_d= os.path.join(args.wd,'RESULTS')
    conda: os.path.join(seq2geno_home, 'install', 'snakemake_env.yml')
    shell:
        '''
        which seq2geno
        seq2geno -f {input.target_yml}
        '''

rule create_genml:
    input: 
        sg_key_output_d= os.path.join(args.wd,'RESULTS')
    output:
        genml= os.path.join(args.sgp_output_d, 'info', 'seq2geno.gml')
    params:
        sgp_home= sgp_home,
        args= args
    run:
        import sys
        import os
        sgp_home= params.sgp_home 
        sys.path.append(os.path.join(sgp_home, 'lib'))
        ## load genml creator module
        import create_genml_from_seq2geno
        create_genml_from_seq2geno.main(params.args)
        

rule geno2pheno:
    ''' 
    ensure Geno2Pheno is executable
    '''
    input:
        #sg_key_output_d= os.path.join(args.wd,'RESULTS'),
        genml= os.path.join(args.sgp_output_d, 'info', 'seq2geno.gml')
    output: directory(args.gp_output_d)
    conda: os.path.join(geno2pheno_home, 'installation', 'requirements.txt')
    params:
        ovrd= args.ovrd
    threads: args.cores
    shell:
        '''
        Geno2Pheno.py --genoparse {input.genml} \
--override {params.ovrd} \
--cores {threads}
        '''
