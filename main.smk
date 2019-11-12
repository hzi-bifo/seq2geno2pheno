
__license__ = "Apache 2"
__version__ = "1.0.0"
__url__='https://github.com/hzi-bifo/seq2geno2pheno.git'

'''
As an interface to link seq2geno, genml creator, and geno2pheno
'''

import os

sgp_home=''
try:
    sgp_home= os.environ['SGP_HOME']
except KeyError as e:
    sys.exit('SGP_HOME was inappropriately set. '
    'Please ensure the correct setting or refer to install/INSTALL.md')
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
opt_f= config['config_f']
args= UserOptionsSGP.read_options(opt_f=opt_f)
## determine time code
## for 
import datetime
import random
time= datetime.datetime.now()
time_code='sgp{}_{}.{}_{}'.format(str(time.day), str(time.month),
str(time.hour), str(time.minute))
info_d= os.path.join(args.sgp_output_d, 'info', time_code)
info_d_copy= info_d
n= 0
while os.path.isdir(info_d):
    n=n+1
    info_d= info_d_copy+'_{}'.format(str(n))

target_data=[os.path.join(info_d,'seq2geno_report')]
if not args.dryrun:
    target_data.append(os.path.join(info_d,'geno2pheno_report'))

rule all:
    input:
        target_data

rule create_options_yaml:
    '''
    Create the input file for seq2geno
    ''' 
    input:
        opt_yml= opt_f
    output:
        target_yml= os.path.join('{info_d}', 'seq2geno_inputs.yml')
    threads: 1
    run:
        import yaml
        yml_fh=open(input.opt_yml, 'r')
        opt_dict= {}
        with open(input.opt_yml, 'r') as yml_fh:
            opt_dict= yaml.safe_load(yml_fh)
        available_functions= [
            'snps', 'expr', 'denovo', 'phylo', 'de', 'ar', 'dryrun']
        for k in available_functions: 
            if not k in opt_dict['features']:
                opt_dict['features'][k]='N'

        target_opt_dict = {x: opt_dict[x] for x in ['features', 'general']}
        target_opt_dict['general']['wd']= args.wd
        del target_opt_dict['general']['sgp_output_d']
        
        with open(output.target_yml, 'w') as out_fh:
            yaml.dump(target_opt_dict, out_fh, default_flow_style=False)

rule seq2geno:
    input: 
        target_yml= os.path.join('{info_d}', 'seq2geno_inputs.yml')
    output: 
        sg_check= os.path.join('{info_d}', 'seq2geno_report')
    conda: os.path.join(seq2geno_home, 'install' , 'snakemake_env.yml')
    log: os.path.join('{info_d}', 'seq2geno_screenshot')
    threads: args.cores
    shell:
        '''
        which seq2geno
        seq2geno -f {input.target_yml} &> {log}
        ## Only when SG was sucessfully done, and new or precomputed data
        ## are available, this step is considered successful
        printf "Seq2Geno: $(date)\ndetails recorded in: {log}" \
> {output.sg_check};
        '''

rule create_genml:
    '''
    Create the input file for geno2pheno
    ''' 
    input: 
        sg_check= os.path.join('{info_d}', 'seq2geno_report'),
        opt_yml= opt_f
    output:
        genml= os.path.join('{info_d}', 'seq2geno.gml')
    params:
        sgp_home= sgp_home,
        args=args
    threads: 1
    run:
        import sys
        import os
        sgp_home= params.sgp_home 
        sys.path.append(os.path.join(sgp_home, 'lib'))
        ## load genml creator module
        import create_genml_from_seq2geno
        create_genml_from_seq2geno.main(params.args, genml_f= output.genml)
        

rule geno2pheno:
    input:
        sg_check= os.path.join(
            '{info_d}', 'seq2geno_report'),
        genml= os.path.join('{info_d}', 'seq2geno.gml')
    output: 
        gp_check= os.path.join('{info_d}', 'geno2pheno_report')
    conda: os.path.join(geno2pheno_home, 'installation', 'requirements.yml')
    params:
        ovrd= args.ovrd
    threads: args.cores
    log: os.path.join('{info_d}', 'geno2pheno_log')
    threads: args.cores
    shell:
        '''
        geno2pheno.py --genoparse {input.genml} \
--override {params.ovrd} \
--cores {threads}
        printf "Geno2Pheno: $(date)\ndetails recorded in: {log}" > {output.gp_check};
        '''
