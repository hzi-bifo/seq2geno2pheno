class arguments:
    '''
    The object of arguments, as argparse is replaced
    '''
    def add_opt(self, **entries):
        self.__dict__.update(entries)
    def print_args(self):
        from pprint import pprint
        pprint(self.__dict__)
    def test_mandatory_args(self):
        import re
        ## ensure a key exists and the values is legal
        obligatory_args= ['dna_reads', 'pred', 'sgp_output_d']
        obligatory_args_values= {'dna_reads': '\w+', 
            'sgp_output_d': '\w+', 'pred': '\w+'}
        for k in obligatory_args:
            if hasattr(self, k):
                v= getattr(self, k)
                if re.search(obligatory_args_values[k], v) is None:
                    raise KeyError('"{}" has unrecognizable value'.format(k))
            else:
                raise KeyError(
                    '"{}" should be included in the options'.format(k))
        return(True)

    def check_args(self):
        import sys
        import os
        # obligatory arguments
        try:
            self.test_mandatory_args()
        except KeyError as e:
            sys.exit('ERROR: {}'.format(str(e)))

        # default values of optional arguments
        optional_args= {'cores':1, 'adaptor': '-', 'rna_reads': '-', 
                'dryrun': True, 
                'phe_table': '-', 'denovo': False, 'snps': False, 'expr': False,
                'phylo': False}
        for k in optional_args:
            if not hasattr(self, k):
                setattr(self, k, optional_args[k])

    def set_auto_filled_args(self):
        import os
        # auto-filled obligatory arguments
        au_obligatory_args= ['wd', 'sg', 'gp_output_d']
        au_obligatory_args_values= {
            'wd': os.path.join(
                getattr(self, 'sgp_output_d'), 'seq2geno'), 
            'sg': os.path.join(
                getattr(self, 'sgp_output_d'), 'seq2geno'),
            'gp_output_d': os.path.join(
                getattr(self, 'sgp_output_d'), 'geno2pheno')}
        for k in au_obligatory_args:
            setattr(self, k, au_obligatory_args_values[k])


def parse_arg_yaml(yml_f):
    '''
    Parse the yaml file where the parameters previously were commandline options
    '''
    import yaml
    import sys
    available_functions= ['snps', 'expr', 'denovo', 'phylo', 'de', 'ar',
    'dryrun']
    with open(yml_f, 'r') as yml_fh:
        opt_dict= yaml.safe_load(yml_fh)
        for k in available_functions: 
            if k in opt_dict['features']:
                opt_dict['features'][k]= (True if opt_dict['features'][k] == 'Y'
                else False)
            else:
                opt_dict['features'][k]=False 

    ## the seq2geno output folder should be fixed
    args= arguments()
    try:
        args.add_opt(**opt_dict['general'])
        args.add_opt(**opt_dict['features'])
        args.add_opt(**opt_dict['prediction'])
    except KeyError as e:
        sys.exit('ERROR: {} not found in the input file'.format(str(e)))
    else:
        args.check_args()
        args.set_auto_filled_args()
        return(args)

def read_options(opt_f, dsply_args):
    import sys
    args= parse_arg_yaml(opt_f)
    if dsply_args:
        print('Your arguments:')
        args.print_args()
        print('Only display options')
    return(args)

def main():
    '''
    Find the yaml file of arguments
    '''
    import yaml 
    import argparse
    import sys

    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='Seq2Geno2Pheno: the pipline tool '
                'for genomic features computation and phenotype predictor '
                'training')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')
    parser.add_argument('-d', dest= 'dsply_args', action= 'store_true',
        help= 'display the arguments in yaml and exit')
    parser.add_argument('-f', dest= 'yml_f', required= True, 
        help= 'the yaml file where the arguments are listed')

    primary_args= parser.parse_args()
    args= read_options(primary_args.yml_f, primary_args.dsply_args)
    return(args)

if __name__=='__main__':    
    main()
