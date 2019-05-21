#!/usr/bin/env python3

'''
As an interface to link seq2geno, genml creator, and geno2pheno
'''

if __name__=='__main__':
    import sys
    import os

    sgp_home= os.path.dirname(os.path.realpath(__file__))
    sys.path.append(os.path.join(sgp_home, 'lib'))

    seq2geno_home= os.path.join(os.environ['SEQ2GENO_HOME'], 'main')
    sys.path.append(seq2geno_home)
    ## accept arguments
    import UserOptionsSGP
    args= UserOptionsSGP.main()

    ## load Seq2Geno module
    try:
        import seq2geno
        seq2geno.main(args)
    except Exception as e:
        sys.exit('ERROR: {}'.format(str(e)))

    ## load genml creator module
    import create_genml_from_seq2geno
    create_genml_from_seq2geno.main(args)
    
    ## call Geno2Pheno
    import os
    import subprocess
    gp_cmd= ['Geno2Pheno.py', 
            '--genoparse', os.path.join(args.sg, 'seq2geno.gml'),
            '--override', str(args.ovrd),
            '--cores', str(args.cores)]
    subprocess.run(gp_cmd)
