#!/usr/bin/env python3
def parse_usr_opts():
    import argparse
    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='create genml file for Geno2Pheno')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')
    parser.add_argument('--seq2geno', dest= 'sg', 
            help= 'seq2geno project folder', required= True)
    parser.add_argument('--out', dest= 'out', 
            help= 'output folder', default= 'geno2pheno.results')

    ## prediction block
    pred_args=  parser.add_argument_group('predict')
    pred_args.add_argument('--pred', dest= 'pred',
            help= 'name of prediction (e.g. classX_vs_classY)', 
            default= 'classX_vs_classY')
    pred_args.add_argument('--opt', dest= 'optimize', 
            help= 'target performance metric to optimize', 
            choices=['scores_f1_1'])
    pred_args.add_argument('--fold_n', dest= 'fold_n', 
            help= 'number of folds during validation', 
            default= 10)
    pred_args.add_argument('--test_perc', dest= 'test_perc', 
            help= 'proportion of samples for testing', 
            default= 0.1)
    pred_args.add_argument('--part', dest= 'part', 
            help= 'method to partition dataset', 
            choices= ['rand'])
    pred_args.add_argument('--models', dest= 'models', 
            help= 'machine learning algorithms', 
            choices= ['svm', 'rf', 'lr'], nargs= '*')
    pred_args.add_argument('--k-mer', dest= 'kmer', 
            help= 'the k-mer size for prediction', 
            default= 6)
    pred_args.add_argument('--cls', dest= 'classes_f', 
            help= 'a two-column file to specify label and prediction group')
    
    args = parser.parse_args()
    return(args)

def main(args):
    from xml.dom import minidom
    import xml.etree.ElementTree as ET
    import os
    import sys
    ## test folder
    if not os.path.isdir(os.path.join(args.sg, 'RESULTS')):
        print('ERROR: Cannot find {}'.format(os.path.join(args.sg, 'RESULTS')))
        sys.exit()

    root = ET.Element("project",
            output=args.gp_output_d,
            name=args.sg)

    ####
    ## genotype block
    geno = ET.SubElement(root, "genotype")
    #### for binary tables
    bin_tab_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS',
        'bin_tables'))
    if os.path.isdir(bin_tab_dir):
        bin_tables = ET.SubElement(
            geno, "tables",
            attrib= {
                'path': bin_tab_dir, 
                'normalization': "binary", 
                'transpose': "False"})
        setattr(bin_tables, 'text', args.sg)
    #### for numeric features
    num_tab_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS',
        'num_tables'))
    if os.path.isdir(num_tab_dir):
        con_tables = ET.SubElement(
            geno, "tables",
            attrib= {
                'path': num_tab_dir, 
                'normalization': "numeric", 
                'transpose': "False"})
        setattr(con_tables, 'text', args.sg)
    #### genome seq
    assem_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS', 
        'assemblies'))
    if os.path.isdir(assem_dir):
        seq= ET.SubElement(
            geno, "sequence",
            attrib={
                'path': assem_dir, 
                'kmer': str(args.kmer)})
        setattr(seq, 'text', args.sg)

    ####
    ## phenotype block
    phe_f= os.path.abspath(os.path.join(args.sg, 'RESULTS', 'phenotype',
        'phenotypes.mat'))
    if os.path.isfile(phe_f):
        pheno = ET.SubElement(
            root, "phenotype",
            attrib={'path': phe_f})
        setattr(pheno, 'text', '\n')
    else:
        sys.exit('No phenotype was found')

    ####
    ## phylogeny block
    phy_f=os.path.abspath(os.path.join(args.sg, 'RESULTS', 'phylogeny',
        'tree.nwk'))
    if os.path.isfile(phy_f):
        phy= ET.SubElement(
            root, "phylogentictree",
            attrib={'path': phy_f})
        setattr(phy, 'text', '\n')
    else:
        sys.exit('No phylogeny was found')


    ####
    ## predict block
    pred= ET.SubElement(root, "predict",
            attrib= {'name': args.pred})
    optimize= ET.SubElement(pred, "optimize")
    setattr(optimize, 'text', str(args.optimize))
    validation= ET.SubElement(pred, "eval",
            attrib= {'folds': str(args.fold_n), 'test': str(args.test_perc)})
    setattr(validation, 'text', str(args.part))

    with open(args.classes_f, 'r') as cls_fh:
        for l in cls_fh:
            d=l.strip().split('\t')
            cls= ET.SubElement(pred, "label",
                    attrib= {'value': str(d[0])})
            setattr(cls, 'text', str(d[1]))

    model= ET.SubElement(pred, "model")
    for m in args.models:
        ET.SubElement(model, m)


    ####
    ## convert to string
    tree = ET.ElementTree(root)
    rough_string= ET.tostring(root, 'utf-8')
    reparsed= minidom.parseString(rough_string)
    indented_string= reparsed.toprettyxml(indent="  ")
    genml_f= os.path.join(args.sg, 'seq2geno.gml')
    with open(genml_f, 'w') as genml_fh:
            genml_fh.write(indented_string)


if __name__== '__main__':
    args= parse_usr_opts()
    main(args)
