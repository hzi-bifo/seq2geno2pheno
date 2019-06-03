#!/usr/bin/env python3

'''
Folder -> check if absent or empty
File -> check if absent 
'''

def parse_usr_opts():
    import argparse
    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='create genml file from Seq2Geno output')

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

def main(args, genml_f):
    from xml.dom import minidom
    import xml.etree.ElementTree as ET
    import os
    import sys
    ## test folder
    if not os.path.isdir(os.path.join(args.sg, 'RESULTS')):
        raise OSError('{} is absent'.format(
            os.path.join(args.sg, 'RESULTS')))

    root = ET.Element("project",
            output=os.path.abspath(args.gp_output_d),
            name=args.sg)

    ####
    ## genotype block
    geno = ET.SubElement(root, "genotype")
    #### for binary tables
    bin_tab_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS',
        'bin_tables'))
    if os.path.isdir(bin_tab_dir) & (len(os.listdir(bin_tab_dir) ) > 0):
        bin_tables = ET.SubElement(
            geno, "tables",
            attrib= {
                'path': bin_tab_dir, 
                'normalization': "binary", 
                'transpose': "False"})
        setattr(bin_tables, 'text', 'binary_tables')
    else:
        raise OSError(
                'GenML: {} (binary features) is absent but required'.format(bin_tables_dir))

    #### for numeric features
    num_tab_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS',
        'num_tables'))
    if os.path.isdir(num_tab_dir) & (len(os.listdir(num_tab_dir) ) > 0):
        con_tables = ET.SubElement(
            geno, "tables",
            attrib= {
                'path': num_tab_dir, 
                'normalization': "zu", 
                'transpose': "False"})
        setattr(con_tables, 'text', 'numeric_tables')
    else:
        raise OSError(
                'GenML: {} (numeric features) is absent but required'.format(num_tables_dir))

    #### genome seq
    assem_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS', 
        'assemblies'))
    if os.path.isdir(assem_dir) & (len(os.listdir(assem_dir) ) > 0):
        seq= ET.SubElement(
            geno, "sequence",
            attrib={
                'path': assem_dir, 
                'kmer': str(args.kmer)})
        setattr(seq, 'text', 'assemblies')
    else:
        raise OSError(
                'GenML: {} (assemblies) is absent but required'.format(assem_dir))

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
        raise OSError(
                'GenML: {} (phenotypes) is absent but required'.format(phe_f))

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
        raise OSError('GenML: {} (phylogeny) is absent but required'.format(phy_f))


    ####
    ## predict block
    pred= ET.SubElement(root, "predict",
            attrib= {'name': args.pred})
    optimize= ET.SubElement(pred, "optimize")
    setattr(optimize, 'text', str(args.optimize[0]))
    validation= ET.SubElement(pred, "eval",
            attrib= {'folds': str(args.fold_n), 'test': str(args.test_perc)})
    setattr(validation, 'text', str(args.part))

    if os.path.isfile(args.classes_f):
        all_cls= ET.SubElement(pred, "labels")
        with open(args.classes_f, 'r') as cls_fh:
            for l in cls_fh:
                d=l.strip().split('\t')
                cls= ET.SubElement(all_cls, "label",
                        attrib= {'value': str(d[1])})
                setattr(cls, 'text', str(d[0]))
    else:
        raise OSError(
                'GenML: {} (prediction classes) is absent but required'.format(args.classes_f))

    model= ET.SubElement(pred, "model")
    if len(args.models) > 0:
        for m in args.models:
            ET.SubElement(model, m)
    else:
        raise ValueError('Prediction models not specified')

    ####
    ## convert to string
    tree = ET.ElementTree(root)
    rough_string= ET.tostring(root, 'utf-8')
    reparsed= minidom.parseString(rough_string)
    indented_string= reparsed.toprettyxml(indent="  ")
    with open(genml_f, 'w') as genml_fh:
            genml_fh.write(indented_string)


if __name__== '__main__':
    args= parse_usr_opts()
    main(args)
