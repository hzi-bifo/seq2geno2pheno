__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "GENO2PHENO of SEQ2GENO2PHENO"
__website__ = ""

import sys

sys.path.append('../')
from utility.genotype_file_utility import GenotypeReader
from utility.file_utility import FileUtility
import os
import tqdm
from scipy import sparse
from multiprocessing import Pool


class IntermediateRepCreate(object):
    '''
    This class is written to save
    the representations
    '''

    def __init__(self, output_path):
        '''
        :param output_path:
        '''
        print('data creator..')
        self.output_path = output_path
        if not os.path.exists(output_path):
            os.makedirs(output_path)

    def create_table(self, path, name, feature_normalization, override=False):

        return GenotypeReader.create_read_tabular_file(path, save_pref=self.output_path + name,
                                                       feature_normalization=feature_normalization, override=override)

    def _get_kmer_rep(self, inp):
        strain, seq_file, k=inp
        seq=FileUtility.read_fasta_sequences(seq_file)
        vec, vocab = GenotypeReader.get_nuc_kmer_distribution(seq, k)
        return strain, vec, vocab



    def create_kmer_table(self, path, k, cores=4, override=False):

        save_path = self.output_path + 'sequence_' + str(k) + 'mer'

        if override  or not os.path.exists('_'.join([save_path, 'feature', 'vect.npz'])):
            files = FileUtility.recursive_glob(path, '*')
            files.sort()
            input_tuples = []
            for file in files:
                input_tuples.append((file.split('/')[-1].split('.')[0], file, k))

            strains = []
            mat = []
            kmers=[]
            pool = Pool(processes=cores)
            for strain, vec, vocab in tqdm.tqdm(pool.imap_unordered(self._get_kmer_rep, input_tuples, chunksize=cores),
                                                total=len(input_tuples)):
                strains.append(strain)
                mat.append(vec)
                kmers=vocab
            pool.close()
            mat = sparse.csr_matrix(mat)

            FileUtility.save_sparse_csr(save_path+'_feature_vect', mat)
            FileUtility.save_list('_'.join([save_path, 'strains', 'list.txt']), strains)
            FileUtility.save_list('_'.join([save_path, 'feature', 'list.txt']), kmers)
        return ('_'.join([save_path]) + ' created')



if __name__ == "__main__":
    IC = IntermediateRepCreate(
        '/net/sgi/metagenomics/projects/pseudo_genomics/results/amr_toolkit/testingpack/intermediate_rep/')
    IC.create_table(
        '/net/sgi/metagenomics/projects/pseudo_genomics/results/PackageTesting/K_pneumoniae/genotables/gpa.uniq.mat',
        'uniqGPA', 'binary')
    IC.create_table(
        '/net/sgi/metagenomics/projects/pseudo_genomics/results/PackageTesting/K_pneumoniae/genotables/non-syn_SNPs.uniq.mat',
        'uniqNonsynSNP', 'binary')
    IC.create_table(
        '/net/sgi/metagenomics/projects/pseudo_genomics/results/PackageTesting/K_pneumoniae/genotables/syn_SNPs.uniq.mat',
        'uniqNonsynSNP', 'binary')
