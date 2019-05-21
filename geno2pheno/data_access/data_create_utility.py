__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "GENO2PHENO of SEQ2GENO2PHENO"
__website__ = ""

import sys
sys.path.append('../')
import codecs
from utility.genotype_file_utility import GenotypeReader
from scipy.sparse import csr_matrix
from utility.file_utility import FileUtility
import pandas as pd
from sklearn.preprocessing import MinMaxScaler



class ABRDataCreate(object):
    '''
    This class is written to save
    data for AMR prediction of Pseudomonas Aeruginosa
    '''

    def __init__(self):
        '''
        To creat the ABD DATA
         self.isolate2label_vec_mapping  ZG02420619 ['1.0', '0.0', '0.0', '', '1.0']
         self.drugs
         self.labeled_isolates: list of sorted isolates
         self.drug2labeledisolates_mapping : {'Ceftazidim': {'BK1496/1': 0, 'BK2061/1': 0, 'CF561_Iso4': 1, 'CF577_Iso7': 1, 'CF590_Iso5': 0.5 ..
        '''

        # init
        self.isolate2label_vec_mapping = dict()
        self.labeled_isolates = list()
        self.drugs = list()
        self.drug2labeled_isolates_mapping = dict()

        # basic loading
        self.make_labels()

    def make_labels(self):
        '''
            This function load labels ZG02420619 ['1.0', '0.0', '0.0', '', '1.0']
        '''
        label_file_address ='/net/sgi/metagenomics/projects/pseudo_genomics/data/K_pneumoniae/v2/list/phenotype_list' #'/mounts/data/proj/asgari/dissertation/git_repos/less_important/amrprediction/data_config/pheno_table_CLSI_S-vs-R.txt'#'/net/sgi/metagenomics/projects/pseudo_genomics/data/MIC/v3/pheno_table_CLSI_S-vs-R.txt'#'/mounts/data/proj/asgari/dissertation/git_repos/amrprediction/data_config/pheno_table_CLSI_S-vs-R.txt'
        #'/mounts/data/proj/asgari/dissertation/git_repos/amrprediction/data_config/mic_bin_without_intermediate.txt'
        #'/net/sgi/metagenomics/projects/pseudo_genomics/data/MIC/v3/pheno_table_CLSI_S-vs-R.txt'
        rows = [l.replace('\n', '') for l in codecs.open(label_file_address, 'r', 'utf-8').readlines()]
        self.isolate2label_vec_mapping = {
            str(entry.split('\t')[0]): [str(x) for idx, x in enumerate(entry.split('\t')[1::])] for entry in rows[1::]}
        self.labeled_isolates = list(self.isolate2label_vec_mapping)
        self.labeled_isolates.sort()
        self.drugs = [x.split('_')[0] for x in rows[0].rstrip().split('\t')[1::]]
        # init
        for drug in self.drugs:
            self.drug2labeled_isolates_mapping[drug] = []
        mapping = {'0': 0, '0.0': 0, '1': 1, '1.0': 1, '': 0.5}
        # only consider non-empty values
        for isolate, resist_vec in self.isolate2label_vec_mapping.items():
            for idx, val in enumerate(resist_vec):
                if val in mapping:
                    self.drug2labeled_isolates_mapping[self.drugs[idx]].append((isolate, mapping[val]))
        # generate dict of labels for each class
        for drug in self.drugs:
            self.drug2labeled_isolates_mapping[drug] = dict(self.drug2labeled_isolates_mapping[drug])

    def get_new_labeling(self, mapping={'0': 0, '0.0': 0, '1': 1, '1.0': 1}):
        '''
        Get new labeling
        Load labels
        :param mapping:
        :return:
        '''
        new_drug2labeled_isolates_mapping = dict()
        for drug in self.drugs:
            new_drug2labeled_isolates_mapping[drug] = []
        # only consider non-empty values
        for isolate, resist_vec in self.isolate2label_vec_mapping.items():
            for idx, val in enumerate(resist_vec):
                if val in mapping:
                    new_drug2labeled_isolates_mapping[self.drugs[idx]].append((isolate, mapping[val]))
        return new_drug2labeled_isolates_mapping

    @staticmethod
    def create_feature_files(path):
        '''
        :param path: output path
        :return:
        '''
        base_path='/net/sgi/metagenomics/projects/pseudo_genomics/'

        #### gene_exp
        GenotypeReader.create_read_tabular_file(base_path+'data/gene_expression/v2/rpg_414_log.txt', save_pref=path+'genexp', feature_normalization='zu')
        GenotypeReader.create_read_tabular_file(base_path+'data/gene_expression/v2/rpg_414_log.txt', save_pref=path+'genexp_percent', feature_normalization='percent')

        #### snp
        GenotypeReader.create_read_tabular_file(base_path+'results/featuresAnalysis/v2/non-syn_snps/non_syn_snps_aa_uq.uniq.txt', save_pref=path+'snps_nonsyn_trimmed', feature_normalization='binary')

        #### gpa
        GenotypeReader.create_read_tabular_file(base_path+'results/featuresAnalysis/v2/gpa/annot.uniq.txt', save_pref=path+'gpa_trimmed', feature_normalization='binary')

        #### gpa - roary
        GenotypeReader.create_read_tabular_file(base_path+'results/assembly/v2/roary/v5/out_95/indels/indel_annot.txt', save_pref=path+'gpa_roary', feature_normalization='binary')


        '''
        /mounts/data/proj/asgari/github_data/data/pseudomonas/data_v3/snps_nonsyn_trimmed  created successfully containing  426  isolates and  73475  features
        /mounts/data/proj/asgari/github_data/data/pseudomonas/data_v3/snps_nonsyn_envclin_trimmed  created successfully containing  442  isolates and  77748  features
        /mounts/data/proj/asgari/github_data/data/pseudomonas/data_v3/snps_all_envclin_trimmed  created successfully containing  442  isolates and  316168  features
        /mounts/data/proj/asgari/github_data/data/pseudomonas/data_v3/snps_all_full_trimmed  created successfully containing  426  isolates and  306527  features
        /mounts/data/proj/asgari/github_data/data/pseudomonas/data_v3/gpa  created successfully containing  508  isolates and  41872  features
        /mounts/data/proj/asgari/github_data/data/pseudomonas/data_v3/genexp_norm01  created successfully containing  426  isolates and  6026  features
        /mounts/data/proj/asgari/github_data/data/pseudomonas/data_v3/genexp_count  created successfully containing  426  isolates and  6026  features
        '''


    def get_multilabel_label_dic(self):
        '''
            drug resistance profile label
        '''
        mapping = {'': 'I', '0': 'S', '0.0': 'S', '1': 'R', '1.0': 'R'}
        return {k: ''.join([mapping[x] for x in list(v)]) for k, v in self.isolate2label_vec_mapping.items()}

    @staticmethod
    def get_common_isolates(list_of_list_of_isolates):
        '''
        :param list_of_list_of_isolates:
        :return:
        '''
        common_isolate = set(list_of_list_of_isolates[0])
        for next_list in list_of_list_of_isolates[1::]:
            common_isolate = common_isolate.intersection(next_list)
        common_isolate = list(common_isolate)
        common_isolate.sort()
        return common_isolate
    
    @staticmethod
    def create_continous_mics():
        '''
        
        '''
        scaler = MinMaxScaler()
        df=pd.read_table("../data_config/Final_MICs_16.06.16.txt")
        res=df[['Isolates','CIP MIC','TOB MIC','COL MIC','CAZ MIC','MEM MIC']]
        matrix=np.array([[float(str(x).replace('<=','').replace('≤','').replace('<=','').replace('≥','').replace('>=','')) for x in row] for row in res[['CIP MIC','TOB MIC','COL MIC','CAZ MIC','MEM MIC']].as_matrix()])
        # find nans [[(idx,idy) for idy,y in enumerate(x) if y] for idx, x in enumerate(np.isnan(matrix))]
        resistances=np.delete(matrix,[509],axis=0)
        isolates=[x[0] for idx, x in enumerate(list(df[['Isolates']].values)) if not idx==509]
        # scale to 0-1
        resistances=scaler.fit_transform(resistances)
        features=['CIP','TOB','COL','CAZ','MEM']
        base_path='/mounts/data/proj/asgari/dissertation/datasets/deepbio/pseudomonas/data_v3/continous_mic_vals'
        resistances=csr_matrix(resistances)
        FileUtility.save_sparse_csr(base_path+'_feature_vect',resistances)
        FileUtility.save_list(base_path+'_isolates_list.txt', isolates)
        FileUtility.save_list(base_path+'_feature_list.txt', features)





if __name__ == "__main__":
    ABRDataCreate.create_feature_files('/net/sgi/metagenomics/projects/pseudo_genomics/results/amr_toolkit/intermediate_reps/')

