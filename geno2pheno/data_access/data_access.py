__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "GENO2PHENO of SEQ2GENO2PHENO"
__website__ = ""

import sys

sys.path.append('../')
from utility.file_utility import FileUtility
from sklearn.feature_extraction.text import TfidfTransformer
import numpy as np
from scipy import sparse
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.model_selection import GroupKFold

class GenotypePhenotypeAccess(object):
    '''
        This class is written to handle genotype phenotype access
    '''

    def __init__(self, project_path, mapping=None):
        '''
        To creat the ABD DATA
         self.strain2labelvector  iso1 ['1.0', '0.0', '0.0', '', '1.0']
         self.drugs
         self.labeled_strains: list of sorted strains
         self.phenotype2labeled_strains_mapping : {'phenotype': {'iso1':0, ..
        '''

        self.project_path = project_path
        self.metadata_path = self.project_path + '/metadata/'
        self.representation_path = self.project_path + '/intermediate_rep/'
        # init
        self.strain2labelvector = dict()
        self.labeled_strains = list()
        self.phenotypes = list()
        self.phenotype2labeled_strains_mapping = dict()

        ## matrix creation variables
        # dictionary of feature types
        self.X = dict()
        # dictionary of feature names
        self.feature_names = dict()
        # dictionary of isolates
        self.strains = dict()

        # basic loading
        self.make_labels(mapping)

    def get_xy_prediction_mats(self, prefix_list, phenotype, mapping=None, features_for_idf=[]):
        '''
        :param phenotype:
        :param mapping:
        :param features_for_idf: if needed..
        :return:
        '''
        ## reset
        # dictionary of feature types
        self.X = dict()
        # dictionary of feature names
        self.feature_names = dict()
        # dictionary of isolates
        self.strains = dict()
        self.load_data(prefix_list)
        ## find a mapping from strains to the phenotypes
        if mapping:
            mapping_isolate_label = dict(self.get_new_labeling(mapping)[phenotype])
        else:
            mapping_isolate_label = self.phenotype2labeled_strains_mapping[phenotype]

        # get common strains
        list_of_list_of_strains = list(self.strains.values())
        list_of_list_of_strains.append(list(mapping_isolate_label.keys()))
        final_strains = GenotypePhenotypeAccess.get_common_strains(list_of_list_of_strains)
        final_strains.sort()

        # feature types
        feature_types = list(self.X.keys())
        # to apply idf if necessary
        if len(features_for_idf) > 0:
            tf = TfidfTransformer(norm=None, use_idf=True, smooth_idf=True)

        feature_names = []
        feature_matrices = []
        for feature_type in feature_types:
            if feature_type in features_for_idf:
                tf.fit(self.X[feature_type])
                temp = tf.transform(self.X[feature_type])
            else:
                temp = self.X[feature_type]
            idx = [self.strains[feature_type].index(strain) for strain in final_strains]
            temp = temp[idx, :]
            feature_matrices.append(temp.toarray())
            feature_names += ['##'.join([feature_type, x]) for x in self.feature_names[feature_type]]

        X = np.concatenate(tuple(feature_matrices), axis=1)
        X = sparse.csr_matrix(X)
        Y = [mapping_isolate_label[strain] for strain in final_strains]
        return X, Y, feature_names, final_strains

    def load_data(self, prefix_list=None):
        '''
        Load list of features
        :param dir:
        :param prefix_list:
        :return:
        '''
        for save_pref in prefix_list:
            print('@@@' + '_'.join([self.representation_path + save_pref, 'feature', 'vect.npz']))
            self.X[save_pref] = FileUtility.load_sparse_csr(
                '_'.join([self.representation_path + save_pref, 'feature', 'vect.npz']))
            self.feature_names[save_pref] = FileUtility.load_list(
                '_'.join([self.representation_path + save_pref, 'feature', 'list.txt']))
            self.strains[save_pref] = FileUtility.load_list(
                '_'.join([self.representation_path + save_pref, 'strains', 'list.txt']))

    def make_labels(self, mapping=None):
        '''
            This function load labels mapping from strain to phenotypes
        '''
        label_file_address = self.metadata_path + 'phenotypes.txt'
        rows = FileUtility.load_list(label_file_address)
        self.strain2labelvector = {
            str(entry.split()[0]): [str(x) for idx, x in enumerate(entry.split('\t')[1::])] for entry in rows[1::]}
        self.labeled_strains = list(self.strain2labelvector)
        self.labeled_strains.sort()

        self.phenotypes = [x for x in rows[0].rstrip().split()[1::]]
        # init
        for phenotype in self.phenotypes:
            self.phenotype2labeled_strains_mapping[phenotype] = []

        # only consider non-empty values
        for strain, phenotype_vec in self.strain2labelvector.items():
            for idx, val in enumerate(phenotype_vec):
                if mapping:
                    if val in mapping:
                        self.phenotype2labeled_strains_mapping[self.phenotypes[idx]].append((strain, mapping[val]))
                else:
                    self.phenotype2labeled_strains_mapping[self.phenotypes[idx]].append((strain, val))
        # generate dict of labels for each class
        for phenotype in self.phenotypes:
            self.phenotype2labeled_strains_mapping[phenotype] = dict(self.phenotype2labeled_strains_mapping[phenotype])

    def get_new_labeling(self, mapping):
        '''
        Get new labeling
        Load labels
        :param mapping:
        :return:
        '''
        new_phenotype2labeled_strain_mapping = dict()
        for drug in self.phenotypes:
            new_phenotype2labeled_strain_mapping[drug] = []
        # only consider non-empty values
        for strain, pheno_vec in self.strain2labelvector.items():
            for idx, val in enumerate(pheno_vec):
                if val in mapping:
                    new_phenotype2labeled_strain_mapping[self.phenotypes[idx]].append((strain, mapping[val]))
        return new_phenotype2labeled_strain_mapping

    def get_multilabel_label_dic(self, mapping=None):
        '''
            Phenotype multilabel
        '''
        if mapping:
            return {k: ''.join([mapping[x] for x in list(v)]) for k, v in self.strain2labelvector.items()}
        else:
            return {k: ''.join([x for x in list(v)]) for k, v in self.strain2labelvector.items()}

    def create_randfold(self, path, cv, test_ratio, phenotype, mapping=None):

        ## find a mapping from strains to the phenotypes
        if mapping:
            mapping_isolate_label = dict(self.get_new_labeling(mapping)[phenotype])
        else:
            mapping_isolate_label = self.phenotype2labeled_strains_mapping[phenotype]

        # get common strains
        list_of_list_of_strains = list(self.strains.values())
        list_of_list_of_strains.append(list(mapping_isolate_label.keys()))
        final_strains = GenotypePhenotypeAccess.get_common_strains(list_of_list_of_strains)
        final_strains.sort()

        # prepare test
        Y = [mapping_isolate_label[strain] for strain in final_strains]
        X_train, X_test, y_train, _ = train_test_split(final_strains, Y, test_size=test_ratio, random_state=0, stratify=Y)
        FileUtility.save_list(path.replace('_folds.txt', '_test.txt'), ['\t'.join(X_test)])

        # prepare train
        spliter=StratifiedKFold(cv)
        folds=['\t'.join([final_strains[x] for x in fold.tolist()]) for _,fold in  list(spliter.split(X_train,y_train))]
        FileUtility.save_list(path, folds)


    def create_treefold(self, path, tree_addr, cv, test_ratio, phenotype, mapping=None):

        ## find a mapping from strains to the phenotypes
        if mapping:
            mapping_isolate_label = dict(self.get_new_labeling(mapping)[phenotype])
        else:
            mapping_isolate_label = self.phenotype2labeled_strains_mapping[phenotype]

        # get common strains
        list_of_list_of_strains = list(self.strains.values())
        list_of_list_of_strains.append(list(mapping_isolate_label.keys()))
        final_strains = GenotypePhenotypeAccess.get_common_strains(list_of_list_of_strains)
        final_strains.sort()

        # prepare test
        Y = [mapping_isolate_label[strain] for strain in final_strains]

        isolate_to_group=dict([tuple(l.split('\t')) for l in FileUtility.load_list(tree_addr.replace(tree_addr.split('/')[-1], 'phylogenetic_nodes_and_clusters.txt'))])

        groups=[int(isolate_to_group[iso]) for iso in final_strains]
        group_kfold = GroupKFold(n_splits=round(1/test_ratio))
        for train_index, test_index in group_kfold.split(final_strains, Y, groups):
            X_test=[final_strains[x] for x in test_index]
            FileUtility.save_list(path.replace('_folds.txt', '_test.txt'), ['\t'.join(X_test)])
            break

        group_kfold = GroupKFold(n_splits=cv)

        folds=[]
        for train_index, test_index in group_kfold.split(train_index, [Y[idx] for idx in train_index],  [groups[idx] for idx in train_index]):
            folds.append(test_index)
        folds=['\t'.join([final_strains[x] for x in fold.tolist()]) for fold in  folds]
        FileUtility.save_list(path, folds)


    @staticmethod
    def get_common_strains(list_of_list_of_strains):
        '''
        :param list_of_list_of_strains:
        :return:
        '''
        common_strains = set(list_of_list_of_strains[0])
        for next_list in list_of_list_of_strains[1::]:
            common_strains = common_strains.intersection(next_list)
        common_strains = list(common_strains)
        common_strains.sort()
        return common_strains


if __name__ == "__main__":
    GPA = GenotypePhenotypeAccess(
        '/net/sgi/metagenomics/projects/pseudo_genomics/results/geno2pheno_package/K_pneumoniae/')
    print(GPA.strain2labelvector)
    print(GPA.phenotypes)
    print(GPA.phenotype2labeled_strains_mapping)
    X, Y, feature_names, final_strains = GPA.get_xy_prediction_mats('all', GPA.phenotypes[0])
    print(X.shape)
    print(len(Y))
    print(len(feature_names))
    print(len(final_strains))
