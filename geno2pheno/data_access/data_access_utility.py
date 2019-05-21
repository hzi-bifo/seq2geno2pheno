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
from scipy import sparse
from sklearn.feature_extraction.text import TfidfTransformer
import numpy as np
from data_access.data_create_utility import ABRDataCreate


class ABRDataAccess(object):
    '''
    This class is written to load data for ABR prediction/analysis of Pseudomonas Aeruginosa
    '''

    def __init__(self, dir, prefix_list):
        '''
        :param dir: directory of features
        :param prefix_list: list of prefixes for features
        '''
        print ('Data access created..')
        # dictionary of feature types
        self.X = dict()
        # dictionary of feature names
        self.feature_names = dict()
        # dictionary of isolates
        self.isolates = dict()
        ## load the features
        self.load_data(dir, prefix_list)
        ## load the labels
        self.BasicDataObj = ABRDataCreate()

    def load_data(self, dir, prefix_list):
        '''
        Load list of features
        :param dir:
        :param prefix_list:
        :return:
        '''
        for save_pref in prefix_list:
            print ('@@@'+'_'.join([dir + save_pref, 'feature', 'vect.npz']))
            self.X[save_pref] = FileUtility.load_sparse_csr('_'.join([dir + save_pref, 'feature', 'vect.npz']))
            self.feature_names[save_pref] = FileUtility.load_list('_'.join([dir + save_pref, 'feature', 'list.txt']))
            self.isolates[save_pref] = FileUtility.load_list('_'.join([dir + save_pref, 'isolates', 'list.txt']))

    def get_xy_prediction_mats(self, drug, mapping={'0': 0, '0.0': 0, '1': 1, '1.0': 1}, features_for_idf=[]):
        '''
        :param drug:
        :param mapping:
        :param features_for_idf: if needed..
        :return:
        '''
        ## find a mapping from isolates to the labels
        mapping_isolate_label = dict(self.BasicDataObj.get_new_labeling(mapping)[drug])

        # get common isolates
        list_of_list_of_isolates = list(self.isolates.values())
        list_of_list_of_isolates.append(list(mapping_isolate_label.keys()))
        final_isolates = ABRDataAccess.common_isolates(list_of_list_of_isolates)
        final_isolates.sort()

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
            idx = [self.isolates[feature_type].index(isol) for isol in final_isolates]
            temp = temp[idx, :]
            feature_matrices.append(temp.toarray())
            feature_names += ['##'.join([feature_type, x]) for x in self.feature_names[feature_type]]

        X = np.concatenate(tuple(feature_matrices), axis=1)
        X = sparse.csr_matrix(X)
        Y = [mapping_isolate_label[isol] for isol in final_isolates]
        return X, Y, feature_names, final_isolates

    def get_xy_multidrug_prediction_mats(self, features_for_idf=[]):
        '''
        :param features_for_idf: if needed..
        :return:
        '''
        ## find a mapping from isolates to the labels
        mapping_isolate_label = dict(self.BasicDataObj.get_multilabel_label_dic())

        # get common isolates
        list_of_list_of_isolates = list(self.isolates.values())
        list_of_list_of_isolates.append(list(mapping_isolate_label.keys()))
        final_isolates = ABRDataAccess.common_isolates(list_of_list_of_isolates)
        final_isolates.sort()

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
            idx = [self.isolates[feature_type].index(isol) for isol in final_isolates]
            temp = temp[idx, :]
            feature_matrices.append(temp.toarray())
            feature_names += ['##'.join([feature_type, x]) for x in self.feature_names[feature_type]]

        X = np.concatenate(tuple(feature_matrices), axis=1)
        X = sparse.csr_matrix(X)
        Y = [mapping_isolate_label[isol] for isol in final_isolates]
        return X, Y, feature_names, final_isolates

    @staticmethod
    def common_isolates(list_of_list_of_isolates):
        '''
        :param list_of_list_of_isolates:
        :return:
        '''
        common_islt = set(list_of_list_of_isolates[0])
        for next_list in list_of_list_of_isolates[1::]:
            common_islt = common_islt.intersection(next_list)
        common_islt = list(common_islt)
        common_islt.sort()
        return common_islt
