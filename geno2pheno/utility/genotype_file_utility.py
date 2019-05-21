__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu ehsaneddin.asgari@helmholtz-hzi.de"

import sys
sys.path.append('../')

import codecs
import numpy as np
from scipy import sparse
from sklearn.preprocessing import MaxAbsScaler
from utility.file_utility import FileUtility
from sklearn import preprocessing
import os
from sklearn.feature_extraction.text import TfidfVectorizer
import itertools

class GenotypeReader(object):
    '''
        Class for reading genotype files
    '''

    def __init__(self):
        print('Genotype reader object created..')

    @staticmethod
    def get_nuc_kmer_distribution(sequences,k):
        vocab = [''.join(xs) for xs in itertools.product('atcg', repeat=k)]
        vocab.sort()
        vectorizer = TfidfVectorizer(use_idf=False, vocabulary=vocab, analyzer='char', ngram_range=(k, k),
                                          norm='l1', stop_words=[], lowercase=True, binary=False)
        return (np.mean(vectorizer.fit_transform(sequences).toarray(), axis=0),vocab)


    @staticmethod
    def create_read_tabular_file(path, save_pref='_', feature_normalization=None, transpose=False, override=False):
        '''
        :param path:
        :param save_pref:
        :param transpose: if isolates are columns
        :param feature_normalization: 'binary': {0,1}, '0-1': [0-1],  'percent': {0,1,..,100}, 'zu': zero mean, unit variance
        :return:
        '''
        print ('Start creating ', save_pref)
        if override  or not os.path.exists('_'.join([save_pref, 'feature', 'vect.npz'])):
            rows = [l.strip() for l in codecs.open(path, 'r', 'utf-8').readlines()]
            tf_vec = sparse.csr_matrix(
                [[GenotypeReader.get_float_or_zero(x) for x in entry.split('\t')[1::]] for entry in rows[1::]])

            if transpose:
                tf_vec = sparse.csr_matrix(tf_vec.toarray().T)
                isolates = [feat.replace(' ', '') for feat in rows[0].rstrip().split('\t')]
                feature_names = [row.split()[0] for row in rows[1::]]
            else:
                isolates = [row.split()[0] for row in rows[1::]]
                feature_names = [feat.replace(' ', '') for feat in rows[0].rstrip().split('\t')]

            # normalizer / discretizer
            if feature_normalization:
                if feature_normalization == 'binary':
                    tf_vec = np.round(MaxAbsScaler().fit_transform(tf_vec))

                elif feature_normalization == '01':
                    tf_vec = MaxAbsScaler().fit_transform(tf_vec)
                elif feature_normalization == 'percent':
                    tf_vec = np.round(MaxAbsScaler().fit_transform(tf_vec) * 100)
                elif feature_normalization == 'zu':
                    tf_vec = sparse.csr_matrix(preprocessing.StandardScaler().fit_transform(tf_vec.toarray()))


                FileUtility.save_sparse_csr('_'.join([save_pref, 'feature', 'vect.npz']), tf_vec)
                FileUtility.save_list('_'.join([save_pref, 'feature', 'list.txt']), feature_names)
                FileUtility.save_list('_'.join([save_pref, 'strains', 'list.txt']), isolates)
                print (save_pref, ' created successfully containing ', str(len(isolates)), ' strains and ', str(len(feature_names)), ' features')
                return (''.join([save_pref, ' created successfully containing ', str(len(isolates)), ' strains and ', str(len(feature_names)), ' features']))
        else:
            print (save_pref, ' already exist ')
            return (''.join([save_pref, ' already exist ']))

    @staticmethod
    def get_float_or_zero(value):
        '''
        to replace any non numeric value with 0
        :param value:
        :return:
        '''
        try:
            return float(value)
        except:
            return 0.0
