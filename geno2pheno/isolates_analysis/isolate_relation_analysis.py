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
import numpy as np
from drug_analysis.drug_relation_analysis import DrugRelation
from utility.visualization_utility import create_mat_plot
from sklearn.preprocessing import MinMaxScaler

class IsolatesRelations(object):
    '''
        Isolates Relations
    '''
    def __init__(self):
        # phenotype relation
        self.DrugRel = DrugRelation()

        self.phenotype_isolates = self.DrugRel.get_isolate_list()
        self.isolate_phenotype_dist = self.DrugRel.get_isolate_profile_kldiv()[0] + self.DrugRel.get_isolate_profile_kldiv()[0].T

        # phylogenetic relation
        self.phyl_isolates=[]
        self.isolate_phyl_distance=[]

        self.make_isolate_phlyogen_dist()

        # common_isolates
        self.common_isolates=list(set(self.phenotype_isolates).intersection(self.phyl_isolates))
        self.common_isolates.sort()
        self.common_dist_phylogen=np.zeros((len(self.common_isolates),len(self.common_isolates)))
        self.common_dist_phenotype=np.zeros((len(self.common_isolates),len(self.common_isolates)))
        self.make_common_distances()
        self.stndrd_common_dist_phylogen=MinMaxScaler().fit_transform(self.common_dist_phylogen.reshape([self.common_dist_phylogen.size, 1])).reshape(self.common_dist_phylogen.shape)
        self.stndrd_common_dist_phenotype=MinMaxScaler().fit_transform(self.common_dist_phenotype.reshape([self.common_dist_phenotype.size, 1])).reshape(self.common_dist_phenotype.shape)

    def create_phenotype_phylog_image(self):
        '''
        :return:
        '''
        create_mat_plot(self.stndrd_common_dist_phylogen, ['']*len(self.common_isolates), 'Phylogenetic distance between isolates', 'results/isolate_analysis/common_phyl', cmap='Purples', filetype='png')
        create_mat_plot(self.stndrd_common_dist_phenotype,['']*len(self.common_isolates), 'Phenotypical distance between isolates', 'results/isolate_analysis/common_phen', cmap='Purples', filetype='png')

    def make_common_distances(self):
        '''
        Produce the common distance matrices
        '''
        for i in range(len(self.common_isolates)):
            for j in range(i+1, len(self.common_isolates)):
                self.common_dist_phylogen  [i,j] = self.get_isolate_phlyogen_dist(self.common_isolates[i], self.common_isolates[j])
                self.common_dist_phylogen  [j,i] = self.common_dist_phylogen [i,j]
                self.common_dist_phenotype [i,j] = self.get_isolate_phenotype_dist(self.common_isolates[i], self.common_isolates[j])
                self.common_dist_phenotype [j,i] = self.common_dist_phenotype [i,j]

    def make_isolate_phlyogen_dist(self):
        '''
        Make the similarity matrix
        :return:
        '''
        self.phyl_isolates=[line.split()[0] for line in codecs.open('../data_config/distance.txt','r','utf-8').readlines()]
        self.isolate_phyl_distance=np.array([[float(x) for x in line.split()[1::]] for line in codecs.open('../data_config/distance.txt','r','utf-8').readlines()])

    def get_isolate_phlyogen_dist(self, isol1, isol2):
        '''
        Get the distance score
        :return:
        '''
        if isol1 not in self.phyl_isolates or isol2 not in self.phyl_isolates:
            return 'NULL'
        else:
            return self.isolate_phyl_distance[self.phyl_isolates.index(isol1), self.phyl_isolates.index(isol2)]


    def get_isolate_phenotype_dist(self, isol1, isol2):
        '''
        Get the similarity score
        :param isol1:
        :param isol2:
        :return:
        '''
        if isol1 not in self.phenotype_isolates or isol2 not in self.phenotype_isolates:
            return 'NULL'
        else:
            return self.isolate_phenotype_dist[self.phenotype_isolates.index(isol1), self.phenotype_isolates.index(isol2)]


if __name__ == "__main__":
    IR = IsolatesRelations()

