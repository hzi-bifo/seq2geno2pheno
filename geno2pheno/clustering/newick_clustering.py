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
import pandas as pd
import itertools
from Bio import Phylo
import numpy as np
from sklearn.cluster import SpectralClustering
import os

def tree2mat_group(tree_file, n_group=20):
    '''
    This function maps the phylgenetic tree to adj matrix and perform spectral clustering on it
    :param tree_file:
    :param n_group:
    :return:
    '''
    if not os.path.exists(tree_file.replace(tree_file.split('/')[-1], 'phylogenetic_nodes_and_clusters.txt')):
        print ('Create phylogenetic information distance and groupings..')
        t = Phylo.read(tree_file, 'newick')
        d = {}
        for x, y in itertools.combinations(t.get_terminals(), 2):
            v = t.distance(x, y)
            d[x.name] = d.get(x.name, {})
            d[x.name][y.name] = v
            d[y.name] = d.get(y.name, {})
            d[y.name][x.name] = v
        for x in t.get_terminals():
            d[x.name][x.name] = 0

        m = pd.DataFrame(d)
        isolates = [x for x in m.axes[0]]
        isolates.sort()
        mat = np.zeros((len(isolates), len(isolates)))
        for x in range(len(isolates)):
            for y in range(len(isolates)):
                mat[x, y] = m[isolates[x]][isolates[y]]
        transferred_mat=np.exp(- mat ** 2 / (2. * 0.08 ** 2))
        clustering = SpectralClustering(n_clusters=n_group, assign_labels="kmeans", random_state=0).fit(transferred_mat)
        np.save(tree_file.replace(tree_file.split('/')[-1], 'phylogenetic_distance_matrix'), transferred_mat)
        FileUtility.save_list(tree_file.replace(tree_file.split('/')[-1], 'phylogenetic_nodes_and_clusters.txt'), ['\t'.join([x, str(clustering.labels_[idx])])for idx,x in enumerate(isolates)])
