__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu ehsaneddin.asgari@helmholtz-hzi.de"

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib


def create_correlation_mat(save_path,df,filetype='pdf', ifmask=False):

    sns.set(style="white")

    df=df.loc[:, (df != 0).any(axis=0)]
    # Compute the correlation matrix
    corr = df.corr()

    font_s=min(max(int(500/len(df.columns.tolist())),1),30)
    params = {
        'legend.fontsize': font_s,
        'xtick.labelsize': font_s,
        'ytick.labelsize': font_s,
        'text.usetex': True,
    }
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
    plt.rc('text', usetex=True)
    plt.rcParams.update(params)

    if np.min(corr.values)>=0:
        selected_cmap='Purples'
    else:
        selected_cmap='coolwarm'


    # Generate a mask for the upper triangle
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))

    if ifmask:
        # Draw the heatmap with the mask and correct aspect ratio
        sns.heatmap(corr, mask=mask, cmap=selected_cmap, vmax=.3, center=0,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5}, xticklabels=[x.replace('_','-') for x in df.columns.tolist()[0:-1]], yticklabels=[x.replace('_','-') for x in ['']+df.columns.tolist()[::-1][0:-1][::-1]])
    else:
        # Draw the heatmap with the mask and correct aspect ratio
        sns.heatmap(corr, cmap=selected_cmap, vmax=.3, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5},xticklabels=[x.replace('_','-') for x in df.columns.tolist()], yticklabels=[x.replace('_','-') for x in df.columns.tolist()])

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)

    plt.tight_layout()

    plt.savefig(save_path+'.'+filetype)

def create_kl_plot(mat, axis_names, title, filename, xlab, ylab, cmap='inferno', filetype='pdf', rx=90, ry=0, font_s=10,
                    annot=True):
    '''
    :param mat: divergence matrix
    :param axis_names: axis_names
    :param title
    :param filename: where to be saved
    :return:
    '''
    font_s=min(max(int(500/len(axis_names)),1),30)

    params = {
        'legend.fontsize': font_s,
        'xtick.labelsize': font_s,
        'ytick.labelsize': font_s,
        'text.usetex': True,
    }
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
    plt.rc('text', usetex=True)
    plt.rcParams.update(params)

    if len(axis_names) == 0:
        ax = sns.heatmap(mat, annot=annot, cmap=cmap, fmt="d")
    else:
        # removed fmt="d",
        ax = sns.heatmap(mat, annot=annot, yticklabels=[x.replace('_','-') for x in axis_names], xticklabels=[x.replace('_','-') for x in axis_names], cmap=cmap)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xticks(rotation=rx)
    plt.yticks(rotation=ry)
    plt.rcParams.update(params)
    plt.savefig(filename + '.' + filetype)
    plt.show()
    plt.clf()

