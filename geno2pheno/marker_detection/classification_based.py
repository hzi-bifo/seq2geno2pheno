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

def generate_top_features(path, classifier_list, topk=200):
    ## TODO: ask as an input topk

    writer = pd.ExcelWriter(path+'/ultimate_outputs/selected_features.xls', engine='xlsxwriter')

    final_results=dict()
    for classifier in classifier_list:
        feature_files = FileUtility.recursive_glob(path+'/feature_selection/','*_'+classifier)
        res=dict()
        for file in feature_files:
            phenotype=file.split('/')[0:-1][-1]
            if not phenotype in res:
                res[phenotype]=[file]
            else:
                if file.split('/')[-1].count('##')>res[phenotype][0].split('/')[-1].count('##'):
                    res[phenotype]=[file]
                elif file.split('/')[-1].count('##')==res[phenotype][0].split('/')[-1].count('##'):
                    res[phenotype].append(file)
        for phenotype in res.keys():
            if phenotype not in final_results:
                final_results[phenotype]=[]
            final_results[phenotype]+=res[phenotype]
    for phenotype,files in final_results.items():
        selected=[{x.split('\t')[0]:1/(idx+1) for idx, x in enumerate(FileUtility.load_list(file)[1:topk])} for file in files]
        res=set(selected[0])
        for set_select in selected[1::]:
             res=res.intersection(set_select)

        geno_val_res=dict()
        for dict_geno_val in selected:
            for x,val in dict_geno_val.items():
                if x not in geno_val_res:
                    geno_val_res[x]=[val,1]
                else:
                    geno_val_res[x][0]+=val
                    geno_val_res[x][1]+=1

        df_dict={'feature_name':[],'mrr':[],'freq_confirmation':[]}
        for name, values in geno_val_res.items():
            rr, nr = values
            df_dict['feature_name'].append(name)
            df_dict['mrr'].append(rr/nr)
            df_dict['freq_confirmation'].append(nr)
        df=pd.DataFrame(df_dict)
        df.sort_values(['freq_confirmation','mrr','feature_name'], ascending=[False, False, False], inplace=True)
        df=df.copy()
        df.to_excel(writer, sheet_name=phenotype, index=False)

