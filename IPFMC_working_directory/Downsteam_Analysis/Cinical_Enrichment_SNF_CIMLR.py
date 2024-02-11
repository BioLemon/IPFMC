import pandas as pd
import numpy as np
import scipy.stats as stats
from snf import make_affinity,snf
from sklearn.cluster import SpectralClustering
from scipy.stats import f_oneway
import os
from collections import Counter
def test_enrichment(clinical_matrix,clinical_param):
    Clinical_list = list(clinical_matrix.columns)
    Clinical_list = Clinical_list[1:]
    Clinical_list = clinical_param & set(Clinical_list)
    p_values = []
    for clinical_param in Clinical_list:
        chi2, p = stats.chi2_contingency(pd.crosstab(clinical_matrix['Cluster'], clinical_matrix[clinical_param]))[:2]
        p_values.append(p)
    return p_values

def Get_Right_Label(cancer,combine,k,method):
    Temp_str = ''
    for data in combine:
        Temp_str = f'{Temp_str}_{data}'
    A_Labels_dir = f'{Cluster_dir}/{Cancer_type}/{Method}/{k}{Temp_str}.csv'  # 第一种格式类型
    Temp_str = ''
    for data in combine:
        Temp_str = f'{Temp_str}{data}_'
    B_Labels_dir = f'{Cluster_dir}/{Cancer_type}/{Method}/{Temp_str}{k}.csv'  # 第一种格式类型

    if os.path.exists(A_Labels_dir):
        label_file = pd.read_csv(A_Labels_dir)
    elif os.path.exists(B_Labels_dir):
        label_file = pd.read_csv(B_Labels_dir)
    else:
        labels = [1]
        return labels
    labels = label_file.iloc[:,1]
    print(len(labels))
    return labels

Cancer_type_list = ['ACC','BRCA','COAD','KIRC','KIRP','LIHC','LUAD','LUSC','THYM']
Clinical_param = {'pathologic_M', 'pathologic_N', 'pathologic_T', 'pathologic_stage','masaoka_stage','clinical_M'}
Data_Combine_list = [['m','mi'],['m','me'],['m','cnv'],['mi','me'],['mi','cnv'],
                     ['me','cnv'],['m','mi','me'],['m','mi','cnv'],['m','me','cnv'],
                     ['mi','me','cnv'],['m','mi','me','cnv']]
Method = 'CIMLR'
Clinical_dir = '../Datas/Clinical/PCB'
Omic_dir = '../Datas/Omics/Raw_Omics'
Cluster_dir = '../Datas/Cluster_Results'
Total_table = []
Num_list = []
for Cancer_type in Cancer_type_list:
    Clinical_data = pd.read_csv(f'{Clinical_dir}/{Cancer_type}_Sub_Clinical.csv')
    Omic_data = pd.read_csv(f'{Omic_dir}/{Cancer_type}/{Cancer_type}_miRNA.csv')
    Omic_data = Omic_data.transpose()
    Omic_IDs = list(Omic_data.index[1:])
    Clinical_IDs = list(Clinical_data['sampleID'])
    Clinical_data.index = Clinical_data['sampleID']
    Clinical_data = Clinical_data.iloc[:, 1:]
    for i in range(len(Omic_IDs)):
        Omic_IDs[i] = Omic_IDs[i].replace('.', '-')
    Same_IDs = set(Omic_IDs) & set(Clinical_IDs)
    Same_IDs = list(Same_IDs)
    Clinical_data = Clinical_data.loc[Same_IDs, :]
    for Combine in Data_Combine_list:
        Total_list = [Cancer_type]
        for k in range(2,9):
            labels = Get_Right_Label(Cancer_type,Combine,k,Method)
            if len(labels) == 1:
                Total_list.append('NA')
                continue
            #print(labels)
            labels = labels.rename('Cluster')
            #print(labels)
            labels.index = Omic_IDs
            labels = labels.loc[Same_IDs]
            Full_Clinical = pd.concat([labels, Clinical_data], axis=1)
            #print(Full_Clinical)
            ps = test_enrichment(Full_Clinical,Clinical_param)
            print(Cancer_type)
            print(ps)
            Num_Sig = 0
            for x in ps:
                if x < 0.05:
                    Num_Sig = Num_Sig + 1
            Total_list.append(Num_Sig)
            Num_list.append(Num_Sig)
            print(Num_Sig)
        Total_table.append(Total_list)
Total_table = np.array(Total_table)
Num_list = np.array(Num_list)
Ave_Num = np.mean(Num_list)
print(Total_table)
Total_table = pd.DataFrame(Total_table,columns=['Cancer',2,3,4,5,6,7,8])
Output_dir = "../Datas/Assessment_Result/"
Total_table.to_excel(Output_dir+f'{Method}_Clinical_Result.xlsx',sheet_name='Complete')


