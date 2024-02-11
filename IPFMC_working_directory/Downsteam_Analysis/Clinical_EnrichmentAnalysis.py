import pandas as pd
import numpy as np
import scipy.stats as stats
from snf import make_affinity,snf
from sklearn.cluster import SpectralClustering
from scipy.stats import f_oneway

def test_enrichment(clinical_matrix,clinical_param):
    Clinical_list = list(clinical_matrix.columns)
    Clinical_list = Clinical_list[1:]
    Clinical_list = clinical_param & set(Clinical_list)
    p_values = []
    for clinical_param in Clinical_list:
        chi2, p = stats.chi2_contingency(pd.crosstab(clinical_matrix['Cluster'], clinical_matrix[clinical_param]))[:2]
        p_values.append(p)
    return p_values

def Fuse_Matrices(DataTypes_list, Sk):
    Matrix_list = []
    for DataType in DataTypes_list:
        Input = '../Datas/Omics/Output_Omics/' + Cancer_type + '/'
        Matrix = np.genfromtxt(Input + Cancer_type + '_' + DataType + '_Matrix.csv', delimiter=',')
        # print(len(Matrix))
        Matrix_list.append(Matrix)
        # print(Matrix)
    fusion_matrix = snf(Matrix_list, K=Sk)
    # print(f'length of the matrix: {len(fusion_matrix)}')
    return fusion_matrix


Cancer_type_list = ['BRCA','COAD','KIRC','LUAD','LUSC','ACC','KIRP','LIHC','THYM']
Clinical_param = {'pathologic_M', 'pathologic_N', 'pathologic_T', 'pathologic_stage','masaoka_stage','clinical_M'}
Data_Combine_list = [['mRNA','miRNA'],['mRNA','Methy'],['mRNA','CNV'],['miRNA','Methy'],['miRNA','CNV'],
                     ['Methy','CNV'],['mRNA','miRNA','Methy'],['mRNA','miRNA','CNV'],['mRNA','Methy','CNV'],
                     ['miRNA','Methy','CNV'],['mRNA','miRNA','Methy','CNV']]

Clinical_dir = '../Datas/Clinical/PCB'
Omic_dir = '../Datas/Omics/Raw_Omics'
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
            sc = SpectralClustering(n_clusters=k, affinity='precomputed')
            Mean_Matrix = Fuse_Matrices(Combine,23)
            labels = sc.fit_predict(Mean_Matrix)
            labels = pd.DataFrame(labels)
            labels.columns = ['Cluster']
            labels.index = Omic_IDs
            labels = labels.loc[Same_IDs,:]
            Full_Clinical = pd.concat([labels, Clinical_data], axis=1)
            ps = test_enrichment(Full_Clinical,Clinical_param)
            Num_Sig = 0
            for x in ps:
                if x < 0.05:
                    Num_Sig = Num_Sig + 1
            Total_list.append(Num_Sig)
            Num_list.append(Num_Sig)
        Total_table.append(Total_list)
Total_table = np.array(Total_table)
Total_table = pd.DataFrame(Total_table,columns=['Cancer',2,3,4,5,6,7,8])
Output_dir = "../Datas/Assessment_Result/"
Total_table.to_excel(Output_dir+'Clinical_Result.xlsx',sheet_name='Complete')