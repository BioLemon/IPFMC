import pandas as pd
import numpy as np
from snf import snf
from IPFMC import direct, separate, analysis

# Filepath of the omics data, ‘LUAD’ is the folder contains omics datas of LUAD cancer
Omic_dir = './Omics/Demo'
# Filepath of the pathway index
BP_dir = './Pathways/Pathway_Index.csv'
# Filepath of the miRNA pathway index
mirBP_dir = 'Pathways/miRNA_Pathway_Index.csv'
datatypes = ['mRNA', 'Methy', 'CNV', 'miRNA']  # The type of data to be used in the experiment
omic_list = []  # A list for storing multiple omics data
BP_data = pd.read_csv(BP_dir, index_col=0)  # The pandas package is used to pass in the pathway data
mirBP_data = pd.read_csv(mirBP_dir, index_col=0)  # Pass in the pathway-mirna relationship data
for datatype in datatypes:
    '''
    We named the omics data <cancer name>_<data type>.csv, for example, LUAD_mRNA.csv
    You can change it according to your habits
    '''
    omicdata = pd.read_csv(f'{Omic_dir}/Demo_{datatype}.csv', index_col=0)
    omic_list.append(omicdata)

represents = []
pathways_list = []
# Only the first three data sets are processed here, and the last data set is miRNA, which needs to be processed separately
for i in range(3):
    represent, pathways = separate.ipfmc_discretize(omic_list[i], BP_data, seed=10)
    represents.append(np.array(represent))
    # print(represent)
    pathways_list.append(pathways)

represent, pathways = separate.ipfmc_discretize(omic_list[3], mirBP_data, seed=10)  # Here processes miRNA dataset
represents.append(np.array(represent))
pathways_list.append(pathways)
represent_final = snf(represents, K=15)  # 'represent_final' is the final multi-omics representation

K = separate.suggest_k(
    represent_final)  # input the final representation, and this function will give a suggested cluster
labels = separate.spec_cluster(omic_list[0], fusion_matrix=represent_final, k=K)
print(labels)

Gene_names = list(omic_list[0].index)  # obtain all gene names occured in the omics dataset
Gene_rank = analysis.gene_occurrence(Full_Genes=Gene_names, Sel_Pathways=pathways_list[0], Full_Pathways=BP_data)
print(Gene_rank)
miRNA_names = list(omic_list[3].index)  # obtain all miRNA names occured in the omics dataset
miRNA_rank = analysis.gene_occurrence(Full_Genes=miRNA_names, Sel_Pathways=pathways_list[3], Full_Pathways=mirBP_data)
print(miRNA_rank)
