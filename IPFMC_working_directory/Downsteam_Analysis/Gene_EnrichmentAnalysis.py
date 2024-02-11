import warnings
import pandas as pd
import numpy as np
from snf import snf
from sklearn.cluster import SpectralClustering
import scipy.stats as stats
import lifelines as ll
import matplotlib.pyplot as plt
from itertools import combinations
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
warnings.filterwarnings(action='ignore')

Cancer_type_list = ['ACC','BRCA','COAD','KIRC','KIRP','LIHC','LUAD','LUSC','THYM']
Data_type_list = ['mRNA','miRNA','Methy','CNV']
Cancer_type = 'LUAD'
Data_type = 'miRNA'

num = 100
Omic_dir = '../Datas/Omics/Raw_Omics'
Full_Path_dir = f'../Datas/Pathways/miRNA_Pathway_Index.csv'
Sel_Path_dir = f'../Datas/Omics/Output_Omics'
Surv_dir = '../Datas/Survival/Subset_Survival'

Omic_data = pd.read_csv(f'{Omic_dir}/{Cancer_type}/{Cancer_type}_{Data_type}.csv', index_col=0)
Surv_data = pd.read_csv(f'{Surv_dir}/{Cancer_type}_Sub_survival.csv')
Full_Pathways = pd.read_csv(Full_Path_dir,index_col=0)
Sel_Pathways = pd.read_excel(f'{Sel_Path_dir}/{Cancer_type}_{Data_type}_Pathway_Sorting.xlsx',index_col=0)
Full_Genes = list(Omic_data.index)  
Sel_Pathways = list(Sel_Pathways['Pathway'])  
Parcial_Pathways = Full_Pathways.loc[Sel_Pathways]
print(Parcial_Pathways)

count_dict = {}

for item in Full_Genes:
    count_dict[item] = 0

for Pathway in Sel_Pathways:
    Genes = Parcial_Pathways.loc[Pathway]
    Genes = list(Genes.dropna())
    for i in Genes:
        if i in Full_Genes:
            count_dict[i] += 1
Gene_Count = pd.DataFrame.from_dict(count_dict, orient='index', columns=['count'])
Gene_Count = Gene_Count.sort_values(by='count', ascending=False)
Gene_Count = Gene_Count.iloc[:num]

print(Gene_Count)
Gene_Count.to_excel(f"../Datas/Assessment_Result/Gene_Enrich/{Data_type}_retopgenes.xlsx")
