import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering
import lifelines
import collections
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.cluster.hierarchy import linkage, cut_tree
from snf import make_affinity, snf
import warnings
from scipy.linalg import eigh
from scipy.sparse.csgraph import laplacian
from scipy.sparse import csr_matrix

warnings.filterwarnings("ignore")
def Count_P(Mean_Matrix,Surv_data,ID):
    kmf = lifelines.KaplanMeierFitter()
    pses = []
    Sses = []
    for k in range(2, 9):
        p_value_list = []
        silhouette_list = []
        sc = SpectralClustering(n_clusters=k, affinity="precomputed")
        for i in range(10):
            labels = sc.fit_predict(Mean_Matrix)
            labels = pd.DataFrame(labels, index=ID)
            SurAna_data = pd.concat([Surv_data, labels], axis=1)
            SurAna_data.rename(columns={'patient_dmfs_e': 'event', 'patient_dmfs_time': 'time', 0: 'label'},
                               inplace=True)
            results = lifelines.statistics.multivariate_logrank_test(SurAna_data['time'], SurAna_data['label'],
                                                                     SurAna_data['event'])
            p_value_list.append(results.p_value)
            silhouette = metrics.silhouette_score(Mean_Matrix, labels)
            silhouette_list.append(silhouette)
        Pcounter = collections.Counter(p_value_list)
        most_common_P = Pcounter.most_common(1)
        p_value = most_common_P[0][0]
        Scounter = collections.Counter(silhouette_list)
        most_common_S = Scounter.most_common(1)
        silhouette = most_common_S[0][0]
        Sses.append(silhouette)
        pses.append(p_value)
        #print(f"Cluster number is {k}, logrank-p value is: {p_value}, silhouette coefficient is {silhouette}")
    return pses,Sses

def Fuse_Matrices(DataTypes_list,Cancer_type):
    Matrix_list = []
    for DataType in DataTypes_list:
        Input = '/home/zhanghaoy/Pytorchdeep/Datas/Omics/Output_Omics/' + Cancer_type + '/'
        Matrix = np.genfromtxt(Input + Cancer_type + '_' + DataType + '_Matrix.csv', delimiter=',')
        #print(len(Matrix))
        Matrix_list.append(Matrix)
        #print(Matrix)
    # fusion_matrix = np.nanmean(Matrix_list, axis=0)
    fusion_matrix = snf(Matrix_list, K=15)

    #print(f'Completed matrix fusion normally, matrix length is {len(fusion_matrix)}')
    return fusion_matrix



# Import similarity matrix, omics data, survival data
Cancer_type_list = ['BRCA','COAD','KIRC','LUAD','LUSC','ACC','KIRP','LIHC','THYM']
#Cancer_type_list = ['LUAD']  # Cancer type list
Data_Combine_list = [['mRNA','miRNA'],['mRNA','Methy'],['mRNA','CNV'],['miRNA','Methy'],['miRNA','CNV'],
                     ['Methy','CNV'],['mRNA','miRNA','Methy'],['mRNA','miRNA','CNV'],['mRNA','Methy','CNV'],
                     ['miRNA','Methy','CNV'],['mRNA','miRNA','Methy','CNV']]
#Data_Combine_list = [['mRNA','miRNA']]
P_table = []
for Cancer_type in Cancer_type_list:
    # input_dir = '../Datas/Omics/Output_Omics/' + Cancer_type + '/'
    # Gold_dir = '../Datas/Omics/Raw_Omics/' + Cancer_type + '/'
    Surv_dir = '../Datas/Survival/Subset_Survival/' + Cancer_type

    Surv_data = pd.read_csv(Surv_dir + '_Sub_survival.csv', index_col=0)  # Import original survival data
    IDs = Surv_data.index
    Combine_Name = ''  # Used for later visualization to show which combination is currently used
    for Combines in Data_Combine_list:
        Mean_Matrix = Fuse_Matrices(Combines,Cancer_type)
        # Mean_Matrix = np.genfromtxt(input_dir + Cancer_type + '_' + Data_type + '_Matrix.csv', delimiter=',')
        P_values, S_values = Count_P(Mean_Matrix,Surv_data,IDs)
        print(f'P value list of {Combines} combination for {Cancer_type} cancer is: \n {P_values}')
        #print(f'Silhouette coefficient value list of {Combines} combination for {Cancer_type} cancer is: \n {S_values}')
        A = 0
        for x in P_values:
            if x < 0.05:
                A = A + 1
        # Use max function to find the maximum value in the list
        max_value = max(S_values)
        # Use index method to find the index of the maximum value
        Best_K = S_values.index(max_value)
        B = Best_K+2
        P_values.append(A)
        P_values.append(B)
        P_values.insert(0, Cancer_type)
        P_table.append(P_values)
P_table = np.array(P_table)
print(P_table)
P_table = pd.DataFrame(P_table,columns=['Cancer',2,3,4,5,6,7,8,'No_Sig','Sug_K'])
print(P_table)
Output_dir = "../Datas/Assessment_Result/"
P_table.to_excel(Output_dir+'Complete_Result.xlsx',sheet_name='Complete')
