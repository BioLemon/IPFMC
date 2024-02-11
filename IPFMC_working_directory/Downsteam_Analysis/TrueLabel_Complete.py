import itertools
import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering, KMeans
import collections
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import adjusted_rand_score, precision_score, f1_score, silhouette_score
from snf import make_affinity, snf
import warnings
from sklearn import metrics
from scipy.linalg import eigh
from scipy.sparse.csgraph import laplacian
from scipy.sparse import csr_matrix

warnings.filterwarnings("ignore")


def Count_NMIARI(Truelabel, Prediclabel):
    Spectral_ARI = adjusted_rand_score(Truelabel, Prediclabel)
    Spectual_NMI = metrics.normalized_mutual_info_score(Truelabel, Prediclabel, average_method='arithmetic')
    return Spectral_ARI, Spectual_NMI


def Count_Precision(y_true, y_pred):
    row_names = np.unique(y_pred)
    col_names = np.unique(y_true)
    cm = np.zeros((len(row_names), len(col_names)), dtype=int)
    for i in range(len(y_true)):
        cm[np.where(row_names == y_pred[i])[0], np.where(col_names == y_true[i])[0]] += 1
    df = pd.DataFrame(cm, index=row_names, columns=col_names)
    original_df = df.copy()
    max_list = []
    position_list = []
    while not df.empty:
        max_value = df.max().max()
        max_positions = df.stack().index[df.stack() == max_value].tolist()

        if len(max_positions) == 1:
            max_list.append(max_value)
            position_list.append(max_positions[0])
            df = df.drop(max_positions[0][0], axis=0)
            df = df.drop(max_positions[0][1], axis=1)
        else:
            max_dict = {}
            for position in max_positions:
                temp_df = df.copy()
                temp_df.loc[position] = 0
                temp_df[position[1]] = 0
                max_dict[position] = temp_df.loc[position[0]].max() + temp_df[position[1]].max()
            min_max = min(max_dict.values())
            min_key = [key for key, value in max_dict.items() if value == min_max][0]
            max_list.append(df.loc[min_key])
            position_list.append(min_key)
            df = df.drop(min_key[0], axis=0)
            df = df.drop(min_key[1], axis=1)
    PrecisionAns = sum(max_list)
    Fullsum = np.sum(original_df.to_numpy())
    PrecisionAns /= Fullsum
    F1_micro = 0
    F_list = []
    Num_Samples = np.sum(np.array(original_df))
    for pos in position_list:
        TP = original_df.loc[pos[0], pos[1]]
        FP = np.sum(original_df.loc[:, pos[1]]) - TP
        FN = np.sum(original_df.loc[pos[0], :]) - TP
        Prcs = TP / (TP + FP)
        Recl = TP / (TP + FN)
        if Prcs == Recl == 0:
            F_part = 0
        else:
            F_part = (2 * Prcs * Recl) / (Prcs + Recl)
        F_part = F_part * (TP + FP) / Num_Samples
        F_list.append(F_part)
    F_final = np.sum(F_list)
    # print(F_final)
    return PrecisionAns, position_list, F_final



def get_combinations(data):
    combinations = []
    for i in itertools.combinations(data, 2):
        temp = list(i)
        temp.sort()
        combinations.append(tuple(temp))
    return combinations


def Fuse_Matrices(DataTypes_list, Cancer_type, Sk):
    Matrix_list = []
    for DataType in DataTypes_list:
        Input = '/home/zhanghaoy/Pytorchdeep/Datas/Omics/Output_Omics/' + Cancer_type + '/'
        Matrix = np.genfromtxt(Input + Cancer_type + '_' + DataType + '_Matrix.csv', delimiter=',')
        # print(len(Matrix))
        Matrix_list.append(Matrix)
        # print(Matrix)
    # fusion_matrix = np.nanmean(Matrix_list, axis=0)
    fusion_matrix = snf(Matrix_list, K=Sk)
    return fusion_matrix


def Find_Best_K(Simular_matrix):
    S_list = []
    for k in range(2, 9):
        S_templist = []
        for i in range(10):
            Sc_label = SpectralClustering(n_clusters=k, affinity='precomputed').fit_predict(Mean_Matrix)
            SScore = silhouette_score(Simular_matrix, Sc_label)
            S_templist.append(SScore)
        Scounter = collections.Counter(S_templist)
        most_common_S = Scounter.most_common(1)
        S_list.append(most_common_S)
    k = S_list.index(max(S_list)) + 2
    return k


def Get_Labels(Simular_matrix, k):
    Sc_label = SpectralClustering(n_clusters=k, affinity='precomputed').fit_predict(Simular_matrix)
    return Sc_label

Cancer_type_list = ['BRCA_GOLD', 'COAD_GOLD']
Data_Combine_list = [['mRNA', 'miRNA'], ['mRNA', 'Methy'], ['mRNA', 'CNV'], ['miRNA', 'Methy'], ['miRNA', 'CNV'],
                     ['Methy', 'CNV'], ['mRNA', 'miRNA', 'Methy'], ['mRNA', 'miRNA', 'CNV'], ['mRNA', 'Methy', 'CNV'],
                     ['miRNA', 'Methy', 'CNV'], ['mRNA', 'miRNA', 'Methy', 'CNV']]
Data_type = 'Fused'
TRUE_Table = []
for Cancer_type in Cancer_type_list:
    input_dir = '/home/zhanghaoy/Pytorchdeep/Datas/Omics/Output_Omics/' + Cancer_type + '/'
    Gold_dir = '/home/zhanghaoy/Pytorchdeep/Datas/Omics/Raw_Omics/' + Cancer_type + '/' + Cancer_type[:4] + '_label.csv'
    Gold_label = np.genfromtxt(Gold_dir, delimiter=',', dtype='str')
    print(Gold_label)
    Gold_label = Gold_label[1:]
    ARI_list = [Cancer_type, 'ARI']
    NMI_list = [Cancer_type, 'NMI']
    k_list = [Cancer_type, 'k']
    Preci_list = [Cancer_type, 'Precision']
    F_list = [Cancer_type, 'F-measure']
    for Combines in Data_Combine_list:
        Mean_Matrix = Fuse_Matrices(Combines, Cancer_type, 16)
        K1 = Find_Best_K(Mean_Matrix)
        Pred_labels = Get_Labels(Mean_Matrix, K1)
        ARI, NMI = Count_NMIARI(Pred_labels, Gold_label)
        ARI_list.append(ARI)
        NMI_list.append(NMI)
        k_list.append(K1)
        if Cancer_type == 'BRCA_GOLD':
            Prepred_labels = Get_Labels(Mean_Matrix, 5)
        else:
            Prepred_labels = Get_Labels(Mean_Matrix, 4)
        preci, posi, F1 = Count_Precision(Gold_label, Prepred_labels)
        F_list.append(F1)
        Preci_list.append(preci)
    TRUE_Table.append(ARI_list)
    TRUE_Table.append(NMI_list)
    TRUE_Table.append(k_list)
    TRUE_Table.append(Preci_list)
    TRUE_Table.append(F_list)
TRUE_Table = np.array(TRUE_Table)
P_table = pd.DataFrame(TRUE_Table, columns=['Cancer', 'Datatype', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
print(P_table)
Output_dir = "/home/zhanghaoy/Pytorchdeep/Datas/Assessment_Result/Trues"
P_table.to_excel(Output_dir + f'True_Result.xlsx', sheet_name='TRUE')

