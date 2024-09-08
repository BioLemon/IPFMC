import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering
import lifelines
import collections
from sklearn import metrics
from snf import snf
import warnings

warnings.filterwarnings("ignore")


def Count_P(Mean_Matrix, Surv_data, ID):
    """
    Calculates p-values and silhouette coefficients for different cluster numbers.

    Args:
        Mean_Matrix (numpy.ndarray): The input matrix for clustering.
        Surv_data (pandas.DataFrame): Survival data.
        ID (list): List of patient IDs.

    Returns:
        list: List of p-values for each cluster number.
        list: List of silhouette coefficients for each cluster number.
    """
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
    return pses, Sses


def Fuse_Matrices(DataTypes_list, Cancer_type):
    """
    Fuses multiple data matrices using Similarity Network Fusion (SNF).

    Args:
        DataTypes_list (list): List of data types to be fused.
        Cancer_type (str): Cancer type.

    Returns:
        numpy.ndarray: Fused matrix.
    """
    Matrix_list = []
    for DataType in DataTypes_list:
        Input = '../Datas/Omics/Output_Omics/' + Cancer_type + '/'
        Matrix = np.genfromtxt(Input + Cancer_type + '_' + DataType + '_Matrix.csv', delimiter=',')
        Matrix_list.append(Matrix)
    fusion_matrix = snf(Matrix_list, K=15)
    return fusion_matrix


# List of cancer types to be processed
Cancer_type_list = ['BRCA', 'COAD', 'KIRC', 'LUAD', 'LUSC', 'ACC', 'KIRP', 'LIHC', 'THYM']

# List of data combinations to be processed
Data_Combine_list = [['mRNA', 'miRNA'], ['mRNA', 'Methy'], ['mRNA', 'CNV'], ['miRNA', 'Methy'], ['miRNA', 'CNV'],
                     ['Methy', 'CNV'], ['mRNA', 'miRNA', 'Methy'], ['mRNA', 'miRNA', 'CNV'], ['mRNA', 'Methy', 'CNV'],
                     ['miRNA', 'Methy', 'CNV'], ['mRNA', 'miRNA', 'Methy', 'CNV']]

# Initialize an empty list to store the results
P_table = []

# Iterate over each cancer type
for Cancer_type in Cancer_type_list:
    # Load the survival data for the current cancer type
    Surv_dir = '../Datas/Survival/Subset_Survival/' + Cancer_type
    Surv_data = pd.read_csv(Surv_dir + '_Sub_survival.csv', index_col=0)
    IDs = Surv_data.index

    # Iterate over each data combination
    for Combines in Data_Combine_list:
        # Fuse the matrices for the current data combination and cancer type
        Mean_Matrix = Fuse_Matrices(Combines, Cancer_type)
        # Calculate the p-values and S-values for the fused matrix and survival data
        P_values, S_values = Count_P(Mean_Matrix, Surv_data, IDs)

        # Count the number of significant p-values (< 0.05)
        A = 0
        for x in P_values:
            if x < 0.05:
                A = A + 1
        # Find the best K value based on the maximum S-value
        max_value = max(S_values)
        Best_K = S_values.index(max_value)
        B = Best_K + 2  # Suggested K is Best_K + 2

        # Append the results to the P_table list
        P_values.append(A)
        P_values.append(B)
        P_values.insert(0, Cancer_type)
        P_table.append(P_values)

# Convert the P_table list to a NumPy array and then to a DataFrame
P_table = np.array(P_table)
P_table = pd.DataFrame(P_table, columns=['Cancer', 2, 3, 4, 5, 6, 7, 8, 'No_Sig', 'Sug_K'])

# Save the results to an Excel file
Output_dir = "../Datas/Assessment_Result/"
P_table.to_excel(Output_dir + 'Complete_Result.xlsx', sheet_name='Complete')