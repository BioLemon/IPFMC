import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering, KMeans
import collections
from sklearn.metrics import adjusted_rand_score, precision_score, f1_score, silhouette_score
from snf import make_affinity, snf
import warnings
from sklearn import metrics
warnings.filterwarnings("ignore")


def Count_NMIARI(Truelabel, Prediclabel):
    # Calculate the Adjusted Rand Index (ARI) between true labels and predicted labels
    Spectral_ARI = adjusted_rand_score(Truelabel, Prediclabel)
    # Calculate the Normalized Mutual Information (NMI) score between true labels and predicted labels
    Spectual_NMI = metrics.normalized_mutual_info_score(Truelabel, Prediclabel, average_method='arithmetic')
    return Spectral_ARI, Spectual_NMI


def Count_Precision(y_true, y_pred):
    # Get a unique list of predicted values to use as row names for the confusion matrix
    row_names = np.unique(y_pred)
    # Get a unique list of true label values to use as column names for the confusion matrix
    col_names = np.unique(y_true)
    # Initialize a zero matrix to store the values of the confusion matrix
    cm = np.zeros((len(row_names), len(col_names)), dtype=int)

    # Iterate over the true labels and predicted labels to update the confusion matrix
    for i in range(len(y_true)):
        # Update the confusion matrix based on the indices of the true and predicted labels
        cm[np.where(row_names == y_pred[i])[0], np.where(col_names == y_true[i])[0]] += 1

    # Convert the confusion matrix to a DataFrame with specified row and column names
    df = pd.DataFrame(cm, index=row_names, columns=col_names)
    original_df = df.copy()  # Keep a copy of the original confusion matrix
    max_list = []  # List to store maximum values from the confusion matrix
    position_list = []  # List to store positions of maximum values

    # While the DataFrame is not empty, find and process the maximum values
    while not df.empty:
        max_value = df.max().max()  # Find the maximum value in the DataFrame
        max_positions = df.stack().index[df.stack() == max_value].tolist()  # Get positions of the maximum value

        if len(max_positions) == 1:  # If there is a single maximum position
            max_list.append(max_value)  # Append the maximum value to the list
            position_list.append(max_positions[0])  # Append the position to the list
            df = df.drop(max_positions[0][0], axis=0)  # Drop the row corresponding to the maximum position
            df = df.drop(max_positions[0][1], axis=1)  # Drop the column corresponding to the maximum position
        else:  # If there are multiple maximum positions
            max_dict = {}  # Dictionary to store temporary maximum values
            for position in max_positions:
                temp_df = df.copy()  # Create a temporary copy of the DataFrame
                temp_df.loc[position] = 0  # Set the row of the maximum position to zero
                temp_df[position[1]] = 0  # Set the column of the maximum position to zero
                # Calculate the new maximum value after dropping the current position
                max_dict[position] = temp_df.loc[position[0]].max() + temp_df[position[1]].max()

            # Find the position with the minimum maximum value
            min_max = min(max_dict.values())
            min_key = [key for key, value in max_dict.items() if value == min_max][0]
            max_list.append(df.loc[min_key])  # Append the maximum value from the original DataFrame
            position_list.append(min_key)  # Append the position to the list
            df = df.drop(min_key[0], axis=0)  # Drop the row corresponding to the minimum maximum position
            df = df.drop(min_key[1], axis=1)  # Drop the column corresponding to the minimum maximum position

    # Calculate the precision based on the maximum values found
    PrecisionAns = sum(max_list)
    Fullsum = np.sum(original_df.to_numpy())  # Total sum of the original confusion matrix
    PrecisionAns /= Fullsum  # Calculate the precision

    F1_micro = 0  # Initialize F1-micro value
    F_list = []  # List to store F1 scores for each position
    Num_Samples = np.sum(np.array(original_df))  # Total number of samples

    # Calculate F1-micro value
    for pos in position_list:
        TP = original_df.loc[pos[0], pos[1]]  # True Positives
        FP = np.sum(original_df.loc[:, pos[1]]) - TP  # False Positives
        FN = np.sum(original_df.loc[pos[0], :]) - TP  # False Negatives
        Prcs = TP / (TP + FP)  # Precision
        Recl = TP / (TP + FN)  # Recall
        if Prcs == Recl == 0:  # If both precision and recall are zero
            F_part = 0  # Set F1 part to zero
        else:
            F_part = (2 * Prcs * Recl) / (Prcs + Recl)  # Calculate F1 part
        F_part = F_part * (TP + FP) / Num_Samples  # Weight the F1 part by the number of samples
        F_list.append(F_part)  # Append the F1 part to the list

    F_final = np.sum(F_list)  # Sum the F1 parts to get the final F1 score
    return PrecisionAns, position_list, F_final  # Return precision, positions, and final F1 score


def Fuse_Matrices(DataTypes_list, Cancer_type, Sk):
    # Initialize a list to hold matrices for each data type
    Matrix_list = []
    for DataType in DataTypes_list:
        # Construct the input file path for each data type
        Input = '../Datas/Omics/Output_Omics/' + Cancer_type + '/'
        # Load the matrix from the CSV file
        Matrix = np.genfromtxt(Input + Cancer_type + '_' + DataType + '_Matrix.csv', delimiter=',')
        Matrix_list.append(Matrix)  # Append the matrix to the list

    # Fuse the matrices using the SNF algorithm
    fusion_matrix = snf(Matrix_list, K=Sk)
    return fusion_matrix  # Return the fused matrix


def Find_Best_K(Simular_matrix):  # Input similarity matrix and return the best k value
    S_list = []  # List to store silhouette scores for different k values
    for k in range(2, 9):  # Test k values from 2 to 8
        S_templist = []  # Temporary list to hold silhouette scores for each k
        for i in range(10):  # Repeat the clustering 10 times for each k
            Sc_label = SpectralClustering(n_clusters=k, affinity='precomputed').fit_predict(Simular_matrix)
            SScore = silhouette_score(Simular_matrix, Sc_label)  # Calculate silhouette score
            S_templist.append(SScore)  # Append the score to the temporary list
        Scounter = collections.Counter(S_templist)  # Count occurrences of each score
        most_common_S = Scounter.most_common(1)  # Find the most common silhouette score
        S_list.append(most_common_S)  # Append the most common score to the list

    k = S_list.index(max(S_list)) + 2  # Find the index of the maximum score and adjust for range
    return k  # Return the best k value


def Get_Labels(Simular_matrix, k):
    # Perform spectral clustering on the similarity matrix to get cluster labels
    Sc_label = SpectralClustering(n_clusters=k, affinity='precomputed').fit_predict(Simular_matrix)
    return Sc_label  # Return the cluster labels


# Import similarity matrices, omics data, and survival data
# List of cancer types to analyze
Cancer_type_list = ['BRCA_GOLD', 'COAD_GOLD']
# List of combinations of data types to be fused
Data_Combine_list = [['mRNA', 'miRNA'], ['mRNA', 'Methy'], ['mRNA', 'CNV'], ['miRNA', 'Methy'], ['miRNA', 'CNV'],
                     ['Methy', 'CNV'], ['mRNA', 'miRNA', 'Methy'], ['mRNA', 'miRNA', 'CNV'], ['mRNA', 'Methy', 'CNV'],
                     ['miRNA', 'Methy', 'CNV'], ['mRNA', 'miRNA', 'Methy', 'CNV']]
# Data type for the fused data
Data_type = 'Fused'
TRUE_Table = []  # Initialize a list to store results
K = 15  # Set the number of clusters for fusion

# Loop through each cancer type in the cancer type list
for Cancer_type in Cancer_type_list:
    # Define input directory for the current cancer type
    input_dir = '../Datas/Omics/Output_Omics/' + Cancer_type + '/'
    # Define the path to the gold standard labels
    Gold_dir = '../Datas/Omics/Raw_Omics/' + Cancer_type + '/' + Cancer_type[:4] + '_label.csv'
    # Load the gold standard labels from the CSV file
    Gold_label = np.genfromtxt(Gold_dir, delimiter=',', dtype='str')
    print(Gold_label)  # Print the gold labels for verification
    Gold_label = Gold_label[1:]  # Exclude the header from the gold labels

    # Initialize lists to store metrics for the current cancer type
    ARI_list = [Cancer_type, 'ARI']
    NMI_list = [Cancer_type, 'NMI']
    k_list = [Cancer_type, 'k']
    Preci_list = [Cancer_type, 'Precision']
    F_list = [Cancer_type, 'F-measure']

    # Loop through each combination of data types to be fused
    for Combines in Data_Combine_list:
        # Generate the similarity matrix based on the current combination of data types
        Mean_Matrix = Fuse_Matrices(Combines, Cancer_type, K)
        # Determine the best k value based on the similarity matrix
        K1 = Find_Best_K(Mean_Matrix)
        # Get predicted labels using the best k value
        Pred_labels = Get_Labels(Mean_Matrix, K1)
        # Calculate ARI and NMI metrics
        ARI, NMI = Count_NMIARI(Pred_labels, Gold_label)
        # Append the calculated metrics to their respective lists
        ARI_list.append(ARI)
        NMI_list.append(NMI)
        k_list.append(K1)

        # Determine the number of clusters for prediction based on cancer type
        if Cancer_type == 'BRCA_GOLD':
            Prepred_labels = Get_Labels(Mean_Matrix, 5)  # Use 5 clusters for BRCA_GOLD
        else:
            Prepred_labels = Get_Labels(Mean_Matrix, 4)  # Use 4 clusters for other types

        # Calculate precision and F1 score based on the predicted labels
        preci, posi, F1 = Count_Precision(Gold_label, Prepred_labels)
        # Append precision and F1 score to their respective lists
        F_list.append(F1)
        Preci_list.append(preci)

    # Append the lists of metrics for the current cancer type to the TRUE_Table
    TRUE_Table.append(ARI_list)
    TRUE_Table.append(NMI_list)
    TRUE_Table.append(k_list)
    TRUE_Table.append(Preci_list)
    TRUE_Table.append(F_list)

# Convert the TRUE_Table to a NumPy array
TRUE_Table = np.array(TRUE_Table)
# Create a DataFrame to organize the results with appropriate column names
P_table = pd.DataFrame(TRUE_Table, columns=['Cancer', 'Datatype', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
print(P_table)  # Print the results DataFrame for verification

# Define output directory for saving the results
Output_dir = "../Datas/Assessment_Result/Trues"
# Save the results DataFrame to an Excel file
P_table.to_excel(Output_dir + f'True_Result_complete.xlsx', sheet_name='TRUE')









