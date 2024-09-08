import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import time
import warnings
warnings.filterwarnings("ignore")
def Eucli_Distance(df):
    # Calculate the Euclidean distance matrix from the input DataFrame
    euclidean_matrix = squareform(pdist(df, metric='euclidean'))  # Compute pairwise Euclidean distances and convert to a square matrix
    return euclidean_matrix  # Return the calculated Euclidean distance matrix

def CSPA(labels_list):
    """
    Cluster-based Similarity Partitioning Algorithm (CSPA)
    Input: a list of arrays of cluster labels, each array has shape (N,)
    Output: a consensus matrix, shape (N, N)
    """
    N = len(labels_list[0])  # Get the number of samples from the first labels array
    S_list = []  # Initialize a list to hold similarity matrices for each label set

    # Iterate over each set of cluster labels in the input list
    for labels in labels_list:
        S = np.zeros((N, N))  # Initialize a similarity matrix of shape (N, N) with zeros
        # Compare each pair of samples to determine similarity
        for i in range(N):
            for j in range(N):
                if labels[i] == labels[j]:  # Check if the labels of the two samples are the same
                    S[i][j] = 1  # If they are the same, set their similarity to 1
        S_list.append(S)  # Append the computed similarity matrix to the list

    # Compute the consensus matrix by averaging all similarity matrices
    C = np.mean(S_list, axis=0)  # Calculate the mean of the similarity matrices along the first axis
    return C, S_list  # Return the consensus matrix and the list of similarity matrices


start_time = time.time()  # Record the start time of the program
Cancer_type_list = ['BRCA', 'COAD', 'KIRC', 'LUAD', 'LUSC', 'ACC', 'KIRP', 'LIHC','THYM']  # List of cancer types to process

# Iterate over each cancer type in the list
for Cancer_type in Cancer_type_list:
    Data_type = 'miRNA'  # Define the type of omics data
    Omic_dir = '../Datas/Omics/Raw_Omics/' + Cancer_type + '/'  # Directory for omics data
    BP_dir = '../Datas/Pathways/miRNA_Pathway_Index.csv'  # Path to miRNA pathway index file

    # Load omics data based on whether the cancer type are for GOLD-standard evaluation purpose
    if Cancer_type[-4:] == 'GOLD':
        Omic_data = pd.read_csv(Omic_dir + Cancer_type[:4] + '_' + Data_type + '.csv', index_col=0)  # Load GOLD data
    else:
        Omic_data = pd.read_csv(Omic_dir + Cancer_type + '_' + Data_type + '.csv', index_col=0)  # Load standard data

    Omic_data = Omic_data.transpose()  # Transpose the omics data for analysis
    Omic_data = Omic_data.rename(columns=lambda x: x.replace('.', '-'))  # Replace dots in column names with dashes
    print(len(Omic_data))  # Print the number of samples in the omics data

    BP_data = pd.read_csv(BP_dir, index_col=0)  # Load pathway data
    Pathways = list(BP_data.index)  # Get the list of pathways
    miRNASymbols = []  # Initialize a list to hold miRNA symbols for pathways

    # Iterate over each pathway to collect miRNA symbols
    for i in Pathways:
        A = BP_data.loc[i]  # Get the miRNA symbols for the current pathway
        A = A.dropna()  # Remove any NaN values
        Temp_list = []  # Temporary list to hold miRNA symbols
        for value in A:
            Temp_list.append(value)  # Append each miRNA symbol to the temporary list
        miRNASymbols.append(Temp_list)  # Add the list of miRNA symbols to the main list

    miRNA_list = list(Omic_data.columns)  # Get all miRNA symbols from omics data
    Valid_Pathways = []  # List to hold valid pathways
    Valid_miRNASymbols = []  # List to hold valid miRNA symbols

    # Check each pathway for valid miRNAs present in omics data
    for i in range(0, len(miRNASymbols)):
        SameGene = set(miRNA_list) & set(miRNASymbols[i])  # Find common miRNAs between omics data and pathway
        SameGene = list(SameGene)  # Convert to list
        Gene_Num = len(SameGene)  # Count the number of common miRNAs
        Patial_Omic = Omic_data.loc[:, SameGene]  # Subset omics data to only include common miRNAs

        # Check for duplicates in the subsetted omics data
        if Patial_Omic.duplicated().any() or Patial_Omic.T.duplicated().any():
            continue  # Skip if duplicates are found

        if Gene_Num > 1:  # Only consider pathways with more than one miRNA
            Valid_Pathways.append(Pathways[i])  # Add valid pathway to the list
            Valid_miRNASymbols.append(miRNASymbols[i])  # Add corresponding miRNA symbols

    print(len(Valid_Pathways))  # Print the number of valid pathways
    print(len(Valid_miRNASymbols))  # Print the number of valid miRNA symbols

    pca = PCA()  # Initialize PCA for dimensionality reduction
    k = 5  # Define the number of clusters for KMeans
    d = round(len(Omic_data) * 0.05)  # Determine number of eigenvectors for PCA
    labels_list = []  # List to hold cluster labels from KMeans
    Sim_Matrix_list = []  # List to hold similarity matrices
    SigPathways_list = Valid_Pathways  # List of significant pathways
    SigGeneSymbols_list = Valid_miRNASymbols  # List of significant miRNA symbols
    X = SigGeneSymbols_list  # Assign significant miRNA symbols to X for processing
    km_pca = KMeans(n_clusters=k, init='k-means++', max_iter=100, n_init=10,
                    random_state=10)  # Initialize KMeans clustering

    # Perform clustering on each set of significant miRNA symbols
    for i in range(0, len(X)):
        SameGene = set(miRNA_list) & set(X[i])  # Find common miRNAs
        SameGene = list(SameGene)  # Convert to list
        Patial_Omic = Omic_data.loc[:, SameGene]  # Subset omics data
        Eucli_Matrix = Eucli_Distance(Patial_Omic)  # Calculate Euclidean distance matrix
        X_pca = pca.fit_transform(Eucli_Matrix)  # Apply PCA to the distance matrix
        labels_pca = km_pca.fit_predict(X_pca[:, :d])  # Fit KMeans and predict cluster labels
        labels_list.append(labels_pca)  # Append labels to the list

    Co_Matrix, Sim_Matrix_list = CSPA(labels_list)  # Compute consensus matrix using CSPA
    sc = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=10)  # Initialize Spectral Clustering

    # Iterate to refine clustering results
    for i in range(6):
        True_labels = sc.fit_predict(Co_Matrix)  # Fit Spectral Clustering and get true labels
        Pathway_Scores = []  # List to hold ARI scores for pathways
        for i in range(len(labels_list)):
            TempARI = adjusted_rand_score(True_labels, labels_list[i])  # Calculate ARI score
            Pathway_Scores.append(TempARI)  # Append score to the list

        Nums = int((len(labels_list)) * 0.6)  # Determine number of top pathways to keep
        top_200 = np.argsort(Pathway_Scores)[::-1][:Nums]  # Get indices of top pathways based on scores

        # Filter matrices and lists to keep only top pathways
        Sim_Matrix_list = np.take(Sim_Matrix_list, top_200, axis=0)  # Filter similarity matrices
        SigPathways_list = np.take(SigPathways_list, top_200, axis=0)  # Filter significant pathways
        labels_list = np.take(labels_list, top_200, axis=0)  # Filter cluster labels
        Co_Matrix = np.mean(Sim_Matrix_list, axis=0)  # Update consensus matrix
        print(f'Current_length:{Nums}')  # Print the current number of top pathways

    final_labels = sc.fit_predict(Co_Matrix)  # Final clustering using the updated consensus matrix

    final_scores = []  # List to hold final ARI scores
    for i in range(len(labels_list)):
        label_vector = labels_list[i]  # Get the cluster labels for the current pathway
        final_score = adjusted_rand_score(final_labels, label_vector)  # Calculate ARI score
        final_scores.append(final_score)  # Append score to the list

    sorted_index = np.argsort(final_scores)[::-1]  # Get indices of final scores sorted in descending order

    sorted_pathways = []  # List to hold sorted pathways
    Ffinal_scores = []  # List to hold final scores in sorted order

    # Collect sorted pathways and their corresponding scores
    for i in sorted_index:
        pathway = SigPathways_list[i]  # Get the pathway at the current index
        sorted_pathways.append(pathway)  # Append to sorted pathways
        Ffinal_scores.append(final_scores[i])  # Append corresponding score

    # Create a DataFrame for output
    df = pd.DataFrame(sorted_pathways, columns=["Pathway"])  # Create DataFrame with pathways
    df["Score"] = Ffinal_scores  # Add scores to the DataFrame
    print(df)  # Print the DataFrame

    output_dir = '../Datas/Omics/Output_Omics/'  # Define output directory
    np.savetxt(output_dir + Cancer_type + '_' + Data_type + '_Matrix.csv', Co_Matrix,
               delimiter=',')  # Save consensus matrix
    df.to_excel(output_dir + Cancer_type + "_" + Data_type + "_" + "Pathway_Sorting.xlsx")  # Save DataFrame to Excel

    end_time = time.time()  # Record the end time
    print(f'Current total running time of the program is {end_time - start_time} seconds')  # Print total running time