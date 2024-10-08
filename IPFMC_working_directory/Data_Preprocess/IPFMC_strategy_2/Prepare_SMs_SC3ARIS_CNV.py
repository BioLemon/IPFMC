import pandas as pd  # Import the pandas library for data manipulation
import numpy as np  # Import the numpy library for numerical operations
from sklearn.cluster import SpectralClustering  # Import the SpectralClustering algorithm from scikit-learn
from sklearn.decomposition import PCA  # Import the PCA algorithm from scikit-learn
from sklearn.cluster import KMeans  # Import the KMeans algorithm from scikit-learn
from sklearn.metrics.pairwise import euclidean_distances  # Import the euclidean_distances function from scikit-learn
import time  # Import the time module for measuring execution time
import warnings  # Import the warnings module

warnings.filterwarnings("ignore")  # Ignore warnings


# Define a function to calculate the Euclidean distance matrix and apply exponential transformation
def Eucli_Distance(Patial_Omic):
    dist_matrix = euclidean_distances(Patial_Omic)  # Calculate the Euclidean distance matrix
    dist_matrix = np.exp(-dist_matrix / np.max(dist_matrix))  # Apply exponential transformation
    return dist_matrix


# Record the start time of the program
start_time = time.time()

# Define a list of cancer types to process
Cancer_type_list = ['BRCA_GOLD', 'COAD_GOLD']

# Iterate over each cancer type in the list
for Cancer_type in Cancer_type_list:
    # Set the data type to CNV (Copy Number Variation)
    Data_type = 'CNV'

    # Define the directories for omics data, pathway data, and gold standard data
    Omic_dir = '../Datas/Omics/Raw_Omics/' + Cancer_type + '/'
    BP_dir = '../Datas/Pathways/Pathway_Index.csv'
    Gold_dir = '../Datas/Omics/Raw_Omics/' + Cancer_type + '/' + Cancer_type[:4] + '_label.csv'

    # Load the omics and gold standard data based on the cancer type
    if Cancer_type[-4:] == 'GOLD':
        Omic_data = pd.read_csv(Omic_dir + Cancer_type[:4] + '_' + Data_type + '.csv', index_col=0)
        Gold_label = np.genfromtxt(Gold_dir, delimiter=',', dtype='str')
        Gold_label = Gold_label[1:]
        unique, inverse = np.unique(Gold_label, return_inverse=True)
    else:
        Omic_data = pd.read_csv(Omic_dir + Cancer_type + '_' + Data_type + '.csv', index_col=0)

    # Transpose the omics data
    Omic_data = Omic_data.transpose()
    print(len(Omic_data))  # Print the length of the omics data

    # Load the pathway data from the CSV file
    BP_data = pd.read_csv(BP_dir, index_col=0)
    Pathways = list(BP_data.index)
    GeneSymbols = []

    # Extract gene symbols for each pathway
    for i in Pathways:
        A = BP_data.loc[i]
        A = A.dropna()
        Temp_list = []
        for value in A:
            Temp_list.append(value)
        GeneSymbols.append(Temp_list)

    FullGene = list(Omic_data.columns)
    var_list = []
    Valid_Pathways = []
    Valid_GeneSymbols = []

    # Filter out pathways with less than 2 genes
    for i in range(0, len(GeneSymbols)):
        SameGene = set(FullGene) & set(GeneSymbols[i])
        SameGene = list(SameGene)
        Gene_Num = len(SameGene)
        Patial_Omic = Omic_data.loc[:, SameGene]
        if Gene_Num > 1:
            Valid_Pathways.append(Pathways[i])
            Valid_GeneSymbols.append(GeneSymbols[i])

    print(len(Valid_Pathways))
    print(len(Valid_GeneSymbols))
    pca = PCA()
    k = 5  # Number of clusters
    d = round(len(Omic_data) * 0.05)  # Number of eigenvectors
    labels_list = []
    Sim_Matrix_list = []
    SigPathways_list = Valid_Pathways
    SigGeneSymbols_list = Valid_GeneSymbols
    X = SigGeneSymbols_list
    km_pca = KMeans(n_clusters=k, init='k-means++', max_iter=100, n_init=10)

    # Calculate the similarity matrices for each pathway
    for i in range(0, len(X)):
        SameGene = set(FullGene) & set(X[i])
        SameGene = list(SameGene)
        Patial_Omic = Omic_data.loc[:, SameGene]
        Eucli_Matrix = Eucli_Distance(Patial_Omic)
        Sim_Matrix_list.append(Eucli_Matrix)

    # Calculate the consensus similarity matrix
    Co_Matrix = np.mean(Sim_Matrix_list, axis=0)

    # Perform spectral clustering on the consensus similarity matrix
    sc = SpectralClustering(n_clusters=k, affinity="precomputed")

    # Iteratively refine the pathway selection
    for i in range(6):
        True_Matrix = Co_Matrix.flatten()
        Pathway_Scores = []
        for i in range(len(Sim_Matrix_list)):
            Comp_Matrix = Sim_Matrix_list[i].flatten()
            Loss = np.sum(abs(Comp_Matrix - True_Matrix))
            Relative_Score = 100 / (np.sqrt(Loss))
            Pathway_Scores.append(Relative_Score)
        Nums = int((len(Sim_Matrix_list)) * 0.6)
        top_200 = np.argsort(Pathway_Scores)[::-1][:Nums]
        Sim_Matrix_list = np.take(Sim_Matrix_list, top_200, axis=0)
        SigPathways_list = np.take(SigPathways_list, top_200, axis=0)
        Co_Matrix = np.mean(Sim_Matrix_list, axis=0)
        print(f'Current_length:{Nums}')

    # Calculate the final pathway scores
    True_Matrix = Co_Matrix.flatten()
    final_scores = []
    for i in range(len(Sim_Matrix_list)):
        Comp_Matrix = Sim_Matrix_list[i].flatten()
        Loss = np.sum(abs(Comp_Matrix - True_Matrix))
        Relative_Score = 100 / (np.sqrt(Loss))
        final_scores.append(Relative_Score)

    # Sort the pathways based on their scores
    sorted_index = np.argsort(final_scores)[::-1]
    sorted_pathways = []
    Ffinal_scores = []
    for i in sorted_index:
        pathway = SigPathways_list[i]
        sorted_pathways.append(pathway)
        Ffinal_scores.append(final_scores[i])

    # Create a DataFrame with the sorted pathways and their scores
    df = pd.DataFrame(sorted_pathways, columns=["Pathway"])
    df["Score"] = Ffinal_scores

    # Define the output directory and save the results
    output_dir = '../Datas/Omics/Output_Omics/' + Cancer_type + '/'
    np.savetxt(output_dir + Cancer_type + '_' + Data_type + '_Matrix.csv', Co_Matrix, delimiter=',')
    df.to_excel(output_dir + Cancer_type + "_" + Data_type + "_" + "Pathway_Sorting.xlsx")

# Record the end time of the program and calculate the total running time
end_time = time.time()
print(f'The total running time of the program is {end_time - start_time} seconds')
