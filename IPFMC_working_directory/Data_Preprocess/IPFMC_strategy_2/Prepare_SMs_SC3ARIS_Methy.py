import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances
import time
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")


# Function to calculate the Euclidean distance matrix and convert it to a similarity matrix
def Eucli_Distance(Patial_Omic):
    # Compute the Euclidean distance matrix for the given data
    dist_matrix = euclidean_distances(Patial_Omic)
    # Convert the distance matrix to a similarity matrix using an exponential decay function
    dist_matrix = np.exp(-dist_matrix / np.max(dist_matrix))
    return dist_matrix


# Start timing the execution of the program
start_time = time.time()

# List of cancer types to analyze
Cancer_type_list = ['BRCA_GOLD', 'COAD_GOLD']

# Loop through each cancer type in the list
for Cancer_type in Cancer_type_list:
    Data_type = 'Methy'  # Define the type of omics data to be used
    # Define file paths for omics data, pathway data, and gold standard labels
    Omic_dir = '../Datas/Omics/Raw_Omics/' + Cancer_type + '/'
    BP_dir = '../Datas/Pathways/Pathway_Index.csv'
    Gold_dir = '../Datas/Omics/Raw_Omics/' + Cancer_type + '/' + Cancer_type[:4] + '_label.csv'

    # Check if the cancer type ends with 'GOLD' to load the appropriate data
    if Cancer_type[-4:] == 'GOLD':
        # Load the omics data from the specified CSV file
        Omic_data = pd.read_csv(Omic_dir + Cancer_type[:4] + '_' + Data_type + '.csv', index_col=0)
        # Load the gold standard labels from the specified CSV file
        Gold_label = np.genfromtxt(Gold_dir, delimiter=',', dtype='str')
        Gold_label = Gold_label[1:]  # Exclude the header
        # Get unique labels and their inverse mapping
        unique, inverse = np.unique(Gold_label, return_inverse=True)
    else:
        # If not 'GOLD', load the omics data from a different file
        Omic_data = pd.read_csv(Omic_dir + Cancer_type + '_' + Data_type + '.csv', index_col=0)

    # Transpose the omics data for easier manipulation
    Omic_data = Omic_data.transpose()
    print(len(Omic_data))  # Print the number of samples in the omics data

    # Load the pathway data from the specified CSV file
    BP_data = pd.read_csv(BP_dir, index_col=0)
    Pathways = list(BP_data.index)  # Get the list of pathways
    GeneSymbols = []  # Initialize a list to hold gene symbols for each pathway

    # Loop through each pathway to extract gene symbols
    for i in Pathways:
        A = BP_data.loc[i]  # Get the gene symbols for the current pathway
        A = A.dropna()  # Remove any NaN values
        Temp_list = []  # Temporary list to hold gene symbols
        for value in A:
            Temp_list.append(value)  # Append each gene symbol to the temporary list
        GeneSymbols.append(Temp_list)  # Add the list of gene symbols to the main list

    FullGene = list(Omic_data.columns)  # Get all gene names from the omics data
    var_list = []  # Initialize an empty list for variance (not used in this code)
    Valid_Pathways = []  # List to store valid pathways with sufficient genes
    Valid_GeneSymbols = []  # List to store gene symbols corresponding to valid pathways

    # Filter pathways to retain only those with more than one gene present in the omics data
    for i in range(0, len(GeneSymbols)):
        SameGene = set(FullGene) & set(GeneSymbols[i])  # Find common genes between omics data and current pathway
        SameGene = list(SameGene)  # Convert the set to a list
        Gene_Num = len(SameGene)  # Count the number of common genes
        Patial_Omic = Omic_data.loc[:, SameGene]  # Subset the omics data to include only the common genes
        if Gene_Num > 1:  # Check if there are more than one common gene
            Valid_Pathways.append(Pathways[i])  # Add the pathway to valid pathways
            Valid_GeneSymbols.append(GeneSymbols[i])  # Add the corresponding gene symbols

    print(len(Valid_Pathways))  # Print the number of valid pathways
    print(len(Valid_GeneSymbols))  # Print the number of valid gene symbols

    # Initialize PCA for dimensionality reduction (not used in this code)
    pca = PCA()
    k = 5  # Define the number of clusters for KMeans and Spectral Clustering
    d = round(len(Omic_data) * 0.05)  # Define the number of eigenvectors (5% of samples)

    labels_list = []  # List to hold labels (not used in this code)
    Sim_Matrix_list = []  # List to store similarity matrices for pathways
    SigPathways_list = Valid_Pathways  # List of significant pathways
    SigGeneSymbols_list = Valid_GeneSymbols  # List of significant gene symbols
    X = SigGeneSymbols_list  # Assign significant gene symbols to X for processing

    # Initialize KMeans clustering
    km_pca = KMeans(n_clusters=k, init='k-means++', max_iter=100, n_init=10)

    # Loop through each set of significant gene symbols
    for i in range(0, len(X)):
        SameGene = set(FullGene) & set(X[i])  # Find common genes
        SameGene = list(SameGene)  # Convert to list
        Patial_Omic = Omic_data.loc[:, SameGene]  # Subset omics data
        Eucli_Matrix = Eucli_Distance(Patial_Omic)  # Calculate the Euclidean distance matrix
        Sim_Matrix_list.append(Eucli_Matrix)  # Append the similarity matrix to the list

    # Calculate the co-association matrix by averaging the similarity matrices
    Co_Matrix = np.mean(Sim_Matrix_list, axis=0)

    # Initialize Spectral Clustering with the specified number of clusters
    sc = SpectralClustering(n_clusters=k, affinity="precomputed")

    # Iteratively refine the similarity matrices and pathway scores
    for i in range(6):  # Perform 6 iterations
        True_Matrix = Co_Matrix.flatten()  # Flatten the co-association matrix for comparison
        Pathway_Scores = []  # List to store scores for each pathway

        # Calculate scores based on the similarity of each matrix to the co-association matrix
        for i in range(len(Sim_Matrix_list)):
            Comp_Matrix = Sim_Matrix_list[i].flatten()  # Flatten the current similarity matrix
            Loss = np.sum(abs(Comp_Matrix - True_Matrix))  # Calculate the loss as the absolute difference
            Relative_Score = 100 / (np.sqrt(Loss))  # Calculate a relative score based on the loss
            Pathway_Scores.append(Relative_Score)  # Append the score to the list

        # Select the top pathways based on their scores
        Nums = int((len(Sim_Matrix_list)) * 0.6)  # Determine the number of top pathways to keep (60%)
        top_200 = np.argsort(Pathway_Scores)[::-1][:Nums]  # Get the indices of the top pathways
        Sim_Matrix_list = np.take(Sim_Matrix_list, top_200, axis=0)  # Keep only the top similarity matrices
        SigPathways_list = np.take(SigPathways_list, top_200, axis=0)  # Keep only the top pathways
        Co_Matrix = np.mean(Sim_Matrix_list, axis=0)  # Recalculate the co-association matrix
        print(f'Current_length: {Nums}')  # Print the current number of pathways

    # Final computation of scores for the remaining pathways
    True_Matrix = Co_Matrix.flatten()  # Flatten the co-association matrix
    final_scores = []  # List to store final scores for pathways

    # Calculate final scores for each pathway
    for i in range(len(Sim_Matrix_list)):
        Comp_Matrix = Sim_Matrix_list[i].flatten()  # Flatten the current similarity matrix
        Loss = np.sum(abs(Comp_Matrix - True_Matrix))  # Calculate the loss
        Relative_Score = 100 / (np.sqrt(Loss))  # Calculate the relative score
        final_scores.append(Relative_Score)  # Append the score to the list

    # Sort pathways based on final scores
    sorted_index = np.argsort(final_scores)[::-1]  # Get indices of pathways sorted by score in descending order
    sorted_pathways = []  # List to hold sorted pathways

    Ffinal_scores = []  # List to hold final scores of sorted pathways
    for i in sorted_index:
        pathway = SigPathways_list[i]  # Get the pathway corresponding to the sorted index
        sorted_pathways.append(pathway)  # Append to the sorted pathways list
        Ffinal_scores.append(final_scores[i])  # Append the corresponding score

    # Create a DataFrame to hold the sorted pathways and their scores
    df = pd.DataFrame(sorted_pathways, columns=["Pathway"])
    df["Score"] = Ffinal_scores  # Add scores to the DataFrame

    # Define output directory for saving results
    output_dir = '../Datas/Omics/Output_Omics/' + Cancer_type + '/'
    # Save the co-association matrix to a CSV file
    np.savetxt(output_dir + Cancer_type + '_' + Data_type + '_Matrix.csv', Co_Matrix, delimiter=',')
    # Save the sorted pathways and their scores to an Excel file
    df.to_excel(output_dir + Cancer_type + "_" + Data_type + "_" + "Pathway_Sorting.xlsx")

    # End timing and print the total running time of the program
    end_time = time.time()
    print(f'The total running time of the program is {end_time - start_time} seconds')
    