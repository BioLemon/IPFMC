import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.cluster import SpectralClustering
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import warnings
warnings.filterwarnings("ignore")

def CSPA(labels_list):
    """
    Cluster-based Similarity Partitioning Algorithm (CSPA)
    Input: a list of arrays of cluster labels, each array has shape (N,)
    Output: a consensus matrix, shape (N, N)
    """
    epoch = 0
    N = len(labels_list[0])  # number of samples
    S_list = []  # a list of similarity matrices
    for labels in labels_list:
        S = np.zeros((N, N))  # initialize a similarity matrix, shape (N, N)
        for i in range(N):
            for j in range(N):
                if labels[i] == labels[j]:  # if two samples have the same label, their similarity is 1
                    S[i][j] = 1
        S_list.append(S)  # append the similarity matrix to the list
        epoch = epoch + 1
        print(f'Calculating Binary similarity metrics of current dataset, progress {(epoch/len(labels_list))*100:.2f}%')
    C = np.mean(S_list, axis=0)  # compute the consensus matrix by averaging all similarity matrices
    return C, S_list



def ipfmc_discretize(dataset, pathwayinfo, k=5, fusetime=6, proportion=0.6, seed=10):
    """
        :param dataset: Your omics dataset.
        :param pathwayinfo: Pathways and their containing genetic information.
        :param k: The number of initial points of kmeans clustering
        :param fusetime: Number of pathway screening and fusion performed
        :param proportion: The proportion of pathways that were retained was fused at each iteration
        :param seed: Random number seed, set to None if no seed is needed
        :return: Final representation of input dataset; pathway rankings of input dataset.
    """
    BP_data = pathwayinfo
    Omic_data = dataset
    Omic_data = Omic_data.transpose()
    Pathways = list(BP_data.index)  # List of all pathway names
    GeneSymbols = []  # List of all pathway gene sets
    # Import the pathway gene set into GeneSymbols in the correct way
    for i in Pathways:
        A = BP_data.loc[i]
        A = A.dropna()
        Temp_list = []
        for value in A:
            Temp_list.append(value)
        GeneSymbols.append(Temp_list)

    # Get all gene names in the omics dataset, saved as FullGene
    FullGene = list(Omic_data.columns)
    Valid_Pathways = []  # Used to store all valid pathways (pathways contained in our dataset)
    Valid_GeneSymbols = []  # Used to store valid pathway gene sets

    # Loop through all pathways to find all valid pathways
    for i in range(0, len(GeneSymbols)):
        SameGene = set(FullGene) & set(GeneSymbols[i])
        SameGene = list(SameGene)
        Gene_Num = len(SameGene)
        Patial_Omic = Omic_data.loc[:, SameGene]
        # Check if Partial_Omic has duplicate rows or columns
        if Patial_Omic.duplicated().any() or Patial_Omic.T.duplicated().any():
            continue  # Skip this pathway
        if Gene_Num > 1:  # If statement checks if the list is not empty, if not empty, it returns True, if statement executes
            Valid_Pathways.append(Pathways[i])
            Valid_GeneSymbols.append(GeneSymbols[i])
    print(f"The effective pathway length of the current dataset is {len(Valid_Pathways)}")
    pca = PCA()
    d = round(len(Omic_data) * 0.05)  # Number of eigenvectors
    labels_list = []  # Used to store all clustering label results
    Sim_Matrix_list = []
    SigPathways_list = Valid_Pathways
    SigGeneSymbols_list = Valid_GeneSymbols
    X = SigGeneSymbols_list
    km_pca = KMeans(n_clusters=k, init='k-means++', max_iter=100, n_init=10, random_state=seed)
    for i in range(0, len(X)):
        SameGene = set(FullGene) & set(X[i])
        SameGene = list(SameGene)
        Patial_Omic = Omic_data.loc[:, SameGene]
        Eucli_Matrix = squareform(pdist(Patial_Omic, metric='euclidean'))
        X_pca = pca.fit_transform(Eucli_Matrix)
        labels_pca = km_pca.fit_predict(X_pca[:, :d])
        labels_list.append(labels_pca)
        print(f"Each pathway representation is being calculated, current progress is {((i+1)/len(X))*100:.2f}%")
    Co_Matrix, Sim_Matrix_list = CSPA(labels_list)
    sc = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=seed)

    # At this time, Sim_Matrix_list, labels_list, SigPathways_list are of the same length
    for i in range(fusetime):
        print(f'Interative Pathway Fusing, current epoch {i}')
        True_labels = sc.fit_predict(Co_Matrix)  # Calculate reference labels
        Pathway_Scores = []
        for i in range(len(labels_list)):
            TempARI = adjusted_rand_score(True_labels, labels_list[i])
            Pathway_Scores.append(TempARI)
        Nums = int((len(labels_list)) * proportion)  # Take the index of the first half and store it in the top_200 variable
        top_200 = np.argsort(Pathway_Scores)[::-1][:Nums]  # Use numpy's argsort function to sort Pathway_Scores in descending order and return the index

        # Use numpy's take function to take the similarity matrix corresponding to top_200 from Sim_Matrix_list
        # Store it in the top_200_sim_matrix variable, which is a three-dimensional array
        Sim_Matrix_list = np.take(Sim_Matrix_list, top_200, axis=0)
        SigPathways_list = np.take(SigPathways_list, top_200, axis=0)
        labels_list = np.take(labels_list, top_200, axis=0)

        # Use numpy's mean function to average top_200_sim_matrix along the first axis
        # Get the final similarity matrix, overwrite Co_Matrix
        Co_Matrix = np.mean(Sim_Matrix_list, axis=0)
        # Codis_Matrix = np.mean(Dis_Matrix_list, axis=0)
        print(f'Length of remaining pathways: {Nums}')

    # Co_Matrix = Codis_Matrix
    # Perform spectral clustering on Co_Matrix again to get the final clustering result
    final_labels = sc.fit_predict(Co_Matrix)

    # Define an empty list to store the ARI scores of each label with the final clustering result
    final_scores = []
    # Traverse each index in top_200
    for i in range(len(labels_list)):
        # Take the corresponding label vector from labels_list
        label_vector = labels_list[i]
        # Calculate the ARI score with the final clustering result and append it to the final_scores list
        final_score = adjusted_rand_score(final_labels, label_vector)
        final_scores.append(final_score)

    # Use numpy's argsort function to sort final_scores in descending order and return the index
    # These indexes are relative to top_200, so you need to use top_200 to index the original Sim_Matrix_list and SigPathways_list
    # Store in sorted_index variable
    sorted_index = np.argsort(final_scores)[::-1]

    # Define an empty list to store the sorted pathways
    sorted_pathways = []

    Ffinal_scores = []
    # Traverse each index in sorted_index
    for i in sorted_index:
        pathway = SigPathways_list[i]
        sorted_pathways.append(pathway)
        Ffinal_scores.append(final_scores[i])

    # Use pandas library to convert sorted_pathways into a DataFrame object and add a column as final_scores
    df = pd.DataFrame(sorted_pathways, columns=["Pathway"])
    df["Score"] = Ffinal_scores
    return Co_Matrix, df


def ipfmc_average(dataset, pathwayinfo, fusetime=6, proportion=0.6):
    """
        :param dataset: Your omics dataset.
        :param pathwayinfo: Pathways and their containing genetic information.
        :param fusetime: Number of pathway screening and fusion performed
        :param proportion: The proportion of pathways that were retained was fused at each iteration
        :return: Final representation of input dataset; pathway rankings of input dataset.
    """
    BP_data = pathwayinfo
    Omic_data = dataset
    Omic_data = Omic_data.transpose()
    Pathways = list(BP_data.index)  # List of all pathway names
    GeneSymbols = []  # List of all pathway gene sets
    # Import the pathway gene set into GeneSymbols in the correct way
    for i in Pathways:
        A = BP_data.loc[i]
        A = A.dropna()
        Temp_list = []
        for value in A:
            Temp_list.append(value)
        GeneSymbols.append(Temp_list)

    # Get all gene names in the omics dataset, saved as FullGene
    FullGene = list(Omic_data.columns)
    Valid_Pathways = []  # Used to store all valid pathways (pathways contained in our dataset)
    Valid_GeneSymbols = []  # Used to store valid pathway gene sets

    # Loop through all pathways to find all valid pathways
    for i in range(0, len(GeneSymbols)):
        SameGene = set(FullGene) & set(GeneSymbols[i])
        SameGene = list(SameGene)
        Gene_Num = len(SameGene)
        Patial_Omic = Omic_data.loc[:, SameGene]
        # Check if Partial_Omic has duplicate rows or columns
        if Patial_Omic.duplicated().any() or Patial_Omic.T.duplicated().any():
            continue  # Skip this pathway
        if Gene_Num > 1:  # If statement checks if the list is not empty, if not empty, it returns True, if statement executes
            Valid_Pathways.append(Pathways[i])
            Valid_GeneSymbols.append(GeneSymbols[i])
    print(f"The effective pathway length of the current dataset is {len(Valid_Pathways)}")
    Sim_Matrix_list = []
    SigPathways_list = Valid_Pathways
    SigGeneSymbols_list = Valid_GeneSymbols
    X = SigGeneSymbols_list
    for i in range(0, len(X)):
        SameGene = set(FullGene) & set(X[i])
        SameGene = list(SameGene)
        Patial_Omic = Omic_data.loc[:, SameGene]
        Eucli_Matrix = squareform(pdist(Patial_Omic, metric='euclidean'))
        Eucli_Matrix = np.exp(-Eucli_Matrix / np.max(Eucli_Matrix))
        Sim_Matrix_list.append(Eucli_Matrix)
        print(f"Calculating pathway representations, current progress {((i + 1) / len(X)) * 100:.2f}%")
    Co_Matrix = np.mean(Sim_Matrix_list, axis=0)
    # At this time, Sim_Matrix_list, labels_list, SigPathways_list are of the same length
    for i in range(fusetime):
        True_Matrix = Co_Matrix.flatten()
        Pathway_Scores = []
        for i in range(len(Sim_Matrix_list)):
            Comp_Matrix = Sim_Matrix_list[i].flatten()
            # Loss = np.linalg.norm(Comp_Matrix-True_Matrix)
            Loss = np.sum(abs(Comp_Matrix - True_Matrix))
            Relative_Score = 100 / (np.sqrt(Loss))
            Pathway_Scores.append(Relative_Score)
        Nums = int((len(Sim_Matrix_list)) * proportion)  # Take the index of the first half and store it in the top_200 variable
        top_200 = np.argsort(Pathway_Scores)[::-1][:Nums]  # Use numpy's argsort function to sort Pathway_Scores in descending order and return the index
        # Use numpy's take function to take the similarity matrix corresponding to top_200 from Sim_Matrix_list
        # Store it in the top_200_sim_matrix variable, which is a three-dimensional array
        Sim_Matrix_list = np.take(Sim_Matrix_list, top_200, axis=0)  # Similarity matrix
        # Dis_Matrix_list = np.take(Dis_Matrix_list, top_200, axis=0)  #
        SigPathways_list = np.take(SigPathways_list, top_200, axis=0)
        # Use numpy's mean function to average top_200_sim_matrix along the first axis
        # Get the final similarity matrix, overwrite Co_Matrix
        Co_Matrix = np.mean(Sim_Matrix_list, axis=0)
        # Codis_Matrix = np.mean(Dis_Matrix_list, axis=0)
        print(f'Length of remaining pathways: {Nums}')

    True_Matrix = Co_Matrix.flatten()
    # Define an empty list to store the ARI scores of each label with the final clustering result
    final_scores = []
    # Traverse each index in top_200
    for i in range(len(Sim_Matrix_list)):
        Comp_Matrix = Sim_Matrix_list[i].flatten()
        Loss = np.sum(abs(Comp_Matrix - True_Matrix))
        Relative_Score = 100 / (np.sqrt(Loss))
        final_scores.append(Relative_Score)

    # These indexes are relative to top_200, so you need to use top_200 to index the original Sim_Matrix_list and SigPathways_list
    # Store in sorted_index variable
    sorted_index = np.argsort(final_scores)[::-1]

    # Define an empty list to store the sorted pathways
    sorted_pathways = []

    Ffinal_scores = []
    # Traverse each index in sorted_index
    for i in sorted_index:
        pathway = SigPathways_list[i]
        sorted_pathways.append(pathway)
        Ffinal_scores.append(final_scores[i])

    # Use pandas library to convert sorted_pathways into a DataFrame object and add a column as final_scores
    df = pd.DataFrame(sorted_pathways, columns=["Pathway"])
    df["Score"] = Ffinal_scores
    return Co_Matrix, df



def spec_cluster(dataset, fusion_matrix, k):
    """
        :param dataset: one of your original data sets, it is used to get sample names.
        :param fusion_matrix: the representation matrix obtained by ipfmc.
        :param k: number of clusters.
        :return: cluster labels with sample names as row index.
    """
    sc = SpectralClustering(n_clusters=k, affinity="precomputed")
    labels = sc.fit_predict(fusion_matrix)
    dataset = dataset.transpose()
    ID = dataset.index
    labels = pd.DataFrame(labels, index=ID)
    return labels

def suggest_k(represent,min_k=2,max_k=8):
    """
        :param represent: the representation matrix obtained by ipfmc.
        :param min_k: minimum number of clusters.
        :param max_k: maximum number of clusters.
        :return: suggested number of clusters.
    """
    silhouette_list = []
    for k in range(min_k, max_k+1):
        sc = SpectralClustering(n_clusters=k, affinity="precomputed")
        labels = sc.fit_predict(represent)
        labels = pd.DataFrame(labels)
        silhouette = metrics.silhouette_score(represent, labels)
        silhouette_list.append(silhouette)
    max_value = max(silhouette_list)
    # 使用index方法找到最大值的索引
    Best_K = silhouette_list.index(max_value)+2
    return Best_K












