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
    euclidean_matrix = squareform(pdist(df, metric='euclidean'))
    #euclidean_matrix = np.exp(-euclidean_matrix / np.max(euclidean_matrix))
    return euclidean_matrix

def CSPA(labels_list):
    """
    Cluster-based Similarity Partitioning Algorithm (CSPA)
    Input: a list of arrays of cluster labels, each array has shape (N,)
    Output: a consensus matrix, shape (N, N)
    """
    epoch = 0
    N = len(labels_list[0]) # number of samples
    S_list = [] # a list of similarity matrices
    for labels in labels_list:
        S = np.zeros((N, N)) # initialize a similarity matrix, shape (N, N)
        for i in range(N):
            for j in range(N):
                if labels[i] == labels[j]: # if two samples have the same label, their similarity is 1
                    S[i][j] = 1
        S_list.append(S) # append the similarity matrix to the list
        epoch = epoch+1
    C = np.mean(S_list, axis=0) # compute the consensus matrix by averaging all similarity matrices
    return C, S_list

start_time = time.time()
Cancer_type_list = ['BRCA','COAD','KIRC','LUAD','LUSC','ACC','KIRP','LIHC','THYM']
for Cancer_type in Cancer_type_list:
    Data_type = 'miRNA'
    Omic_dir = '../../Datas/Omics/Raw_Omics/' + Cancer_type + '/'
    BP_dir = '../../Datas/Pathways/miRNA_Pathway_Index.csv'
    Gold_dir = '../../Datas/Omics/Raw_Omics/' + Cancer_type + '/'+Cancer_type[:4]+'_label.csv'
    if Cancer_type[-4:] == 'GOLD':
        Omic_data = pd.read_csv(Omic_dir + Cancer_type[:4] + '_' + Data_type + '.csv', index_col=0)
        Gold_label = np.genfromtxt(Gold_dir, delimiter=',', dtype='str')
        Gold_label = Gold_label[1:]
        unique, inverse = np.unique(Gold_label, return_inverse=True)
    else:
        Omic_data = pd.read_csv(Omic_dir+Cancer_type+'_'+Data_type+'.csv', index_col=0)
    Omic_data = Omic_data.transpose()
    Omic_data = Omic_data.rename(columns=lambda x: x.replace('.', '-'))
    print(len(Omic_data))

    BP_data = pd.read_csv(BP_dir, index_col=0)
    Pathways = list(BP_data.index)
    miRNASymbols = []

    for i in Pathways:
        A = BP_data.loc[i]
        A = A.dropna()
        Temp_list = []
        for value in A:
            Temp_list.append(value)
        miRNASymbols.append(Temp_list)

    miRNA_list = list(Omic_data.columns)
    Valid_Pathways = []
    Valid_miRNASymbols = []

    for i in range(0, len(miRNASymbols)):
        SameGene = set(miRNA_list) & set(miRNASymbols[i])
        SameGene = list(SameGene)
        Gene_Num = len(SameGene)
        Patial_Omic = Omic_data.loc[:, SameGene]
        if Patial_Omic.duplicated().any() or Patial_Omic.T.duplicated().any():
            continue
        if Gene_Num > 1:
            Valid_Pathways.append(Pathways[i])
            Valid_miRNASymbols.append(miRNASymbols[i])

    print(len(Valid_Pathways))
    print(len(Valid_miRNASymbols))
    pca = PCA()
    k = 5
    d = round(len(Omic_data) * 0.05)
    labels_list = []
    Sim_Matrix_list = []
    SigPathways_list = Valid_Pathways
    SigGeneSymbols_list = Valid_miRNASymbols
    X = SigGeneSymbols_list
    km_pca = KMeans(n_clusters=k, init='k-means++', max_iter=100, n_init=10, random_state=10)
    for i in range(0, len(X)):
        SameGene = set(miRNA_list) & set(X[i])
        SameGene = list(SameGene)
        Patial_Omic = Omic_data.loc[:, SameGene]
        Eucli_Matrix = Eucli_Distance(Patial_Omic)
        X_pca = pca.fit_transform(Eucli_Matrix)
        labels_pca = km_pca.fit_predict(X_pca[:, :d])
        labels_list.append(labels_pca)
    Co_Matrix, Sim_Matrix_list = CSPA(labels_list)
    sc = SpectralClustering(n_clusters=k, affinity="precomputed",random_state=10)

    for i in range(6):
        True_labels = sc.fit_predict(Co_Matrix)
        Pathway_Scores = []
        for i in range(len(labels_list)):
            TempARI = adjusted_rand_score(True_labels, labels_list[i])
            Pathway_Scores.append(TempARI)
        Nums = int((len(labels_list)) * 0.6)
        top_200 = np.argsort(Pathway_Scores)[::-1][:Nums]

        Sim_Matrix_list = np.take(Sim_Matrix_list, top_200, axis=0)
        SigPathways_list = np.take(SigPathways_list, top_200, axis=0)
        labels_list = np.take(labels_list, top_200, axis=0)
        Co_Matrix = np.mean(Sim_Matrix_list, axis=0)
        print(f'Current_length:{Nums}')

    final_labels = sc.fit_predict(Co_Matrix)

    final_scores = []
    for i in range(len(labels_list)):
        label_vector = labels_list[i]
        final_score = adjusted_rand_score(final_labels, label_vector)
        final_scores.append(final_score)

    sorted_index = np.argsort(final_scores)[::-1]

    sorted_pathways = []

    Ffinal_scores = []
    for i in sorted_index:
        pathway = SigPathways_list[i]
        sorted_pathways.append(pathway)
        Ffinal_scores.append(final_scores[i])

    df = pd.DataFrame(sorted_pathways, columns=["Pathway"])
    df["Score"] = Ffinal_scores
    print(df)


    output_dir = '../../Datas/Omics/Output_Omics/'
    np.savetxt(output_dir+Cancer_type+'_'+Data_type+'_Matrix.csv', Co_Matrix, delimiter=',')
    df.to_excel(output_dir+Cancer_type+"_"+Data_type+"_"+"Pathway_Sorting.xlsx")
    end_time = time.time()
    print(f'The total running time of the program is {end_time-start_time} seconds')
