import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import warnings

warnings.filterwarnings("ignore")
"""
BP_dir = '/home/zhanghaoy/Pytorchdeep/Datas/Pathways/'
#Omic_dir = '/home/zhanghaoy/Pytorchdeep/Datas/Omics/Raw_Omics/'+Cancer_type+'/'
BioProcesses = pd.read_json(BP_dir+'GoTerm_NCP.json')
"""


def Create_pathwayindex(GoTerminfo, min_genes, max_genes):
    '''
    :param GoTerminfo: Pathway information in json format downloaded from GSEA.
    :param min_genes: Set the minimum number of genes in the pathway index created
    :param max_genes: Set the maximum number of genes in the pathway index created
    :return: pathway_index
    '''
    Process_lists = []
    Pronames_lists = []
    print(GoTerminfo)
    for i in GoTerminfo:
        a = GoTerminfo[i]['geneSymbols']
        if (len(a) > min_genes and len(a) < max_genes):
            Process_lists.append(a)
            Pronames_lists.append(i)
    Process_lists = pd.DataFrame(Process_lists)
    Pronames_lists = pd.DataFrame(Pronames_lists)
    Pronames_lists = Pronames_lists.rename(columns={0: 'Pathways'})
    A = [Pronames_lists, Process_lists]
    A = pd.concat(A, axis=1)
    A.set_index(A.columns[0], inplace=True)
    return A

def Create_miRNA_pathwayindex(miRNA_Genes_dict, pathwayinfo):
    # Only read the pathway name and gene name columns
    BP_data = pathwayinfo
    Pathways = list(BP_data.index)
    GeneSymbols = []
    # Import the pathway gene set into GeneSymbols in the correct way
    for i in Pathways:
        A = BP_data.loc[i]
        A = A.dropna()
        Temp_list = []
        for value in A:
            Temp_list.append(value)
        GeneSymbols.append(Temp_list)

    Valid_Pathways = []  # Used to store all valid pathways (pathways included in our dataset)
    Valid_miRNASymbols = []  # Used to store the miRNA set of valid pathways

    # Convert the miRNA target gene list to a two-dimensional array
    miRNA_Genes_array = np.array(list(miRNA_Genes_dict.values()))
    print(len(miRNA_Genes_array))
    print(miRNA_Genes_array)
    # Convert the miRNA name list to a one-dimensional array
    miRNA_Names_array = np.array(list(miRNA_Genes_dict.keys()))
    # Get the total number of elements in the population
    N = 20000

    # Traverse all pathways
    for i in range(len(GeneSymbols)):
        # Get the current pathway gene set
        Pathway_Genes = GeneSymbols[i]
        # Get the number of genes in the current pathway
        M = len(Pathway_Genes)
        # Calculate the number of intersections between each miRNA target gene set and the current pathway gene set, return an integer array
        intersection_count = np.array(
            [len(set(Pathway_Genes) & set(Target_Genes)) for Target_Genes in miRNA_Genes_array])
        # Calculate the number of each miRNA target gene set, return an integer array
        target_count = np.array([len(Target_Genes) for Target_Genes in miRNA_Genes_array])
        # Calculate the p-value of the hypergeometric distribution test of each miRNA and the current pathway, return a floating-point number array
        p_values = hypergeom.sf(intersection_count, N, M, target_count)
        # Determine whether the p-value is less than 0.01, return a boolean array
        is_significant = p_values < 0.005
        # If any p-value is less than 0.05, it means that the pathway is valid, add it to the valid pathway list, and add the corresponding miRNA to the valid miRNA list
        if is_significant.any():
            Valid_Pathways.append(Pathways[i])
            Valid_miRNASymbols.append(miRNA_Names_array[is_significant])
        print(
            f'Current progress:{i}/{len(GeneSymbols)}，Number of mirnas in the current pathway:{len(miRNA_Names_array[is_significant])}')
    print(len(Valid_Pathways))
    print(len(Valid_miRNASymbols))
    Valid_Pathways = pd.DataFrame(Valid_Pathways)
    Valid_miRNASymbols = pd.DataFrame(Valid_miRNASymbols)
    Valid_Pathways = Valid_Pathways.rename(columns={0: 'Pathways'})
    A = [Valid_Pathways, Valid_miRNASymbols]
    A = pd.concat(A, axis=1)
    A.set_index(A.columns[0], inplace=True)
    print('Successfully created miRNA-pathway index')
    return A


def gene_occurrence(Full_Genes, Sel_Pathways, Full_Pathways, top=50):
    Sel_Pathways = list(Sel_Pathways['Pathway'])  # Get all the pathway names after screening, and overwrite the original Sel_Pathways variable
    Parcial_Pathways = Full_Pathways.loc[Sel_Pathways]
    print(Parcial_Pathways)
    # Define a dictionary to store the count of each element
    count_dict = {}
    # Traverse the first list and initialize the count of each element to 0
    for item in Full_Genes:
        count_dict[item] = 0
    # Traverse the second list, if the element is in the first list, add 1 to the corresponding count
    for Pathway in Sel_Pathways:
        Genes = Parcial_Pathways.loc[Pathway]
        Genes = list(Genes.dropna())
        for i in Genes:
            if i in Full_Genes:
                count_dict[i] += 1
    Gene_Count = pd.DataFrame.from_dict(count_dict, orient='index', columns=['count'])
    Gene_Count = Gene_Count.sort_values(by='count', ascending=False)
    Gene_Count = Gene_Count.iloc[:top]
    return Gene_Count

def create_mirTar_Dict(MTI_info, mir_set):
    Let_list = []
    Mir_list = []
    All_Names = mir_set
    target_gene_data = MTI_info
    LET_PREFIX = 'hsa-let'
    MIR_PREFIX = 'hsa-mir'
    PREFIX_LEN = 7
    for miRNA in All_Names:
        if miRNA[:PREFIX_LEN] == LET_PREFIX:
            Let_list.append(miRNA)
        elif miRNA[:PREFIX_LEN] == MIR_PREFIX:
            Mir_list.append(miRNA)
    # print(Let_list)
    # Convert our target gene to a dictionary
    target_gene_dict = {}
    for i in range(len(target_gene_data['miRNA'])):
        MiR = target_gene_data['miRNA'][i]
        Gene = target_gene_data['Target Gene'][i]
        if MiR not in target_gene_dict:
            target_gene_dict[MiR] = set()
        target_gene_dict[MiR].add(Gene)
    print(len(target_gene_dict))
    print(target_gene_dict['hsa-let-7a-5p'])
    print(len(target_gene_dict['hsa-let-7a-5p']))

    miRNA_Genes_dict = {}
    for miRNA in Let_list:
        Genes = set()
        if len(miRNA) == 10:  # Here we extract miRNAs similar to hsa-let-7c
            # Directly find the corresponding target gene set from the dictionary and take the union
            Genes = Genes.union(target_gene_dict.get(miRNA + '-5p', set()))
            Genes = Genes.union(target_gene_dict.get(miRNA + '-3p', set()))
            Genes = Genes.union(target_gene_dict.get(miRNA, set()))
        elif miRNA[10:] == '-1':
            # Directly find the corresponding target gene set from the dictionary and take the union
            Genes = Genes.union(target_gene_dict.get(miRNA[:10] + '-5p', set()))
            Genes = Genes.union(target_gene_dict.get(miRNA[:10] + '-3p', set()))
            Genes = Genes.union(target_gene_dict.get(miRNA, set()))
        else:  # Other RNAs, that is, miRNAs similar to hsa-let-7a-2
            # Directly find the corresponding target gene set from the dictionary and take the union
            Genes = Genes.union(target_gene_dict.get(miRNA[:10] + '-5p', set()))
            Genes = Genes.union(target_gene_dict.get(miRNA + '-3p', set()))
            Genes = Genes.union(target_gene_dict.get(miRNA, set()))
        miRNA_Genes_dict[miRNA] = Genes

    for miRNA in Mir_list:
        Genes = set()
        # According to the naming rules of miRNA, convert the precursor miRNA into mature miRNA
        # The format of precursor miRNA is generally hsa-mir-xxx or hsa-mir-xxx-y
        # The format of mature miRNA is generally hsa-miR-xxx-5p or hsa-miR-xxx-3p or hsa-miR-xxx-y-5p or hsa-miR-xxx-y-3p
        # Therefore, we need to change mir to miR, and add -5p or -3p at the end
        # For example, hsa-mir-122 corresponds to mature miRNAs hsa-miR-122-5p and hsa-miR-122-3p
        # For example, hsa-mir-199a-1 corresponds to mature miRNAs hsa-miR-199a-1-5p and hsa-miR-199a-1-3p
        mature_miRNA_5p = miRNA.replace("mir", "miR") + "-5p"  # Convert to 5' end arm mature miRNA
        mature_miRNA_3p = miRNA.replace("mir", "miR") + "-3p"  # Convert to 3' end arm mature miRNA
        mature_miRNA = miRNA.replace("mir", "miR")
        # Directly find the corresponding target gene set from the dictionary and take the union
        Genes = Genes.union(target_gene_dict.get(mature_miRNA_5p, set()))
        Genes = Genes.union(target_gene_dict.get(mature_miRNA_3p, set()))
        Genes = Genes.union(target_gene_dict.get(mature_miRNA, set()))
        miRNA_Genes_dict[miRNA] = Genes
        # print(miRNA)
        # print(Genes)

    print(len(miRNA_Genes_dict))
    for miRNA in miRNA_Genes_dict:
        Genes = set()
        if miRNA_Genes_dict[miRNA] == set():
            mature_miRNA = miRNA.replace("mir", "miR")
            mature_miRNA = mature_miRNA[:-2]
            Genes = Genes.union(target_gene_dict.get(mature_miRNA + '-5p', set()))
            Genes = Genes.union(target_gene_dict.get(mature_miRNA + '-3p', set()))
            Genes = Genes.union(target_gene_dict.get(mature_miRNA, set()))
            miRNA_Genes_dict[miRNA] = Genes
    return miRNA_Genes_dict
