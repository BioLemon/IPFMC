import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import warnings

warnings.filterwarnings("ignore")

def Create_pathwayindex(GoTerminfo, min_genes, max_genes):
    """
    Create a pathway index from GSEA pathway information.

    Args:
        GoTerminfo (dict): Pathway information in JSON format downloaded from GSEA.
        min_genes (int): Minimum number of genes required for a pathway to be included in the index.
        max_genes (int): Maximum number of genes allowed for a pathway to be included in the index.

    Returns:
        pathway_index (DataFrame): DataFrame containing the pathway index.
    """
    process_lists = []  # List to store gene sets for each pathway
    pronames_lists = []  # List to store pathway names

    # Iterate through the pathway information
    for pathway_name, pathway_info in GoTerminfo.items():
        gene_symbols = pathway_info['geneSymbols']  # Get the gene symbols for the current pathway
        if min_genes < len(gene_symbols) < max_genes:  # Check if the number of genes is within the specified range
            process_lists.append(gene_symbols)  # Add the gene set to the process_lists
            pronames_lists.append(pathway_name)  # Add the pathway name to the pronames_lists

    # Create DataFrames from the lists
    process_lists = pd.DataFrame(process_lists)
    pronames_lists = pd.DataFrame(pronames_lists)
    pronames_lists = pronames_lists.rename(columns={0: 'Pathways'})

    # Concatenate the DataFrames
    pathway_index = pd.concat([pronames_lists, process_lists], axis=1)

    # Set the index of the DataFrame to the pathway names
    pathway_index.set_index(pathway_index.columns[0], inplace=True)

    return pathway_index  # Return the pathway index

def Create_miRNA_pathwayindex(miRNA_Genes_dict, pathwayinfo):
    """
    Create a miRNA-pathway index by finding significant associations between miRNAs and pathways.

    Args:
        miRNA_Genes_dict (dict): Dictionary mapping miRNAs to their target genes.
        pathwayinfo (DataFrame): DataFrame containing pathway information.

    Returns:
        miRNA_pathway_index (DataFrame): DataFrame containing the miRNA-pathway index.
    """
    # Read the pathway name and gene name columns from the pathway information
    bp_data = pathwayinfo
    pathways = list(bp_data.index)
    gene_symbols = []

    # Import the pathway gene sets into the gene_symbols list
    for pathway in pathways:
        pathway_genes = bp_data.loc[pathway]  # Get the genes associated with the current pathway
        pathway_genes = pathway_genes.dropna()  # Remove NaN values
        temp_list = []
        for gene in pathway_genes:
            temp_list.append(gene)
        gene_symbols.append(temp_list)

    valid_pathways = []  # List to store valid pathways (pathways included in the dataset)
    valid_mirna_symbols = []  # List to store the miRNAs associated with valid pathways

    # Convert the miRNA target gene list to a 2D NumPy array
    mirna_genes_array = np.array(list(miRNA_Genes_dict.values()))

    # Convert the miRNA name list to a 1D NumPy array
    mirna_names_array = np.array(list(miRNA_Genes_dict.keys()))

    # Get the total number of genes in the population
    n = 20000

    # Iterate through the pathways
    for i, pathway_genes in enumerate(gene_symbols):
        # Get the number of genes in the current pathway
        m = len(pathway_genes)

        # Calculate the number of intersections between each miRNA target gene set and the current pathway gene set
        intersection_count = np.array([len(set(pathway_genes) & set(target_genes)) for target_genes in mirna_genes_array])

        # Calculate the number of genes targeted by each miRNA
        target_count = np.array([len(target_genes) for target_genes in mirna_genes_array])

        # Calculate the p-values using the hypergeometric distribution test
        p_values = hypergeom.sf(intersection_count, n, m, target_count)

        # Determine if the p-values are significant (less than 0.005)
        is_significant = p_values < 0.005

        # If any p-value is significant, add the pathway to the valid_pathways list and the corresponding miRNAs to the valid_mirna_symbols list
        if is_significant.any():
            valid_pathways.append(pathways[i])
            valid_mirna_symbols.append(mirna_names_array[is_significant])

        # Print progress and the number of miRNAs associated with the current pathway
        print(f'Current progress: {i}/{len(gene_symbols)}, Number of miRNAs in the current pathway: {len(mirna_names_array[is_significant])}')

    # Create DataFrames from the valid_pathways and valid_mirna_symbols lists
    valid_pathways = pd.DataFrame(valid_pathways)
    valid_mirna_symbols = pd.DataFrame(valid_mirna_symbols)
    valid_pathways = valid_pathways.rename(columns={0: 'Pathways'})

    # Concatenate the DataFrames
    mirna_pathway_index = pd.concat([valid_pathways, valid_mirna_symbols], axis=1)

    # Set the index of the DataFrame to the pathway names
    mirna_pathway_index.set_index(mirna_pathway_index.columns[0], inplace=True)

    print('Successfully created miRNA-pathway index')
    return mirna_pathway_index  # Return the miRNA-pathway index

def gene_occurrence(Full_Genes, Sel_Pathways, Full_Pathways, top=50):
    """
    Calculate the occurrence count of genes in selected pathways.

    Args:
        Full_Genes (list): List of all genes.
        Sel_Pathways (DataFrame): DataFrame containing selected pathways.
        Full_Pathways (DataFrame): DataFrame containing full pathway information.
        top (int, optional): Number of top genes to return. Defaults to 50.

    Returns:
        Gene_Count (DataFrame): DataFrame containing the top genes and their occurrence counts.
    """
    Sel_Pathways = list(Sel_Pathways['Pathway'])  # Get all the pathway names after screening, and overwrite the original Sel_Pathways variable
    Parcial_Pathways = Full_Pathways.loc[Sel_Pathways]  # Get partial pathways based on the selected pathways

    # Define a dictionary to store the count of each gene
    count_dict = {}

    # Initialize the count of each gene to 0
    for gene in Full_Genes:
        count_dict[gene] = 0

    # Count the occurrence of each gene in the selected pathways
    for pathway in Sel_Pathways:
        genes = Parcial_Pathways.loc[pathway]  # Get genes associated with the current pathway
        genes = list(genes.dropna())  # Remove NaN values
        for gene in genes:
            if gene in Full_Genes:
                count_dict[gene] += 1  # Increment the count for the gene

    # Create a DataFrame from the count dictionary
    Gene_Count = pd.DataFrame.from_dict(count_dict, orient='index', columns=['count'])

    # Sort the DataFrame by the count column in descending order
    Gene_Count = Gene_Count.sort_values(by='count', ascending=False)

    # Select the top genes based on the specified number
    Gene_Count = Gene_Count.iloc[:top]

    return Gene_Count  # Return the DataFrame containing the top genes and their counts

def create_mirTar_Dict(MTI_info, mir_set):
    """
    Create a dictionary mapping miRNAs to their target genes.

    Args:
        MTI_info (DataFrame): DataFrame containing miRNA-target gene information.
        mir_set (list): List of miRNAs.

    Returns:
        miRNA_Genes_dict (dict): Dictionary mapping miRNAs to their target genes.
    """
    let_list = []  # List to store let-7 miRNAs
    mir_list = []  # List to store miR miRNAs
    all_names = mir_set  # Alias for mir_set
    target_gene_data = MTI_info  # Alias for MTI_info
    let_prefix = 'hsa-let'  # Prefix for let-7 miRNAs
    mir_prefix = 'hsa-mir'  # Prefix for miR miRNAs
    prefix_len = 7  # Length of the prefix

    # Separate let-7 and miR miRNAs into their respective lists
    for miRNA in all_names:
        if miRNA[:prefix_len] == let_prefix:
            let_list.append(miRNA)
        elif miRNA[:prefix_len] == mir_prefix:
            mir_list.append(miRNA)

    # Convert the target gene data into a dictionary
    target_gene_dict = {}
    for i in range(len(target_gene_data['miRNA'])):
        miR = target_gene_data['miRNA'][i]
        gene = target_gene_data['Target Gene'][i]
        if miR not in target_gene_dict:
            target_gene_dict[miR] = set()
        target_gene_dict[miR].add(gene)

    # Create a dictionary mapping miRNAs to their target genes
    miRNA_Genes_dict = {}
    for miRNA in let_list:
        genes = set()
        if len(miRNA) == 10:  # For miRNAs similar to hsa-let-7c
            genes = genes.union(target_gene_dict.get(miRNA + '-5p', set()))
            genes = genes.union(target_gene_dict.get(miRNA + '-3p', set()))
            genes = genes.union(target_gene_dict.get(miRNA, set()))
        elif miRNA[10:] == '-1':  # For miRNAs similar to hsa-let-7a-1
            genes = genes.union(target_gene_dict.get(miRNA[:10] + '-5p', set()))
            genes = genes.union(target_gene_dict.get(miRNA[:10] + '-3p', set()))
            genes = genes.union(target_gene_dict.get(miRNA, set()))
        else:  # For other miRNAs, such as hsa-let-7a-2
            genes = genes.union(target_gene_dict.get(miRNA[:10] + '-5p', set()))
            genes = genes.union(target_gene_dict.get(miRNA + '-3p', set()))
            genes = genes.union(target_gene_dict.get(miRNA, set()))
        miRNA_Genes_dict[miRNA] = genes

    return miRNA_Genes_dict  # Return the dictionary mapping miRNAs to their target genes