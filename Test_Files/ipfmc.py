import argparse
import pandas as pd
import numpy as np
from snf import snf
from IPFMC import separate, direct, analysis
import sys
import pickle
import pandas as pd
import random
def read_omics(omics_file):
    """Read omics data file and return a dictionary of omics types and corresponding data paths."""
    with open(omics_file, 'r') as f:
        lines = f.readlines()

    # Read omics types from the first line
    omics_types = lines[0].strip().split()
    omics_paths = [line.strip() for line in lines[1:]]

    # Ensure that the number of omics types matches the number of paths
    if len(omics_types) != len(omics_paths):
        raise ValueError("The number of omics types and paths do not match.")

    # Check if there are fewer than 2 omics types
    if len(omics_types) < 2:
        raise ValueError("Fewer than 2 omics types provided, at least two are required.")

    # Create a dictionary to store omics data
    omics_data = {}
    for omics_type, path in zip(omics_types, omics_paths):
        omics_data[omics_type] = pd.read_csv(path, index_col=0)  # Assuming CSV file format

    return omics_data

def read_pathways(pathway_file):
    """Read pathway information file and return pathway_index and miRNA_pathway_index."""
    with open(pathway_file, 'r') as f:
        lines = f.readlines()

    # Read the first two lines for paths
    if len(lines) < 2:
        raise ValueError("The pathway file must contain at least two lines.")

    pathway_index_path = lines[0].strip()
    mirna_pathway_index_path = lines[1].strip()

    # Read the corresponding CSV files
    pathway_index = pd.read_csv(pathway_index_path, index_col=0)
    mirna_pathway_index = pd.read_csv(mirna_pathway_index_path, index_col=0)

    return pathway_index, mirna_pathway_index


def cluster_st1(omics_data, pathway_index, mirna_pathway_index, output_dir, cancer_type, num_cluster=None,
                pathway_prop=100, seed=None):
    # Check if miRNA omics data exists in the provided omics_data dictionary
    if 'miRNA' in omics_data:
        mirna_data = omics_data.pop('miRNA')  # Extract miRNA data and remove it from the omics_data dictionary
        print("Extracting miRNA omics data, processing other omics data")
    else:
        mirna_data = None  # Set miRNA data to None if it doesn't exist
        print("No miRNA omics data found")

    # Prepare to process the remaining omics data
    omic_list = [data for data in omics_data.values()]  # Create a list of remaining omics data
    represents = []  # List to hold representations of each omics type
    pathways_list = {}  # Dictionary to hold pathways associated with each omics type

    # Process each type of omics data
    for omic_type, omic in omics_data.items():
        # Use the ipfmc_discretize method to discretize the omics data based on the pathway index
        represent, pathways = separate.ipfmc_discretize(omic, pathway_index, preselect=int(pathway_prop),
                                                        seed=seed)  # Apply Strategy 1
        represents.append(np.array(represent))  # Convert the representation to a NumPy array and add to the list
        pathways_list[omic_type] = pathways  # Store the pathways for the current omics type
        pathways_df = pd.DataFrame(pathways)  # Convert pathways to a DataFrame for saving
        # Save the pathways associated with the current omics type to a CSV file
        pathways_df.to_csv(f"{output_dir}/{cancer_type}_{omic_type}_pathways.csv", index=False)

    # Process miRNA data if it exists
    if mirna_data is not None:
        # Discretize the miRNA data using the miRNA pathway index
        represent, pathways = separate.ipfmc_discretize(mirna_data, mirna_pathway_index,
                                                        preselect=int(pathway_prop))  # Use Strategy 2 for miRNA data
        represents.append(np.array(represent))  # Convert the representation to a NumPy array and add to the list
        pathways_list['miRNA'] = pathways  # Store the pathways for miRNA
        pathways_df = pd.DataFrame(pathways)  # Convert pathways to a DataFrame for saving
        # Save the pathways associated with miRNA to a CSV file
        pathways_df.to_csv(f"{output_dir}/{cancer_type}_miRNA_pathways.csv", index=False)

    # Final aggregated representation of all omics data
    represent_final = snf(represents,
                          K=15)  # Perform spectral clustering on the representations, K value can be adjusted as needed

    # Determine the number of clusters to use for clustering
    if num_cluster is None:
        K = separate.suggest_k(
            represent_final)  # Get the suggested number of clusters based on the final representation
    else:
        K = num_cluster  # Use the user-specified number of clusters if provided

    # Perform spectral clustering to obtain labels for each sample in the first omics data type
    labels = separate.spec_cluster(omic_list[0], fusion_matrix=represent_final,
                                   k=K)  # Clustering labels based on the final representation
    labels.columns = ["cluster"]  # Rename the column to "cluster"

    return labels, K, pathways_list  # Return the clustering labels, number of clusters, and pathways list


def cluster_st2(omics_data, pathway_index, mirna_pathway_index, output_dir, cancer_type, num_cluster=None):
    # Check if miRNA omics data exists in the provided omics_data dictionary
    if 'miRNA' in omics_data:
        mirna_data = omics_data.pop('miRNA')  # Extract miRNA data and remove it from the omics_data dictionary
        print("Extracting miRNA omics data, processing other omics data")
    else:
        mirna_data = None  # Set miRNA data to None if it doesn't exist
        print("No miRNA omics data found")

    # Prepare to process the remaining omics data
    omic_list = [data for data in omics_data.values()]  # Create a list of remaining omics data
    represents = []  # List to hold representations of each omics type
    pathways_list = {}  # Dictionary to hold pathways associated with each omics type

    # Process each type of omics data
    for omic_type, omic in omics_data.items():
        # Use the ipfmc_average method to average the omics data based on the pathway index
        represent, pathways = separate.ipfmc_average(omic, pathway_index)  # Apply Strategy 2
        represents.append(np.array(represent))  # Convert the representation to a NumPy array and add to the list
        pathways_list[omic_type] = pathways  # Store the pathways for the current omics type
        pathways_df = pd.DataFrame(pathways)  # Convert pathways to a DataFrame for saving
        # Save the pathways associated with the current omics type to a CSV file
        pathways_df.to_csv(f"{output_dir}/{cancer_type}_{omic_type}_pathways.csv", index=False)

    # Process miRNA data if it exists
    if mirna_data is not None:
        # Average the miRNA data using the miRNA pathway index
        represent, pathways = separate.ipfmc_average(mirna_data, mirna_pathway_index)  # Use Strategy 2 for miRNA data
        represents.append(np.array(represent))  # Convert the representation to a NumPy array and add to the list
        pathways_list['miRNA'] = pathways  # Store the pathways for miRNA
        pathways_df = pd.DataFrame(pathways)  # Convert pathways to a DataFrame for saving
        # Save the pathways associated with miRNA to a CSV file
        pathways_df.to_csv(f"{output_dir}/{cancer_type}_miRNA_pathways.csv", index=False)

    # Final aggregated representation of all omics data
    represent_final = snf(represents,
                          K=15)  # Perform spectral clustering on the representations, K value can be adjusted as needed

    # Determine the number of clusters to use for clustering
    if num_cluster is None:
        K = separate.suggest_k(
            represent_final)  # Get the suggested number of clusters based on the final representation
    else:
        K = num_cluster  # Use the user-specified number of clusters if provided

    # Perform spectral clustering to obtain labels for each sample in the first omics data type
    labels = separate.spec_cluster(omic_list[0], fusion_matrix=represent_final,
                                   k=K)  # Clustering labels based on the final representation
    labels.columns = ["cluster"]  # Rename the column to "cluster"

    return labels, K, pathways_list  # Return the clustering labels, number of clusters, and pathways list
def create_pathway_index(BP_dir):
    """
    Create a pathway index from a JSON file containing biological pathway information.

    Args:
        BP_dir (str): Path to the JSON file containing biological pathway information.

    Returns:
        pathway_index (pandas.DataFrame): A DataFrame containing the pathway index.
    """
    Bio_Pathways = pd.read_json(BP_dir)
    pathway_index = analysis.Create_pathwayindex(Bio_Pathways, min_genes=10, max_genes=200)
    print("Completed constructing pathway index.")
    return pathway_index

def create_miRTar_index(mirTar_dir, mir_set_dir, pathway_dir):
    """
    Create a miRNA-pathway index from miRNA target information, miRNA set, and biological pathway information.

    Args:
        mirTar_dir (str): Path to the Excel file containing miRNA target information.
        mir_set_dir (str): Path to the CSV file containing the miRNA set.
        pathway_dir (str): Path to the CSV file containing biological pathway information.

    Returns:
        mirPathway_index (dict): A dictionary containing the miRNA-pathway index.
    """
    mirTar_info = pd.read_excel(mirTar_dir)
    mirset = pd.read_csv(mir_set_dir, header=None)
    Bio_Pathways = pd.read_csv(pathway_dir, index_col=0)
    print(Bio_Pathways)
    print(mirset)
    mirset = mirset[0].to_list()
    mirTar_index = analysis.create_mirTar_Dict(mirTar_info, mirset)
    mirPathway_index = analysis.Create_miRNA_pathwayindex(mirTar_index, Bio_Pathways)
    return mirPathway_index

def main():
    seed = 10
    random.seed(seed)
    parser = argparse.ArgumentParser(description="IPFMC")

    # Add simplified arguments
    parser.add_argument('-p', '--pathway', required=True, help='Path to the pathway information file')
    parser.add_argument('-m', '--mirset', default=None, help='Path to the miRNA set file')
    parser.add_argument('-o', '--omics', default=None, help='Path to the omics information file')
    parser.add_argument('-a', '--action', type=str, choices=['c', 'p', 'm'], required=True, help='Action to perform: c, p, m')
    parser.add_argument('-s', '--strategy', type=int, default=1, help='Strategy to use: 1 or 2')
    parser.add_argument('-n', '--num_cluster', type=int, help='Number of clusters, default is the suggested number')
    parser.add_argument('-d', '--output', type=str, default='./results', help='Output directory, default is results folder in the current directory')
    parser.add_argument('-co', '--count', type=int, default=0, help='Output directory, default is results folder in the current directory')
    parser.add_argument('-can', '--cancer', type=str, default="Cancer", help="Cancer type of input datasets")
    parser.add_argument('-se', '--select', default=100, help='Proportion of retained pathways')
    args = parser.parse_args()

    # Read the arguments
    pathway_file = args.pathway
    omics_file = args.omics
    action = args.action
    strategy = args.strategy
    num_cluster = args.num_cluster
    output_dir = args.output
    compute_count = args.count
    cancer_type = args.cancer
    mir_set_file = args.mirset
    pathway_proportion = args.select

    # Check if the strategy is valid
    if strategy != 1 and strategy != 2:
        print("Please set strategy to 1 or 2!")  # Prompt the user to set a valid strategy
        sys.exit(1)  # Terminate the program and return status code 1 to indicate an error

    # Execute different operations based on the action specified
    if action == "c":
        # Read omics data from the specified file
        omics_data = read_omics(omics_file)  # Load omics data into a DataFrame
        origin_omicsdata = omics_data.copy()  # Create a copy of the original omics data for later use

        # Read pathway data from the specified file
        pathway_index, mirna_pathway_index = read_pathways(pathway_file)  # Load pathway data into DataFrames

        # Perform clustering based on the selected strategy
        if strategy == 1:
            # Execute clustering using strategy 1
            Labels, K, Pathways_list = cluster_st1(omics_data, pathway_index, mirna_pathway_index, output_dir,
                                                   cancer_type, num_cluster, pathway_proportion, seed=seed)
            # Save the clustering results to a CSV file
            Labels.to_csv(f"{output_dir}/{cancer_type}_{K}.csv")
        elif strategy == 2:
            # Execute clustering using strategy 2
            Labels, K, Pathways_list = cluster_st2(omics_data, pathway_index, mirna_pathway_index, output_dir,
                                                   cancer_type, num_cluster)
            # Save the clustering results to a CSV file
            Labels.to_csv(f"{output_dir}/{cancer_type}_{K}.csv")

        # If compute_count is set to 1, calculate gene occurrence in pathways
        if compute_count == 1:
            for omic_type, omic in origin_omicsdata.items():
                # Process each omics type except miRNA
                if omic_type != "miRNA":
                    Full_Genes = list(omic.index)  # Get the list of genes from the omics data
                    # Calculate gene occurrence in pathways for the current omics type
                    Gene_counts = analysis.gene_occurrence(Full_Genes, Pathways_list[omic_type], pathway_index)
                    # Save the gene occurrence counts to a CSV file
                    Gene_counts.to_csv(f"{output_dir}/{cancer_type}_{omic_type}_counts.csv")
                else:
                    # Process miRNA separately
                    Full_Genes = list(omic.index)  # Get the list of genes from the miRNA data
                    # Calculate gene occurrence in pathways for miRNA
                    Gene_counts = analysis.gene_occurrence(Full_Genes, Pathways_list[omic_type], mirna_pathway_index)
                    # Save the gene occurrence counts to a CSV file
                    Gene_counts.to_csv(f"{output_dir}/{cancer_type}_{omic_type}_counts.csv")

    elif action == "p":
        # If action is "p", create a pathway index
        pathway_index = create_pathway_index(pathway_file)  # Create the pathway index from the specified file
        # Save the pathway index to a CSV file
        pathway_index.to_csv(f"{output_dir}/Pathway_index.csv")

    elif action == "m":
        # If action is "m", create a miRNA-pathway index
        mirna_index = create_miRTar_index(omics_file, mir_set_file, pathway_file)  # Create the miRNA-pathway index
        # Save the miRNA-pathway index to a CSV file
        mirna_index.to_csv(f"{output_dir}/miRNA_Pathway_index.csv")

    else:
        # If the action is not recognized, prompt the user for a correct operation
        print("Please provide a valid operation.")

# Entry point of the script
if __name__ == "__main__":
        main()  # Call the main function to execute the program





