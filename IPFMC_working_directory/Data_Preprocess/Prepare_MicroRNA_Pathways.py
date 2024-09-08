import pandas as pd  # Import the pandas library for data manipulation
import numpy as np  # Import the numpy library for numerical calculations
from scipy.stats import hypergeom  # Import the hypergeom function from scipy.stats for hypergeometric distribution
import time  # Import the time module for timing the execution
import pickle  # Import the pickle module for serializing and deserializing objects

start_time = time.time()  # Record the start time of the program
BP_dir = '../Datas/Pathways/Pathway_Index.csv'  # Path to the biological process data
BP_data = pd.read_csv(BP_dir, index_col=0)  # Read the biological process data from a CSV file
Pathways = list(BP_data.index)  # Get a list of all pathway names
GeneSymbols = []  # Initialize an empty list to store gene symbols

# Iterate over each pathway to extract gene symbols
for i in Pathways:
    A = BP_data.loc[i]  # Get the genes corresponding to the current pathway
    A = A.dropna()  # Remove any missing values
    Temp_list = []  # Initialize a temporary list to hold gene symbols
    for value in A:  # Iterate over the genes in the current pathway
        Temp_list.append(value)  # Append each gene to the temporary list
    GeneSymbols.append(Temp_list)  # Add the temporary list to the GeneSymbols list

# Load the miRNA-gene dictionary from a pickle file
f = open("../Datas/Pathways/miRNA_Genes_dict.pkl", "rb")  # Open the miRNA-gene dictionary file
miRNA_Genes_dict = pickle.load(f)  # Load the miRNA-gene dictionary from the file
f.close()  # Close the file
print(len(miRNA_Genes_dict))  # Print the length of the miRNA-gene dictionary

Valid_Pathways = []  # Initialize an empty list to store valid pathways
Valid_miRNASymbols = []  # Initialize an empty list to store valid miRNA symbols

# Convert the miRNA-gene dictionary values to a numpy array
miRNA_Genes_array = np.array(list(miRNA_Genes_dict.values()))  # Create an array of miRNA gene sets
print(len(miRNA_Genes_array))  # Print the length of the miRNA-gene array
print(miRNA_Genes_array)  # Print the miRNA-gene array
miRNA_Names_array = np.array(list(miRNA_Genes_dict.keys()))  # Create an array of miRNA names

N = 20000  # Set the total number of genes

# Iterate over each pathway to compute significance
for i in range(len(GeneSymbols)):
    Pathway_Genes = GeneSymbols[i]  # Get the genes for the current pathway
    M = len(Pathway_Genes)  # Get the number of genes in the current pathway
    # Calculate the intersection count between pathway genes and each miRNA target gene set
    intersection_count = np.array([len(set(Pathway_Genes) & set(Target_Genes)) for Target_Genes in miRNA_Genes_array])
    target_count = np.array(
        [len(Target_Genes) for Target_Genes in miRNA_Genes_array])  # Get the number of target genes for each miRNA
    p_values = hypergeom.sf(intersection_count, N, M,
                            target_count)  # Calculate the p-values using the hypergeometric distribution
    is_significant = p_values < 0.005  # Determine if p-values are less than 0.005 (significant)

    if is_significant.any():  # If there are any significant miRNAs
        Valid_Pathways.append(Pathways[i])  # Add the current pathway to the valid pathways list
        Valid_miRNASymbols.append(
            miRNA_Names_array[is_significant])  # Add significant miRNAs to the valid miRNA symbols list

    # Print the progress of the computation
    print(
        f'Computing: {i}/{len(GeneSymbols)}, miRNA length of current pathway: {len(miRNA_Names_array[is_significant])}')

print(len(Valid_Pathways))  # Print the number of valid pathways
print(len(Valid_miRNASymbols))  # Print the number of valid miRNAs
Valid_Pathways = pd.DataFrame(Valid_Pathways)  # Convert the valid pathways list to a DataFrame
Valid_miRNASymbols = pd.DataFrame(Valid_miRNASymbols)  # Convert the valid miRNAs list to a DataFrame
Valid_Pathways = Valid_Pathways.rename(columns={0: 'Pathways'})  # Rename the DataFrame column to 'Pathways'
A = [Valid_Pathways, Valid_miRNASymbols]  # Combine the two DataFrames into a list
A = pd.concat(A, axis=1)  # Concatenate the DataFrames along the columns
print(A)  # Print the combined DataFrame
A.set_index(A.columns[0], inplace=True)  # Set the first column as the index of the DataFrame

end_time = time.time()  # Record the end time of the program
print(f'The total running time of the program is {end_time - start_time} seconds')  # Print the total running time
Output_dir = '../Datas/Pathways/'  # Set the output directory
pd.DataFrame(A).to_csv(Output_dir + 'miRNA_Pathway_Index.csv')  # Save the DataFrame to a CSV file




