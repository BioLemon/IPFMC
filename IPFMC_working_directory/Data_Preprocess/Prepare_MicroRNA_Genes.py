import pandas as pd  # Import the pandas library for data manipulation and analysis

# pd.set_option('display.max_rows', None)  # Uncomment to display all rows in DataFrame output
# Set basic parameters for cancer types
Cancer_type_list = ['ACC','BRCA','COAD','KIRC','KIRP','LIHC','LUAD','LUSC','THYM']  # List of cancer types
miRNA_dir = '../Datas/Omics/Raw_Omics/'  # Directory containing raw miRNA omics data
Target_dir = '../Datas/Pathways/hsa_MTI.xlsx'  # Path to the miRNA target data Excel file
Output_dir = '../Datas/Pathways/'  # Directory to store processed pathway index files

target_gene_data = pd.read_excel(Target_dir)  # Read the miRNA target data from the Excel file
print(target_gene_data)  # Print the target gene data for verification

All_Names = set()  # Initialize an empty set to store all unique miRNA names

# Iterate over each cancer type to process corresponding miRNA data
for Cancer_type in Cancer_type_list:
    Temp_miRNA = pd.read_csv(miRNA_dir + Cancer_type + "/" + Cancer_type + "_miRNA.csv", index_col=0)  # Read miRNA data for the current cancer type
    Names = set(list(Temp_miRNA.index))  # Extract unique miRNA names from the DataFrame index
    All_Names = All_Names.union(Names)  # Update the set of all miRNA names
    print(Names)  # Print the names extracted for the current cancer type

print(All_Names)  # Print the complete set of all unique miRNA names

Let_list = []  # Initialize a list to store 'let' miRNA names
Mir_list = []  # Initialize a list to store 'mir' miRNA names

LET_PREFIX = 'hsa-let'  # Prefix for 'let' miRNA names
MIR_PREFIX = 'hsa-mir'  # Prefix for 'mir' miRNA names
PREFIX_LEN = 7  # Length of the miRNA prefix

# Classify miRNA names into 'let' and 'mir' lists based on their prefixes
for miRNA in All_Names:
    if miRNA[:PREFIX_LEN] == LET_PREFIX:  # Check if the miRNA name starts with 'hsa-let'
        Let_list.append(miRNA)  # Add to 'let' list
    elif miRNA[:PREFIX_LEN] == MIR_PREFIX:  # Check if the miRNA name starts with 'hsa-mir'
        Mir_list.append(miRNA)  # Add to 'mir' list

target_gene_dict = {}  # Initialize a dictionary to map miRNA to their target genes

# Populate the target gene dictionary from the target gene data
for i in range(len(target_gene_data['miRNA'])):  # Iterate over each row in target gene data
    MiR = target_gene_data['miRNA'][i]  # Get the miRNA name from the current row
    Gene = target_gene_data['Target Gene'][i]  # Get the target gene name from the current row
    if MiR not in target_gene_dict:  # If the miRNA is not already in the dictionary
        target_gene_dict[MiR] = set()  # Initialize an empty set for the miRNA
    target_gene_dict[MiR].add(Gene)  # Add the target gene to the miRNA's set of target genes

print(len(target_gene_dict))  # Print the total number of unique miRNAs in the dictionary
print(target_gene_dict['hsa-let-7a-5p'])  # Print the target genes for a specific 'let' miRNA
print(len(target_gene_dict['hsa-let-7a-5p']))  # Print the number of target genes for that specific miRNA

miRNA_Genes_dict = {}  # Initialize a dictionary to store genes associated with 'let' miRNAs

# Populate the miRNA-Gene dictionary for 'let' miRNAs
for miRNA in Let_list:
    Genes = set()  # Initialize an empty set for genes associated with the current miRNA
    if len(miRNA) == 10:  # Check if the miRNA name has a length of 10
        Genes = Genes.union(target_gene_dict.get(miRNA + '-5p', set()))  # Add genes for 5p arm
        Genes = Genes.union(target_gene_dict.get(miRNA + '-3p', set()))  # Add genes for 3p arm
        Genes = Genes.union(target_gene_dict.get(miRNA, set()))  # Add genes for the miRNA itself
    elif miRNA[10:] == '-1':  # Handle special case for miRNA names ending with '-1'
        Genes = Genes.union(target_gene_dict.get(miRNA[:10] + '-5p', set()))  # Add genes for 5p arm
        Genes = Genes.union(target_gene_dict.get(miRNA[:10] + '-3p', set()))  # Add genes for 3p arm
        Genes = Genes.union(target_gene_dict.get(miRNA, set()))  # Add genes for the miRNA itself
    else:  # Handle other cases for miRNA names
        Genes = Genes.union(target_gene_dict.get(miRNA[:10] + '-5p', set()))  # Add genes for 5p arm
        Genes = Genes.union(target_gene_dict.get(miRNA + '-3p', set()))  # Add genes for 3p arm
        Genes = Genes.union(target_gene_dict.get(miRNA, set()))  # Add genes for the miRNA itself
    miRNA_Genes_dict[miRNA] = Genes  # Store the set of genes for the current miRNA

# Populate the miRNA-Gene dictionary for 'mir' miRNAs
for miRNA in Mir_list:
    Genes = set()  # Initialize an empty set for genes associated with the current miRNA
    mature_miRNA_5p = miRNA.replace("mir", "miR") + "-5p"  # Create the mature 5p miRNA name
    mature_miRNA_3p = miRNA.replace("mir", "miR") + "-3p"  # Create the mature 3p miRNA name
    mature_miRNA = miRNA.replace("mir", "miR")  # Create the mature miRNA name
    Genes = Genes.union(target_gene_dict.get(mature_miRNA_5p, set()))  # Add genes for 5p arm
    Genes = Genes.union(target_gene_dict.get(mature_miRNA_3p, set()))  # Add genes for 3p arm
    Genes = Genes.union(target_gene_dict.get(mature_miRNA, set()))  # Add genes for the miRNA itself
    miRNA_Genes_dict[miRNA] = Genes  # Store the set of genes for the current miRNA
    #print(miRNA)  # Uncomment to print the current miRNA
    #print(Genes)  # Uncomment to print the associated genes

print(len(miRNA_Genes_dict))  # Print the number of unique miRNAs with associated genes
for miRNA in miRNA_Genes_dict:  # Iterate over each miRNA in the dictionary
    Genes = set()  # Initialize an empty set for genes
    if miRNA_Genes_dict[miRNA] == set():  # Check if the current miRNA has no associated genes
        mature_miRNA = miRNA.replace("mir", "miR")  # Convert 'mir' to 'miR' in the miRNA name
        mature_miRNA = mature_miRNA[:-2]  # Remove the last two characters from the miRNA name
        #print(mature_miRNA)  # Uncomment to print the mature miRNA name
        Genes = Genes.union(target_gene_dict.get(mature_miRNA + '-5p', set()))  # Add genes for 5p arm
        Genes = Genes.union(target_gene_dict.get(mature_miRNA + '-3p', set()))  # Add genes for 3p arm
        Genes = Genes.union(target_gene_dict.get(mature_miRNA, set()))  # Add genes for the miRNA itself
        #print(Genes)  # Uncomment to print the associated genes
        miRNA_Genes_dict[miRNA] = Genes  # Update the miRNA-Gene dictionary with the found genes
#print(miRNA_Genes_dict)  # Uncomment to print the final miRNA-Gene dictionary
