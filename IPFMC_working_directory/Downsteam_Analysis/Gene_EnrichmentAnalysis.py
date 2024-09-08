import warnings
import pandas as pd

# Suppress warnings for cleaner output
warnings.filterwarnings(action='ignore')

# List of cancer types to analyze
Cancer_type_list = ['ACC', 'BRCA', 'COAD', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'THYM']
# List of data types (omics) to analyze
Data_type_list = ['mRNA', 'miRNA', 'Methy', 'CNV']

# Specify the cancer type and data type for analysis
Cancer_type = 'LUAD'
Data_type = 'miRNA'

# Set the number of top genes to retrieve
num = 100

# Define directories for data
Omic_dir = '../Datas/Omics/Raw_Omics'  # Directory for raw omics data
Full_Path_dir = f'../Datas/Pathways/miRNA_Pathway_Index.csv'  # Pathway index file for miRNA
Sel_Path_dir = f'../Datas/Omics/Output_Omics'  # Directory for selected pathways output
Surv_dir = '../Datas/Survival/Subset_Survival'  # Directory for survival data

# Load omics data for the specified cancer type and data type
Omic_data = pd.read_csv(f'{Omic_dir}/{Cancer_type}/{Cancer_type}_{Data_type}.csv', index_col=0)
# Load survival data for the specified cancer type
Surv_data = pd.read_csv(f'{Surv_dir}/{Cancer_type}_Sub_survival.csv')
# Load full pathways data
Full_Pathways = pd.read_csv(Full_Path_dir, index_col=0)
# Load selected pathways data for the specified cancer type and data type
Sel_Pathways = pd.read_excel(f'{Sel_Path_dir}/{Cancer_type}_{Data_type}_Pathway_Sorting.xlsx', index_col=0)

# Create a list of all genes from the omics data
Full_Genes = list(Omic_data.index)
# Create a list of selected pathways
Sel_Pathways = list(Sel_Pathways['Pathway'])
# Extract the pathways corresponding to the selected pathways
Parcial_Pathways = Full_Pathways.loc[Sel_Pathways]

# Print the selected pathways and their corresponding genes
print(Parcial_Pathways)

# Initialize a dictionary to count occurrences of each gene in selected pathways
count_dict = {}

# Initialize the count for each gene to zero
for item in Full_Genes:
    count_dict[item] = 0

# Count the occurrences of each gene in the selected pathways
for Pathway in Sel_Pathways:
    Genes = Parcial_Pathways.loc[Pathway]  # Get genes associated with the current pathway
    Genes = list(Genes.dropna())  # Remove any NaN values from the gene list
    for i in Genes:
        if i in Full_Genes:  # Check if the gene is in the full gene list
            count_dict[i] += 1  # Increment the count for the gene

# Convert the count dictionary to a DataFrame for easier manipulation
Gene_Count = pd.DataFrame.from_dict(count_dict, orient='index', columns=['count'])
# Sort the DataFrame by count in descending order
Gene_Count = Gene_Count.sort_values(by='count', ascending=False)
# Select the top 'num' genes based on count
Gene_Count = Gene_Count.iloc[:num]

# Print the resulting gene counts
print(Gene_Count)

# Save the gene counts to an Excel file
Gene_Count.to_excel(f"../Datas/Assessment_Result/Gene_counts/{Data_type}_retopgenes.xlsx")