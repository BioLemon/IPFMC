import pandas as pd  # Import the pandas library for data manipulation

# Define the directory containing biological process data
BP_dir = '../Datas/Pathways/'

# Read biological process data from a JSON file
BioProcesses = pd.read_json(BP_dir + 'GoTerm_NCP.json')

# Initialize two empty lists to store gene symbols and pathway names that meet the criteria
Process_lists = []
Pronames_lists = []

# Print the loaded biological process data
print(BioProcesses)

# Iterate over each entry in BioProcesses
for i in BioProcesses:
    # Retrieve the gene symbols for the current entry
    a = BioProcesses[i]['geneSymbols']

    # Check if the number of gene symbols is between 10 and 200
    if (len(a) > 10 and len(a) < 200):
        # If the condition is met, append the gene symbols and pathway name to the respective lists
        Process_lists.append(a)
        Pronames_lists.append(i)

# Convert the list of gene symbols that meet the criteria into a DataFrame
Process_lists = pd.DataFrame(Process_lists)

# Convert the list of pathway names into a DataFrame and rename the column to 'Pathways'
Pronames_lists = pd.DataFrame(Pronames_lists)
Pronames_lists = Pronames_lists.rename(columns={0: 'Pathways'})

# Merge the two DataFrames into one, concatenating them by columns
A = [Pronames_lists, Process_lists]
A = pd.concat(A, axis=1)

# Print the merged DataFrame
print(A)

# Set the first column of the merged DataFrame as the index
A.set_index(A.columns[0], inplace=True)

# Print the DataFrame after setting the index
print(A)

# Define the output directory
Output_dir = '../Datas/Pathways/'

# Save the final DataFrame as a CSV file
pd.DataFrame(A).to_csv(Output_dir + 'Pathway_Index.csv')