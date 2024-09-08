import pandas as pd  # Import the pandas library for data manipulation
import numpy as np  # Import the numpy library for numerical operations

# List of cancer types to process
Cancer_type_list = ['BRCA', 'COAD', 'KIRC', 'LUAD', 'LUSC', 'ACC', 'KIRP', 'LIHC', 'THYM', 'BRCA_GOLD', 'COAD_GOLD']

# Define input and output directories for survival and omics data
input_dir = '../Datas/Survival/Raw_Survival/'
Omic_dir = '../Datas/Omics/Raw_Omics/'
Output_dir = '../Datas/Survival/Subset_Survival/'

# Iterate over each cancer type in the list
for Cancer_type in Cancer_type_list:
    # Read the survival data for the current cancer type, using the first column as the index
    BSurv_data = pd.read_csv(input_dir + Cancer_type[:4] + '_survival.csv', index_col=0)

    # Read the omics data for the current cancer type, using the first column as the index
    Omic_data = pd.read_csv(Omic_dir + Cancer_type + '/' + Cancer_type[:4] + '_miRNA.csv', index_col=0)

    # Transpose the omics data to switch rows and columns
    Omic_data = Omic_data.transpose()

    # Get the index (IDs) of the transposed omics data
    IDs = Omic_data.index
    print(IDs)  # Print the IDs for debugging purposes

    # Convert the index to a numpy array
    IDs = IDs.to_numpy()

    # Ensure the IDs are of type string
    IDs = IDs.astype(str)

    # Replace any '.' characters in the IDs with '-' characters
    IDs = np.char.replace(IDs, ".", "-")

    # Filter the survival data to include only the rows that match the modified IDs
    Surv_data = BSurv_data.loc[IDs]

    # Save the filtered survival data to a new CSV file with specified headers
    pd.DataFrame(Surv_data).to_csv(Output_dir + Cancer_type + '_Sub_survival.csv', header=['event', 'time'])