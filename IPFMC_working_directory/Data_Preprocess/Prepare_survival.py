import pandas as pd
import numpy as np

Cancer_type_list = ['BRCA', 'COAD', 'KIRC', 'LUAD', 'LUSC', 'ACC', 'KIRP', 'LIHC', 'THYM', 'BRCA_GOLD', 'COAD_GOLD']
input_dir = '../Datas/Survival/Raw_Survival/'
Omic_dir = '../Datas/Omics/Raw_Omics/'
Output_dir = '../Datas/Survival/Subset_Survival/'
for Cancer_type in Cancer_type_list:
    BSurv_data = pd.read_csv(input_dir + Cancer_type[:4] + '_survival.csv', index_col=0)
    Omic_data = pd.read_csv(Omic_dir + Cancer_type + '/' + Cancer_type[:4] + '_miRNA.csv', index_col=0)
    Omic_data = Omic_data.transpose()
    IDs = Omic_data.index
    print(IDs)
    IDs = IDs.to_numpy()
    IDs = IDs.astype(str)
    IDs = np.char.replace(IDs, ".", "-")
    Surv_data = BSurv_data.loc[IDs]
    pd.DataFrame(Surv_data).to_csv(Output_dir + Cancer_type + '_Sub_survival.csv', header=['event', 'time'])
