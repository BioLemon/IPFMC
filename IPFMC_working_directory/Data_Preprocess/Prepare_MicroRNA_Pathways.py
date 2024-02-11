import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import time
import pickle

start_time = time.time()
BP_dir = '../Datas/Pathways/Pathway_Index.csv'
BP_data = pd.read_csv(BP_dir,index_col=0)
Pathways = list(BP_data.index)
GeneSymbols = []
for i in Pathways:
    A = BP_data.loc[i]
    A = A.dropna()
    Temp_list = []
    for value in A:
        Temp_list.append(value)
    GeneSymbols.append(Temp_list)


f = open("../Datas/Pathways/miRNA_Genes_dict.pkl", "rb")
miRNA_Genes_dict = pickle.load(f)
f.close()
print(len(miRNA_Genes_dict))

Valid_Pathways = []
Valid_miRNASymbols = []

miRNA_Genes_array = np.array(list(miRNA_Genes_dict.values()))
print(len(miRNA_Genes_array))
print(miRNA_Genes_array)
miRNA_Names_array = np.array(list(miRNA_Genes_dict.keys()))

N = 20000

for i in range(len(GeneSymbols)):
    Pathway_Genes = GeneSymbols[i]
    M = len(Pathway_Genes)
    intersection_count = np.array([len(set(Pathway_Genes) & set(Target_Genes)) for Target_Genes in miRNA_Genes_array])
    target_count = np.array([len(Target_Genes) for Target_Genes in miRNA_Genes_array])
    p_values = hypergeom.sf(intersection_count, N, M, target_count)
    is_significant = p_values < 0.005
    if is_significant.any():
        Valid_Pathways.append(Pathways[i])
        Valid_miRNASymbols.append(miRNA_Names_array[is_significant])
    print(f'Computing:{i}/{len(GeneSymbols)}，miRNA length of current pathway：{len(miRNA_Names_array[is_significant])}')

print(len(Valid_Pathways))
print(len(Valid_miRNASymbols))
Valid_Pathways= pd.DataFrame(Valid_Pathways)
Valid_miRNASymbols = pd.DataFrame(Valid_miRNASymbols)
Valid_Pathways = Valid_Pathways.rename(columns={0: 'Pathways'})
A = [Valid_Pathways, Valid_miRNASymbols]
A = pd.concat(A,axis=1)
print(A)
A.set_index(A.columns[0],inplace=True)
print(A)

end_time = time.time()
print(f'The total running time of the program is {end_time-start_time} seconds')
Output_dir = '../Datas/Pathways/'
pd.DataFrame(A).to_csv(Output_dir+'miRNA_Pathway_Index.csv')





