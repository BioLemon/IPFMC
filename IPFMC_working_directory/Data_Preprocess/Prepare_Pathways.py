import pandas as pd
BP_dir = '../Datas/Pathways/'
BioProcesses = pd.read_json(BP_dir+'GoTerm_NCP.json')
Process_lists = []
Pronames_lists = []
print(BioProcesses)

for i in BioProcesses:
    a = BioProcesses[i]['geneSymbols']
    if(len(a)>10 and len(a)<200):
        Process_lists.append(a)
        Pronames_lists.append(i)
Process_lists = pd.DataFrame(Process_lists)
Pronames_lists = pd.DataFrame(Pronames_lists)
Pronames_lists = Pronames_lists.rename(columns={0: 'Pathways'})
# print(Pronames_lists)
A = [Pronames_lists,Process_lists]
A = pd.concat(A,axis=1)

print(A)
A.set_index(A.columns[0],inplace=True)
print(A)

Output_dir = '../Datas/Pathways/'
pd.DataFrame(A).to_csv(Output_dir+'Pathway_Index.csv')
