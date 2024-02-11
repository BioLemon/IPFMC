import pandas as pd
import numpy as np

Cancer_type_list = ['ACC', 'BRCA', 'COAD', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'THYM']
Clinical_List = ['gender', 'age_at_initial_pathologic_diagnosis', 'pathologic_M', 'pathologic_N', 'pathologic_T',
                 'pathologic_stage']
Clinical_dir = '../Datas/Clinical/PCB'
for Cancer_type in Cancer_type_list:
    Clinical_data = pd.read_csv(f'{Clinical_dir}/{Cancer_type}_Clinical', delimiter='\t')
    if Cancer_type == 'THYM':
        Exist_Columns = set(Clinical_data.columns) & set(Clinical_List)
        Exist_Columns = list(set(Exist_Columns) | {'masaoka_stage'})
        print(Exist_Columns)
        Exist_Columns.insert(0, 'sampleID')
        sub_Clinical = Clinical_data[Exist_Columns]
        sub_Clinical = sub_Clinical.dropna()
        print(sub_Clinical)
        sub_Clinical.to_csv(f'{Clinical_dir}/{Cancer_type}_Sub_Clinical.csv', index=0)
        continue
    elif Cancer_type == 'ACC':
        Exist_Columns = set(Clinical_data.columns) & set(Clinical_List)
        Exist_Columns = list(set(Exist_Columns) | {'clinical_M'})
        print(Exist_Columns)
        Exist_Columns.insert(0, 'sampleID')
        sub_Clinical = Clinical_data[Exist_Columns]
        sub_Clinical = sub_Clinical.dropna()
        print(sub_Clinical)
        sub_Clinical.to_csv(f'{Clinical_dir}/{Cancer_type}_Sub_Clinical.csv', index=0)
        continue
    print(Cancer_type)
    #print(Clinical_data)
    Exist_Columns = set(Clinical_data.columns) & set(Clinical_List)
    Exist_Columns = list(Exist_Columns)
    Exist_Columns.insert(0,'sampleID')
    sub_Clinical = Clinical_data[Exist_Columns]
    sub_Clinical = sub_Clinical.dropna()
    print(sub_Clinical)
    sub_Clinical.to_csv(f'{Clinical_dir}/{Cancer_type}_Sub_Clinical.csv',index=0)







