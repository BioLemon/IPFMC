# IPFMC

## Usage of IPFMC

In this section, we will introduce how to use our python package published on PyPI to perform clustering and biological interpretation of cancer multi-omics data. **We will show the process using the LUAD cancer datasets as an example.**

### Step 1: Prepare datasets
Before you start, you can download our repository by clicking the green ‘code’ button in the upper right corner and selecting “Download ZIP” from the drop-down menu. We will use the data in the “Test_Files” folder.

1. **Omics datasets**

   We apologize for the inconvenience caused by the large size of the cancer omics data file, which prevents us from uploading it to the github repository. Therefore, we provide the source of the omics data used in our experiments and show the downloading method using the omics data of LUAD cancer as an example. Please follow our guidance to download the LUAD cancer datasets.

   Our omics data is derived from the comprehensive study conducted by Duan et al., as published in *PLOS Computational Biology* in 2021. Their research evaluated and compared multi-omics data integration methods for cancer subtyping. You can find their original paper here: ([Evaluation and comparison of multi-omics data integration methods for cancer subtyping | PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009224)). 

   Their datasets can be downloaded in their github repository: [GaoLabXDU/MultiOmicsIntegrationStudy: Experimental evaluation and comparison of multi-omics data integration methods for cancer subtyping (github.com)](https://github.com/GaoLabXDU/MultiOmicsIntegrationStudy/).

   After accessing their github repository, go to the ‘Availability of datasets’ section in their README file and locate the ‘Dataset #1 LUAD Complete’ dataset. Choose any of the download links they offer and download it. You will get a rar file with these four files inside:

   (1) LUAD_mRNA.csv
   (2) LUAD_miRNA.csv
   (3) LUAD_Methy.csv
   (4) LUAD_CNV.csv

   You can also download their other Complete datasets to perform experiments on other cancer.

   **Additional Options**

   In case downloading the complete datasets is inconvenient, we also offer several demo datasets. These can be found in our `Test_Files/Demo` folder. Feel free to modify the file paths in our sample code or rename these demo datasets to suit your needs.

2. **Pathway datasets**

   IPFMC is a multi-omics integrated clustering method based on iterative fusion of pathways, and we provide information on the correspondence between processed pathways and genes and the correspondence between pathways and mirnas that we used in our experiments. These files can be found in ‘IPFMC/Test_Files’ folder, containing two neccessary pathway files:
   (1) Pathway_Index.csv
   (2) miRNA_Pathway_Index.csv

With the above data, you need to build the following file structure to easily replicate our method:

```python
.
├── Omics
│   └── LUAD
│       ├── LUAD_mRNA.csv
│       ├── LUAD_miRNA.csv
│       ├── LUAD_Methy.csv
│       └── LUAD_CNV.csv
├── Pathways
│   ├── Pathway_Index.csv
│   └── miRNA_Pathway_Index.csv
└── script.py
```

The python script currently in use for test is script.py. **To ensure the success of the following steps, please build this file structure.** 

### Step 2: Import neccessary packages

To run the following code, you need to install pandas, numpy, snfpy, and IPFMC packages in advance by typing the following command in the terminal:

```python
pip install pandas numpy snfpy IPFMC
```

Then use the following lines of code to import the package:

```python
import pandas as pd
import numpy as np
from snf import snf
from IPFMC import direct,separate,analysis
```

### Step 3: Input datasets

Add the following lines to the script to import the datasets correctly:

```python
# Filepath of the omics data, ‘LUAD’ is the folder contains omics datas of LUAD cancer
Omic_dir = './Omics/LUAD'  
# Filepath of the pathway index
BP_dir = './Pathways/Pathway_Index.csv'
# Filepath of the miRNA pathway index
mirBP_dir = './Pathways/miRNA_Pathway_Index.csv'
datatypes = ['mRNA','Methy','CNV','miRNA']  # The type of data to be used in the experiment
omic_list = []  # A list for storing multiple omics data
BP_data = pd.read_csv(BP_dir,index_col=0)  # The pandas package is used to pass in the pathway data
mirBP_data = pd.read_csv(mirBP_dir,index_col=0)  # Pass in the pathway-mirna relationship data
for datatype in datatypes:
    omicdata = pd.read_csv(f'{Omic_dir}/LUAD_{datatype}.csv',index_col=0)
    omic_list.append(omicdata)
```

### Step 4: Acquisition of single/multi-omics data representation

After obtaining all the necessary data, we can input them into IPFMC for multi-omics data integration. This will produce the multi-omics integrated representation and the ranking of the filtered retained pathway for each omics. In this step, IPFMC offers two modalities, each with two strategies. We use strategy 1 of IPFMC as an example to illustrate its usage (The corresponding method in the IPFMC package is 'ipfmc discretize()'). We showed two approaches (direct integration and separate computation) to obtain the multi-omics representation.

#### Approach 1: directly input the multi-omics data list and obtain the multi-omics representation 

You can choose to use a direct multi-omics integration strategy. Here's the code (The ‘omic_list’, ‘BP_data’ and ‘mirBP_data’ variable obtained earlier are used in this step):

```python
represent, pathways = direct.ipfmc_discretize(omic_list,BP_data,mirna=True,mirtarinfo=mirBP_data)
"""
	represent: Integrated representation of multi-omics data calculated by IPFMC
	pathways: The pathway ranking of each omics calculated by IPFMC (each omics has a pathway ranking), in the same order as the order of the omics in the input omic_list. 
"""
```

Detailed Parameters of ‘direct.ipfmc_discretize()’ are listed below:

```python
"""
    :param datasets: List of your multi-omics datasets, each element of the list should be a pandas dataframe.
    :param pathwayinfo: Pathways and their containing genetic information.
    :param k: The number of initial points of kmeans clustering
    :param fusetime: Number of pathway screening and fusion performed
    :param proportion: The proportion of pathways that were retained was fused at each iteration
    :param snfk: Number of SNF neighborhoods when multiple data sets are fused
    :param seed: Random number seed, set to None if no seed is needed
    :param mirtarinfo: miRNA target gene information, valid only if miRNA data is included in the dataset
    :param mirna: Set to True if your dataset contains mirna data, and False otherwise
    :return: Final representation of input datasets; a list of pathway rankings of each dataset.
"""
```

By default, all variables except `datasets`, `pathwayinfo`, `mirtarinfo`, and `mirna` have preset values and do not need to be manually configured. We set the `seed` parameter to 10 by default because our method involves k-means clustering and spectral clustering, both of which are influenced by random factors, such as the initial point selection in k-means clustering. Setting the seed to 10 helps ensure that you can reproduce results similar to those in our paper. While results may still vary due to differences in Python and package versions, they should be comparable. Alternatively, you can set the seed to None or any other number to run our method. The resulting data will differ, but the performance should be comparable to the results showcased in our paper.

**If your datasets contains miRNA expression data, please make sure the ‘mirna’ parameter is set to ‘True’, and the miRNA expression data must be the last element of ‘omic_list’ variable, ‘mirtarinfo’ must be set to the variable that contains miRNA-pathway relationship data.**

#### Approach 2: Compute the representation of each single omics separately

You can also choose to obtain single omics representation for each omics and then using SNF integration. 

```python
represents = []
pathways_list = [] # A list to store pahtway rankings of each omics
# Only the first three data sets are processed here, and the last data set is miRNA, which needs to be processed separately
for i in range(3):  
    represent, pathways = separate.ipfmc_discretize(omic_list[i], BP_data)
    represents.append(np.array(represent))
    print(represent)
    pathways_list.append(pathways)

represent, pathways = separate.ipfmc_discretize(omic_list[3], mirBP_data)  # Here processes miRNA dataset
represents.append(np.array(represent))
pathways_list.append(pathways)
represent_final = snf(represents, K=15)  # 'represent_final' is the final multi-omics representation
# print(pathways_list)
```

We recommend using this approach because computing the representation of each single-omics separately is more flexible in performing downstream tasks and has fewer parameters to consider.

### Step 5: Clustering using multi-omics representation

You can directly select number of clusters and use the code below to obtain cluster labels:

```python
labels = separate.spec_cluster(omic_list[0],fusion_matrix=represent_final,k=4)  # Here we set number of clusters to 4
# 'labels' is the cluster labels of input multi-omics datasets.
```

(The first parameter can be any element in ‘omic_list’. It is used to retrieve the sample name)

Or you can use the function we provide to recommend a suggested number of clusters.

```python
K = separate.suggest_k(represent_final)  # input the final representation, and this function will give a suggested cluster
labels = separate.spec_cluster(omic_list[0],fusion_matrix=represent_final,k=K)
```

Then you can use the obtained cluster labels to perform all kinds of analysis.

### (Optional) Step 6: Compute Gene occurence

This is an analysis covered in our paper, and we also provide the corresponding implementation method in the IPFMC package.

1. For datasets using genes as features, the top 50 genes with frequency of occurrence can be counted

```python
Gene_names = list(omic_list[0].index)  # obtain all gene names occured in the omics dataset
Gene_rank = analysis.gene_occurrence(Full_Genes=Gene_names, Sel_Pathways=pathways_list[0], Full_Pathways=BP_data)
# print(Gene_rank)
```

2. For datasets using miRNA as features, the top 50 miRNA with frequency of occurrence can be counted

```python
miRNA_names = list(omic_list[3].index)  # obtain all miRNA names occured in the omics dataset
miRNA_rank = analysis.gene_occurrence(Full_Genes=miRNA_names, Sel_Pathways=pathways_list[3], Full_Pathways=mirBP_data)
# print(miRNA_rank)
```
## Apply to your own datasets

To apply IPFMC to your own datasets, you can check the data format of the sample file we provided in the previous section and convert your own data to the appropriate format.

For pathway information data, you can either use the constructed data that we applied in our experiments (in the ‘IPFMC/Test_Flies/Pathways/’) or convert your own downloaded pathway information into the format we used with the user-friendly methods in the IPFMC package. However, we currently only support miRNA target gene data from mirTarBase([miRTarBase: the experimentally validated microRNA-target interactions database (cuhk.edu.cn)](https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/download.php)) and pathway data from GSEA([https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2](#C2)). 

The following is an explanation of the imported files that appear in the example code below.

(1) GoTerm_NCP.json: pathway information file downloaded from GSEA in the json format.

(2) hsa_MTI.xlsx: miRNA target gene information file downloaded from mirTarBase.

In addition, the ‘All_Names’ variable that appears in the code is the set of all miRNA names that you need to construct the relationship between miRNA and pathway. You can build your own code to construct this variable according to the characteristics of your data.


```python
from IPFMC import analysis
import pandas as pd
import pickle
#Cancer_type = 'LUAD'
miRNA_dir = '../Datas/Omics/Raw_Omics/'
BP_dir = '../Datas/Pathways/'
Target_dir = '../Datas/Pathways/hsa_MTI.xlsx'
BioProcesses = pd.read_json(BP_dir+'GoTerm_NCP.json')
# 'pathway_index' is constructed pathway-gene relationship data, you can save it as a file 
pathway_index = analysis.Create_pathwayindex(BioProcesses,min_genes=10,max_genes=200)
# pandas reads target gene data, which are obtained from miRTarBase
target_gene_data = pd.read_excel(Target_dir)
print(target_gene_data)
All_Names = set() # All_Names variable is here
Cancer_type_list = ['ACC','BRCA','COAD','KIRC','KIRP','LIHC','LUAD','LUSC','THYM']
for Cancer_type in Cancer_type_list:
    Temp_miRNA = pd.read_csv(miRNA_dir+Cancer_type+"/"+Cancer_type+"_miRNA.csv", index_col=0)
    Names = set(list(Temp_miRNA.index))
    All_Names = All_Names.union(Names)
    print(Names)
print(All_Names)
miRNA_Genes_dict = analysis.create_mirTar_Dict(target_gene_data,All_Names)

print(pathway_index)
'''
f = open("../Datas/Pathways/miRNA_Genes_dict.pkl", "rb")  
miRNA_Genes_dict = pickle.load(f)  
f.close() 
'''
print(len(miRNA_Genes_dict))  
# 'mir_index' is final constructed miRNA-pathway relationship data, you can save it as a file 
mir_index = analysis.Create_miRNA_pathwayindex(miRNA_Genes_dict,pathway_index)
print(mir_index)
```

## Evaluation codes of our experiments

You can check ‘IPFMC/IPFMC_working_directory’for all our evaluation code. 

