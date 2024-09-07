# IPFMC

## Brief description of each module

This section provides a brief description of each module, for detailed description of parameters of each method, see function descriptions.

### ‘direct’ module

This module provides methods to directly perform integrated cancer multi-omics clustering using IPFMC.

(1) **direct.ipfmc_discretize()**: Implementation of strategy 1 for IPFMC.

(2) **direct.ipfmc_average()**: Implementation of strategy 2 for IPFMC.

(3) **direct.spec_cluster()**: Generate spectral clustering results in cluster labels with sample indexes.

(4) **direct.suggest_k()**: Gives a suggested number of clusters according to the silhouette coefficient.

### ‘separate’ module

This module has roughly the same function as the direct module, but the two strategies of ipfmc accept single omics data and return the single omics representation and pathway ranking. Users can use similarity network fusion(SNF) to fuse each single omics representation to obtain multi-omics representation.

(1) **separate.ipfmc_discretize()**: Implementation of strategy 1 for IPFMC.

(2) **separate.ipfmc_average()**: Implementation of strategy 2 for IPFMC.

(3) **separate.spec_cluster()**: Generate spectral clustering results in cluster labels with sample indexes.

(4) **separate.suggest_k()**: Gives a suggested number of clusters according to the silhouette coefficient.

### ‘analysis’ module

This module provides some functions for pathway data processing and downstream analysis.

## Simple Test Case

This section provides sample code for multi-omics data integration clustering and biological interpretation using the package, and you can change some of the variables to apply it to your own dataset.

### Import neccessary packages

```python
import pandas as pd
import numpy as np
from snf import snf
from IPFMC import direct
from IPFMC import separate
```

### Input datasets

1. Omics data

   All standard input omics data should be a csv file with one feature in each row and one sample in each column. The first row should be the sample name and the first column should be the gene name. (For other omics data besides miRNA and mRNA expression data, such as methylation, copy number variation, etc., the features should be mapped to genes and converted to gene names before being used as IPFMC input data).

2. Pathway data

   In addition to omics data, it is also necessary to input the gene information data contained in the general pathway. If your omics data includes miRNA omics, you also need to input the corresponding relationship data between miRNA and pathway.

Code is as follows:

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
    '''
    We named the omics data <cancer name>_<data type>.csv, for example, LUAD_mRNA.csv
    You can change it according to your habits
    '''
    omicdata = pd.read_csv(f'{Omic_dir}/LUAD_{datatype}.csv',index_col=0)
    omic_list.append(omicdata)
```

The file structure used in the sample code is as follows:

```bash
.
├── Omics
│   └── LUAD
│       ├── LUAD_mRNA.csv
│       ├── LUAD_miRNA.csv
│       ├── LUAD_Methy.csv
│       └── LUAD_CNV.csv
└── Pathways
    ├── Pathway_Index.csv
    └── miRNA_Pathway_Index.csv
└── script.py
```

Where script.py is the python script currently in use. You can also personalize the data by changing the path of each file, but the key is to use the read_csv provided by pandas and make sure that the row index of omics data is the feature name, the column index is the sample name, and the row index of pathway data is the pathway name.

### Acquisition of single/multi-omics data representation

After obtaining all the necessary data, we can input them into IPFMC for multi-omics data integration. This will produce the multi-omics integrated representation and the ranking of the filtered retained pathway for each omics. In this step, IPFMC offers two modalities, each with two strategies. We use strategy 1 of IPFMC as an example to illustrate its usage. We showed two approaches (direct integration and separate computation) to obtain the multi-omics representation.

#### directly input the multi-omics data list and obtain the multi-omics representation 

You can choose to use a direct multi-omics integration strategy. This requires importing the direct module. Here's the code (The ‘omic_list’, ‘BP_data’ and ‘mirBP_data’ variable obtained earlier are used in this step):

```python
represent, pathways = direct.ipfmc_discretize(omic_list,BP_data,mirna=True,mirtarinfo=mirBP_data)
"""
	represent: Integrated representation of multi-omics data calculated by IPFMC
	pathways: The pathway ranking of each omics calculated by IPFMC (each omics has a path ranking), in the same order as the order of the omics in the input omic_list
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

**If your datasets contains miRNA expression data, please make sure the ‘mirna’ parameter is set to ‘True’, and the miRNA expression data must be the last element of ‘omic_list’ variable, ‘mirtarinfo’ must be set to the variable that contains miRNA-pathway relationship data.**

#### Compute the representation of each single omics separately

You can also choose to obtain single omics representation for each omics and then using SNF integration. 

```python
represents = []
pathways_list = []
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
```

We recommend using this approach because computing the representation of each single-omics separately is more flexible in performing downstream tasks and has fewer parameters to consider.

### Clustering using multi-omics representation

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
