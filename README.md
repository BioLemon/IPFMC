# IPFMC
<P align="center">
    <a>
    <img alt="Licence" src="https://img.shields.io/pypi/l/IPFMC" />
	</a>
    <a>
    <img alt="Licence" src="https://img.shields.io/pypi/v/IPFMC" />
	</a>
    <a>
    <img alt="Licence" src="https://img.shields.io/pypi/pyversions/IPFMC" />
	</a>
</p>

## Overview

This README file serves as a comprehensive guide for understanding and utilizing our IPFMC method. It is divided into three main sections to address different user needs and provide clear instructions:

### Section 1: IPFMC Simple Use Case

This section offers a straightforward tutorial on running our method using command-line interface. It provides guidance on setting parameters to meet your specific requirements, ensuring a seamless integration of IPFMC into your workflow.

### Section 2: Detailed Usage of IPFMC

For users who need to incorporate our method into their own codes, this section delves into the intricacies of using IPFMC as a general Python package. It offers comprehensive instructions and examples to help you seamlessly integrate our tool with your existing codebase.

### Section 3: Instructions for Paper Data Reproduction

This section is dedicated to reproducing all experimental results presented in our research paper. By following the instructions provided, you can validate our findings and explore the performance of IPFMC on the same datasets used in our study.

**Note: Please ensure that you update IPFMC to the latest version before using it, as this guidance is based on the most recent version. Otherwise, you may encounter some errors.**

We hope this README file serves as a valuable resource for understanding and leveraging the capabilities of our IPFMC method.
## Section 1: IPFMC Simple Use Case

### Create Environment

Download the IPFMC/Test_Files, and install the Python packages listed in the requirements.txt file. **Ensure that the ipfmc package is updated to the latest version**.

### Run on Demo datasets

You can run the Demo dataset by typing the following command from a python terminal, and all results will be saved by default in the Test_files/results folder:

```python
python ipfmc.py -p pathway_info.txt -o omics_info.txt -a c -s 1
```

In this example command, `-p` should be set to a text file containing the file paths of the pathway index and the miRNA pathway index, and the order cannot be changed. The `-o` parameter needs to be set to a text file where the first line contains the required input omics types separated by spaces, and each subsequent line contains the file path for one type of omics data (the order of the omics file paths must correspond to the order of the input omics types). 

You can also modify the parameters of our method to perform analysis according to your needs, following is an introduction of all available parameters.

#### Optional Parameters

- `-p, --pathway`: Path to the pathway information file. 
  
- `-m, --mirset`: Path to the miRNA set file. If provided, the tool will use the miRNA set to build the miRNA-pathway index.
  
- `-o, --omics`: Path to the txt file containing omics’ file path information. 

- `-a, --action`: Specifies the action to be performed. Acceptable values are 'c', 'p', and 'm', each representing different operations. 
  - 'c' indicates the clustering mode, which performs clustering analysis using the provided multi-omics data.
  - 'p' indicates the pathway index construction mode, which is used to build a pathway index similar to `Test_Files/Pathways/Pathway_Index.csv` based on the pathway information file from MSigDB. In this mode, the `-p` parameter should be set to the path of the pathway information JSON file downloaded from MSigDB.
  - 'm' indicates the construction of the relationship between miRNAs and pathways, similar to `Test_Files/Pathways/miRNA_Pathway_Index.csv`. In this case, you need to set the `-p` parameter to the path of `Pathway_Index.csv`, the `-m` parameter to a CSV file containing all miRNAs for which the relationships are to be constructed (similar to `Test_Files/Pathways/raw/mirset_example.csv`), and the `-o` parameter to the miRNA target gene information file downloaded from miRTarBase, such as `Test_Files/Pathways/raw/hsa_MTI.xlsx`.

- `-s, --strategy`: Specifies whether to use strategy one or strategy two to perform IPFMC (Integrated Pathway and Functional Modulation Clustering). Acceptable values are 1 or 2, with a default of 1.

- `-n, --num_cluster`: Specifies the number of clusters. If not provided, the tool will use a suggested number of clusters.

- `-d, --output`: Specifies the output directory. Defaults to a folder named 'results' in the current directory.

- `-co, --count`: Determines whether to output gene count rankings. The default value is 0, which means no output. Setting it to 1 will output the gene count rankings for important pathways.

- `-can, --cancer`: Specifies the cancer type of the input datasets. Defaults to 'Cancer'. This parameter only affects the names of the output files.

- `-se, --select`: Specifies the proportion of retained pathways. Defaults to 100. If set to a lower proportion, it will filter a subset of pathways with higher absolute median deviations, thereby increasing the speed of execution, but it may affect the clustering results. We recommend that this proportion should not be set below 40%.

### Example of running with m and p mode

1. building the pathway index

```python
python ipfmc.py -p ./Pathways/raw/GoTerm_NCP.json -a p
```

2. building the miRNA pathway index

Due to the large size of the raw miRNA target data, you will need to download it from the following link:
https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/hsa_MTI.xlsx 

After downloading, place the file in the `Test_Files/Pathways/raw` directory. Once this is done, you can execute the following command to build the miRNA-pathway index.
```python
python ipfmc.py -p ./Pathways/Pathway_index.csv -o ./Pathways/raw/hsa_MTI.xlsx -m ./Pathways/raw/mirset_example.csv -a m
```

Index obtained by above commands will be stored in `results` folder by default. You can set `-d` parameter to specify the output directory.

## Section 2: Detailed Usage of IPFMC package

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
represent, pathways = direct.ipfmc_discretize(omic_list,BP_data,mirna=True,mirtarinfo=mirBP_data,seed=10)
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

By default, all variables except `datasets`, `pathwayinfo`, `mirtarinfo`, and `mirna` have preset values and do not need to be manually configured. We set the `seed` parameter to 10 in the example codes because our method involves k-means clustering and spectral clustering, both of which are influenced by random factors, such as the initial point selection in k-means clustering. Setting the seed to 10 helps ensure that you can reproduce results similar to those in our paper. While results may still vary due to differences in Python and package versions, they should be comparable. Alternatively, you can set the seed to None or any other number to run our method. The resulting data will differ, but the performance should be comparable to the results showcased in our paper.

**If your datasets contains miRNA expression data, please make sure the ‘mirna’ parameter is set to ‘True’, and the miRNA expression data must be the last element of ‘omic_list’ variable, ‘mirtarinfo’ must be set to the variable that contains miRNA-pathway relationship data.**

#### Approach 2: Compute the representation of each single omics separately

You can also choose to obtain single omics representation for each omics and then using SNF integration. 

```python
represents = []
pathways_list = [] # A list to store pahtway rankings of each omics
# Only the first three data sets are processed here, and the last data set is miRNA, which needs to be processed separately
for i in range(3):  
    represent, pathways = separate.ipfmc_discretize(omic_list[i], BP_data,seed=10)
    represents.append(np.array(represent))
    print(represent)
    pathways_list.append(pathways)

represent, pathways = separate.ipfmc_discretize(omic_list[3], mirBP_data,seed=10)  # Here processes miRNA dataset
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

In addition, the ‘miRNA_Names’ variable that appears in the code is the set of all miRNA names that you need to construct the relationship between miRNA and pathway. You can build your own code to construct this variable according to the characteristics of your data.


```python
from IPFMC import analysis
import pandas as pd
import pickle

# Define directories for data sources
miRNA_dir = '../Datas/Omics/Raw_Omics/'
BP_dir = '../Datas/Pathways/'
Target_dir = '../Datas/Pathways/hsa_MTI.xlsx'

# Load pathway data from JSON file
Pathway_data = pd.read_json(BP_dir + 'GoTerm_NCP.json')

"""
(1) Constructs a pathway-gene relationship index. The resulting 'pathway_index' can be saved as a file for future use.
(2) The parameters 'min_genes' and 'max_genes' define the range of gene counts to retain for pathways.
"""
pathway_index = analysis.Create_pathwayindex(Pathway_data, min_genes=10, max_genes=200)

# Read target gene data from an Excel file sourced from miRTarBase
target_gene_data = pd.read_excel(Target_dir)

"""
The following code extracts a comprehensive set of miRNA names that require mapping to pathways. 
It specifically retrieves indices from all cancer miRNA omics data and stores all unique values. 
Alternative methods may also be employed to obtain the set of miRNA names.
"""
miRNA_Names = set()
Cancer_type_list = ['ACC', 'BRCA', 'COAD', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'THYM']

for Cancer_type in Cancer_type_list:
    Temp_miRNA = pd.read_csv(miRNA_dir + Cancer_type + "/" + Cancer_type + "_miRNA.csv", index_col=0)
    Names = set(list(Temp_miRNA.index))
    miRNA_Names = miRNA_Names.union(Names)

# Create a dictionary that maps miRNAs to their corresponding target genes.
miRNA_Genes_dict = analysis.create_mirTar_Dict(target_gene_data, miRNA_Names)

# The 'mir_index' represents the final constructed miRNA-pathway relationship index, which can be saved as a file for further analysis.
mir_index = analysis.Create_miRNA_pathwayindex(miRNA_Genes_dict, pathway_index)

# Output the final miRNA-pathway index
print(mir_index)
```

## Section 3: Instructions for Paper Data Reproduction

### Step 1: Download Evaluation Codes and Omics Datasets

(1) Download Evaluation Codes

First, please download all files from the `IPFMC/IPFMC_working_directory/` directory in our repository to automatically obtain the code and folder structure required to reproduce our paper data.

(2) Download Cancer Multi-Omics Datasets

Next, you need to acquire all the evaluation omics datasets we used.

Our omics datasets are derived from the comprehensive study conducted by Duan et al., as published in *PLOS Computational Biology* in 2021. Their research evaluated and compared multi-omics data integration methods for cancer subtyping. You can find their original paper here: ([Evaluation and comparison of multi-omics data integration methods for cancer subtyping | PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009224)).

Their datasets can be downloaded from their github repository: [GaoLabXDU/MultiOmicsIntegrationStudy: Experimental evaluation and comparison of multi-omics data integration methods for cancer subtyping (github.com)](https://github.com/GaoLabXDU/MultiOmicsIntegrationStudy/).

After accessing their github repository, go to the `Availability of datasets` section in their README file. Download all Complete datasets. These datasets have a unified naming convention: "Dataset #1 XXXX Complete". Here, XXXX represents the cancer type, including BRCA, COAD, KIRC, LUAD, LUSC, ACC, KIRP, LIHC, and THYM. Please download the Complete datasets for all nine cancer types, as these datasets are used in the survival analysis experiments in our paper. Additionally, you need to download the two most recent datasets: Dataset #3 BRCA Complete and Dataset #3 COAD Complete, which are used in the gold standard dataset evaluation experiments in our paper.

After downloading, each dataset contains multiple omics files. You need to place these omics files in the respective cancer folders under the `IPFMC_working_directory/Datas/Omics/Raw_Omics` directory, matching the cancer types. The placeholder files already present in the folders can be deleted or left as is, as they will not affect the reproduction process.

Specifically, the two gold standard evaluation datasets (Dataset #3 BRCA Complete and Dataset #3 COAD Complete) should be placed in the corresponding cancer+GOLD named folders.

The remaining data, such as pathway information and patient survival data, are already included in `IPFMC/IPFMC_working_directory` and do not require separate downloading or processing.

### Step 2: Deploy to the server and create  ipfmc conda environment

Upload the entire working directory, which now contains the downloaded and processed omics data, to a Linux server (preferably one with more than four CPUs). Subsequently, create and activate the `ipfmc` Conda environment using the `ipfmc.yml` file with the following commands:

```bash
conda env create -f ipfmc.yml
conda activate ipfmc
```

Once the environment is created, you can proceed to the next step.

### Step 3:  Run evaluations

Next, set the current directory to the `Evaluation_bash` folder and execute the following steps in order:

1. Run strategy 1 to obtain the four omics representations and pathway rankings for each of the nine cancer types:

```bash
bash Evaluation_S1.sh
```

This step is time-consuming (about 1 hour running on four Intel(R) Xeon(R) Platinum 8375C CPU @ 2.90GHz). Please monitor the background command execution and ensure that all four programs activated by the `sh` file have completed before executing the next step.

2. Run strategy 2 to obtain the four omics representations and pathway rankings for the two gold standard datasets:

```bash
bash Evaluation_S2.sh
```
As in the first step, wait for this step to complete before proceeding to the next step.

3. Run the survival analysis and gold standard comparison analysis to obtain the survival analysis results for nine cancer datasets and true label comparison experiment results for the two gold standard datasets:

```bash
bash Downana.sh
```

After completing the above three steps, you can find all the evaluation data from our paper in `IPFMC_working_directory/Datas/Assessment_Result`.
