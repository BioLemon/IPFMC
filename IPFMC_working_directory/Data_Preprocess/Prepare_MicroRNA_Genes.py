import pandas as pd

# pd.set_option('display.max_rows', None)
Cancer_type_list = ['ACC','BRCA','COAD','KIRC','KIRP','LIHC','LUAD','LUSC','THYM']
miRNA_dir = '../Datas/Omics/Raw_Omics/'
Target_dir = '../Datas/Pathways/hsa_MTI.xlsx'
Output_dir = '../Datas/Pathways/'

target_gene_data = pd.read_excel(Target_dir)
print(target_gene_data)

All_Names = set()

for Cancer_type in Cancer_type_list:
    Temp_miRNA = pd.read_csv(miRNA_dir+Cancer_type+"/"+Cancer_type+"_miRNA.csv", index_col=0)
    Names = set(list(Temp_miRNA.index))
    All_Names = All_Names.union(Names)
    print(Names)
print(All_Names)

Let_list = []
Mir_list = []

LET_PREFIX = 'hsa-let'
MIR_PREFIX = 'hsa-mir'
PREFIX_LEN = 7
for miRNA in All_Names:
    if miRNA[:PREFIX_LEN] == LET_PREFIX:
        Let_list.append(miRNA)
    elif miRNA[:PREFIX_LEN] == MIR_PREFIX:
        Mir_list.append(miRNA)

target_gene_dict = {}
for i in range(len(target_gene_data['miRNA'])):
    MiR = target_gene_data['miRNA'][i]
    Gene = target_gene_data['Target Gene'][i]
    if MiR not in target_gene_dict:
        target_gene_dict[MiR] = set()
    target_gene_dict[MiR].add(Gene)
print(len(target_gene_dict))
print(target_gene_dict['hsa-let-7a-5p'])
print(len(target_gene_dict['hsa-let-7a-5p']))

miRNA_Genes_dict = {}
for miRNA in Let_list:
    Genes = set()
    if len(miRNA) == 10:
        Genes = Genes.union(target_gene_dict.get(miRNA+'-5p', set()))
        Genes = Genes.union(target_gene_dict.get(miRNA+'-3p', set()))
        Genes = Genes.union(target_gene_dict.get(miRNA, set()))
    elif miRNA[10:] == '-1':
        Genes = Genes.union(target_gene_dict.get(miRNA[:10]+'-5p', set()))
        Genes = Genes.union(target_gene_dict.get(miRNA[:10]+'-3p', set()))
        Genes = Genes.union(target_gene_dict.get(miRNA, set()))
    else:
        Genes = Genes.union(target_gene_dict.get(miRNA[:10]+'-5p', set()))
        Genes = Genes.union(target_gene_dict.get(miRNA+'-3p', set()))
        Genes = Genes.union(target_gene_dict.get(miRNA, set()))
    miRNA_Genes_dict[miRNA] = Genes


for miRNA in Mir_list:
    Genes = set()
    mature_miRNA_5p = miRNA.replace("mir", "miR") + "-5p"
    mature_miRNA_3p = miRNA.replace("mir", "miR") + "-3p"
    mature_miRNA = miRNA.replace("mir", "miR")
    Genes = Genes.union(target_gene_dict.get(mature_miRNA_5p, set()))
    Genes = Genes.union(target_gene_dict.get(mature_miRNA_3p, set()))
    Genes = Genes.union(target_gene_dict.get(mature_miRNA, set()))
    miRNA_Genes_dict[miRNA] = Genes
    #print(miRNA)
    #print(Genes)

print(len(miRNA_Genes_dict))
for miRNA in miRNA_Genes_dict:
    Genes = set()
    if miRNA_Genes_dict[miRNA] == set():
        mature_miRNA = miRNA.replace("mir", "miR")
        mature_miRNA = mature_miRNA[:-2]
        #print(mature_miRNA)
        Genes = Genes.union(target_gene_dict.get(mature_miRNA + '-5p', set()))
        Genes = Genes.union(target_gene_dict.get(mature_miRNA + '-3p', set()))
        Genes = Genes.union(target_gene_dict.get(mature_miRNA, set()))
        #print(Genes)
        miRNA_Genes_dict[miRNA] = Genes
#print(miRNA_Genes_dict)

