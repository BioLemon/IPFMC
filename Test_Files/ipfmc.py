import argparse
import pandas as pd
import numpy as np
from snf import snf
from IPFMC import separate, direct, analysis
import sys
import pickle
import pandas as pd

def read_omics(omics_file):
    """读取组学数据文件，返回组学类型和对应的数据路径字典"""
    with open(omics_file, 'r') as f:
        lines = f.readlines()

    # 读取组学类型
    omics_types = lines[0].strip().split()
    omics_paths = [line.strip() for line in lines[1:]]

    # 确保组学类型和路径一一对应
    if len(omics_types) != len(omics_paths):
        raise ValueError("组学类型和路径数量不匹配。")

    # 检查组学类型数量是否少于2个
    if len(omics_types) < 2:
        raise ValueError("组学类型数量少于2个，请提供至少两个组学类型。")

    # 创建一个字典来保存组学数据
    omics_data = {}
    for omics_type, path in zip(omics_types, omics_paths):
        omics_data[omics_type] = pd.read_csv(path, index_col=0)  # 假设CSV文件格式

    return omics_data

def read_pathways(pathway_file):
    """读取通路信息文件，返回pathway_index和miRNA_pathway_index"""
    with open(pathway_file, 'r') as f:
        lines = f.readlines()

    # 读取前两行路径
    if len(lines) < 2:
        raise ValueError("通路文件至少需要两行。")

    pathway_index_path = lines[0].strip()
    mirna_pathway_index_path = lines[1].strip()

    # 读取对应的CSV文件
    pathway_index = pd.read_csv(pathway_index_path,index_col=0)
    mirna_pathway_index = pd.read_csv(mirna_pathway_index_path,index_col=0)

    return pathway_index, mirna_pathway_index


def cluster_st1(omics_data, pathway_index, mirna_pathway_index, output_dir, cancer_type, num_cluster=None, pathway_prop=100):
    """聚类算法的实现（策略1）"""
    # 检查是否存在miRNA组学数据
    if 'miRNA' in omics_data:
        mirna_data = omics_data.pop('miRNA')  # 提取miRNA数据并从omics_data中删除
        print("提取miRNA组学数据，处理其他组学数据")
    else:
        mirna_data = None
        print("未找到miRNA组学数据")

    # 处理其他组学数据
    omic_list = [data for data in omics_data.values()]  # 获取剩余组学数据列表
    represents = []
    pathways_list = {}

    # 处理每个组学数据
    for omic_type, omic in omics_data.items():
        represent, pathways = separate.ipfmc_discretize(omic, pathway_index,preselect=int(pathway_prop))  # 使用策略1
        represents.append(np.array(represent))
        pathways_list[omic_type] = pathways
        pathways_df = pd.DataFrame(pathways)
        pathways_df.to_csv(f"{output_dir}/{cancer_type}_{omic_type}_pathways.csv", index=False)

    # 处理miRNA数据（如果需要）
    if mirna_data is not None:
        represent, pathways = separate.ipfmc_discretize(mirna_data, mirna_pathway_index,preselect=int(pathway_prop))  # 使用策略2处理miRNA数据
        represents.append(np.array(represent))
        pathways_list['miRNA'] = pathways
        pathways_df = pd.DataFrame(pathways)
        pathways_df.to_csv(f"{output_dir}/{cancer_type}_miRNA_pathways.csv", index=False)

    # 最终聚合表示
    represent_final = snf(represents, K=15)  # K值可以根据需要调整

    # 使用建议的聚类数
    if num_cluster is None:
        K = separate.suggest_k(represent_final)  # 获取建议的聚类数
    else:
        K = num_cluster

    labels = separate.spec_cluster(omic_list[0], fusion_matrix=represent_final, k=K)  # 聚类标签
    labels.columns = ["cluster"]

    return labels, K, pathways_list

def cluster_st2(omics_data, pathway_index, mirna_pathway_index, output_dir, cancer_type, num_cluster=None):
    """聚类算法的实现（策略2）"""
    # 检查是否存在miRNA组学数据
    if 'miRNA' in omics_data:
        mirna_data = omics_data.pop('miRNA')  # 提取miRNA数据并从omics_data中删除
        print("提取miRNA组学数据，处理其他组学数据")
    else:
        mirna_data = None
        print("未找到miRNA组学数据")

    # 处理其他组学数据
    omic_list = [data for data in omics_data.values()]  # 获取剩余组学数据列表
    represents = []
    pathways_list = {}

    # 处理每个组学数据
    for omic_type, omic in omics_data.items():
        represent, pathways = separate.ipfmc_average(omic, pathway_index)  # 使用策略2
        represents.append(np.array(represent))
        pathways_list[omic_type] = pathways
        pathways_df = pd.DataFrame(pathways)
        pathways_df.to_csv(f"{output_dir}/{cancer_type}_{omic_type}_pathways.csv", index=False)

    # 处理miRNA数据（如果需要）
    if mirna_data is not None:
        represent, pathways = separate.ipfmc_average(mirna_data, mirna_pathway_index)  # 使用策略2处理miRNA数据
        represents.append(np.array(represent))
        pathways_list['miRNA'] = pathways
        pathways_df = pd.DataFrame(pathways)
        pathways_df.to_csv(f"{output_dir}/{cancer_type}_miRNA_pathways.csv", index=False)

    # 最终聚合表示
    represent_final = snf(represents, K=15)  # K值可以根据需要调整

    # 使用建议的聚类数
    if num_cluster is None:
        K = separate.suggest_k(represent_final)  # 获取建议的聚类数
    else:
        K = num_cluster

    labels = separate.spec_cluster(omic_list[0], fusion_matrix=represent_final, k=K)  # 聚类标签
    labels.columns = ["cluster"]

    return labels, K, pathways_list

def create_pathway_index(BP_dir):
    Bio_Pathways = pd.read_json(BP_dir)
    pathway_index = analysis.Create_pathwayindex(Bio_Pathways, min_genes=10, max_genes=200)
    print("Completed constructing pathway index.")
    return pathway_index
def create_miRTar_index(mirTar_dir, mir_set_dir,pathway_dir):
    mirTar_info = pd.read_excel(mirTar_dir)
    mirset = pd.read_csv(mir_set_dir,header=None)
    Bio_Pathways = pd.read_csv(pathway_dir,index_col=0)
    print(Bio_Pathways)
    print(mirset)
    mirset = mirset[0].to_list()
    mirTar_index = analysis.create_mirTar_Dict(mirTar_info,mirset)
    mirPathway_index = analysis.Create_miRNA_pathwayindex(mirTar_index,Bio_Pathways)
    return mirPathway_index

def main():
    parser = argparse.ArgumentParser(description="处理组学数据的命令行工具")

    # 添加简化的参数
    parser.add_argument('-p', '--pathway', required=True, help='路径到通路信息文件')
    parser.add_argument('-m', '--mirset', default=None, help='miRNA集合文件路径')
    parser.add_argument('-o', '--omics', default=None, help='路径到组学信息文件')
    parser.add_argument('-a', '--action', type=str, choices=['c', 'p', 'm'], required=True, help='执行的操作: c, p, m')
    parser.add_argument('-s', '--strategy', type=int, default=1, help='策略选择: 1或2')
    parser.add_argument('-n', '--num_cluster', type=int, help='聚类数量，默认为建议数量')
    parser.add_argument('-d', '--output', type=str, default='./results', help='输出目录，默认为当前目录下的results文件夹')
    parser.add_argument('-co', '--count', type=int, default=0,
                        help='输出目录，默认为当前目录下的results文件夹')
    parser.add_argument('-can', '--cancer',type=str,default="Cancer",help="Cancer type of input datasets")
    parser.add_argument('-se', '--select', default=100, help='Proportion of retained pathways')
    args = parser.parse_args()

    # 读取参数
    pathway_file = args.pathway
    omics_file = args.omics
    action = args.action
    strategy = args.strategy
    num_cluster = args.num_cluster
    output_dir = args.output
    compute_count = args.count
    cancer_type = args.cancer
    mir_set_file = args.mirset
    pathway_proportion = args.select

    if strategy != 1 and strategy != 2:
        print("请设置strategy为1或2！")
        sys.exit(1)  # 终止程序，返回状态码1表示错误
    # 根据 action 执行不同的操作
    if action == "c":
        omics_data = read_omics(omics_file)  # 读取组学数据
        origin_omicsdata = omics_data.copy()
        pathway_index, mirna_pathway_index = read_pathways(pathway_file)  # 读取通路数据
        if strategy == 1:
            Labels,K,Pathways_list = cluster_st1(omics_data, pathway_index, mirna_pathway_index, output_dir, cancer_type, num_cluster, pathway_proportion)
            Labels.to_csv(f"{output_dir}/{cancer_type}_{K}.csv")
        elif strategy == 2:
            Labels, K, Pathways_list = cluster_st2(omics_data, pathway_index, mirna_pathway_index, output_dir, cancer_type, num_cluster)
            Labels.to_csv(f"{output_dir}/{cancer_type}_{K}.csv")
        if compute_count == 1:
            for omic_type, omic in origin_omicsdata.items():
                if omic_type != "miRNA":
                    Full_Genes = list(omic.index)
                    Gene_counts = analysis.gene_occurrence(Full_Genes,Pathways_list[omic_type],pathway_index)
                    Gene_counts.to_csv(f"{output_dir}/{cancer_type}_{omic_type}_counts.csv")
                else:
                    Full_Genes = list(omic.index)
                    Gene_counts = analysis.gene_occurrence(Full_Genes, Pathways_list[omic_type], mirna_pathway_index)
                    Gene_counts.to_csv(f"{output_dir}/{cancer_type}_{omic_type}_counts.csv")
    elif action == "p":
        pathway_index = create_pathway_index(pathway_file)
        pathway_index.to_csv(f"{output_dir}/Pathway_index.csv")
    elif action == "m":
        mirna_index = create_miRTar_index(omics_file,mir_set_file,pathway_file)
        mirna_index.to_csv(f"{output_dir}/miRNA_Pathway_index.csv")
    else:
        print("请提供正确的操作")

if __name__ == "__main__":
    main()





