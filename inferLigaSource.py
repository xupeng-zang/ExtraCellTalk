import warnings
import pandas as pd
import numpy as np
import calculateLRScore

warnings.filterwarnings("ignore")


"""
提供一个推测配体可能来源的选项，
根据配体蛋白丰度log2FC矩阵和输入的配体可能来源细胞基因表达矩阵，推测配体蛋白来源。

Provides an option to infer the possible source of the ligand, 
based on the ligand protein abundance log2FC matrix and 
the inputted ligand possible source cell gene expression matrix, 
infer the source of the ligand protein.
"""


def inferLigaSource(LR_pair,
                    liga_data,
                    liga_info,
                    liga_control,
                    liga_exper,
                    liga_source_data,
                    liga_source_info,
                    liga_source_control,
                    liga_source_exper,
                    liga_source_expr_cutoff):
    """
    利用配体在输入的可能来源细胞中的表达情况推测LR关系对中配体的来源。

    LR_pair: ExtraCellTalk计算得到的配体-受体关系对结果,
    liga_data: 配体蛋白丰度矩阵,
    liga_info: 配体蛋白样本分组信息,
    liga_control: 配体蛋白对照组名称,
    liga_exper: 配体蛋白实验组名称,
    liga_source_data: 配体可能来源细胞基因表达矩阵,
    liga_source_info: 配体可能来源细胞类型分组注释,
    liga_source_control: 配体可能来源细胞对照组名称,
    liga_source_exper: 配体可能来源细胞实验组名称,
    liga_source_expr_cutoff: 0.25,
        基因在该配体可能来源细胞中表达数量大于整体细胞数量*expr_cutoff时，
        才认为该基因在该类型细胞中表达，
        默认数值为 0.25。

    The source of the ligand in the LR relationship pair was inferred 
    by the expression of the ligand in the input potential source cells.

    LR_pair: the ligand-receptor relationship pairs calculated by ExtraCellTalk,
    liga_data: ligand protein abundance matrix,
    liga_info: ligand protein sample grouping information,
    liga_control: ligand protein control group name,
    liga_exper: ligand protein experimental group name,
    liga_source_data: ligand possible source cell gene expression matrix,
    liga_source_info: ligand possible source cell type grouping information,
    liga_source_control: ligand possible source cell control group name,
    liga_source_exper: ligand possible source cell experimental group name,
    liga_source_expr_cutoff: 0.25,
        The gene is considered to be expressed in ligand possible source cell 
        only when the cell numbers of gene expressed is greater than 
        the overall number of cells*expr_cutoff.
        The default value is 0.25.
    """

    # 计算两个分组中配体蛋白丰度的log2FC值
    # Calculate the log2FC value of ligand protein abundance in both groups.
    liga_mean = calculateLRScore.calculateProtLog2FC(liga_data,
                                                     liga_info,
                                                     liga_control,
                                                     liga_exper)

    # 计算配体可能来源细胞中两个分组基因表达的log2FC值
    # Calculation of log2FC values for the expression of 
    # two grouped genes in cells of possible origin of the ligand.
    liga_source_mean = \
        calculateLRScore.calculateGeneLog2FC(liga_source_data,
                                             liga_source_info,
                                             liga_source_control,
                                             liga_source_exper,
                                             liga_source_expr_cutoff)

    # 创建储存结果的标准格式数据框
    # Create a standard format dataframe to store the results
    res = pd.DataFrame(columns=['ligands','receptors','interaction_name',\
        'LR_score','celltype','Pvalue','FDR','liga_source','liga_source_score'])


    def _is_in_index(x):
        """
        判断配体是否存在于可能来源的细胞基因表达矩阵中。

        Determine whether the ligand is present 
        in the cell gene expression matrix of possible sources.
        """

        return x in liga_source_mean.index
    

    # 对LR关系对中的配体迭代
    # Iterate over ligands in LR relationship pairs
    for i in list(filter(_is_in_index,set(LR_pair['ligands']))):

        # 获取LR关系对中相同配体的不同关系对
        # Obtain different relationship pairs for the same ligand 
        # in the LR relationship pairs
        temp_LR_pair = LR_pair[LR_pair['ligands']==i]

        # 计算配体来源分数
        # Calculate ligand source score
        temp_score = \
            liga_source_mean.loc[i,'log2FC'] * liga_mean.loc[i,'log2FC']
        
        # 判断配体是否为单一细胞来源,
        # 一种细胞来源时需将自动转换的 pd.Series 类型数据转换为 pd.DataFrame
        # Determining whether the ligand is of single cell origin, 
        # convert the automatically converted pd.Series type data 
        # to pd.DataFrame when there is one cell origin.
        if (type(temp_score) == np.float_):
            
            # 转换 pd.Series 类型为 pd.DataFrame
            # Convert the pd.Series type data to pd.DataFrame
            temp = pd.DataFrame(liga_source_mean.loc[i].to_dict(),index=[0])
            
            # 筛选计算得到的配体来源分数，当分数大于0时代表可能配体来源
            # Screen the calculated ligand source score, when the ligand source 
            # score > 0, it means a possible ligand source.
            if temp_score > 0:
                temp_liga_source = \
                    pd.DataFrame(
                        {'liga_source':temp['celltype'].values,\
                         'liga_source_score':temp_score})
                temp_liga_source['ligands'] = i
                
                # 为LR关系对添加配体来源
                # Add ligand source for LR relationship pairs
                temp_res = pd.merge(temp_LR_pair,temp_liga_source,on='ligands')
                
                # 存储得到的临时结果
                # Store the temporary result obtained
                res = pd.concat([res,temp_res])
            
        else:
            temp_liga_source = \
                pd.DataFrame(
                    {'liga_source':\
                    liga_source_mean.loc[i,'celltype'][temp_score>0].values,
                    'liga_source_score':temp_score[temp_score>0].values})
            temp_liga_source['ligands'] = i
            
            # 为LR关系对添加配体来源
            # Add ligand source for LR relationship pairs
            temp_res = pd.merge(temp_LR_pair,temp_liga_source,on='ligands')
            
            # 存储得到的临时结果
            # Store the temporary result obtained
            res = pd.concat([res,temp_res])


    # 重建索引
    # Rebuild index     
    res.index = res['interaction_name']

    # 按LR分数从大到小对最终结果进行排序
    # Sort the final results by LR score from largest to smallest
    res.sort_values(by='LR_score',inplace=True,ascending=False)
    
    # 重新排列结果各列
    # Rearrange result columns
    res = res[['liga_source','liga_source_score','ligands','receptors',\
        'interaction_name','LR_score','celltype','Pvalue','FDR']]


    return res


