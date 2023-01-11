import warnings
import loadDB
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import fdrcorrection

warnings.filterwarnings("ignore")


"""
根据实验组和对照组配体蛋白丰度和受体细胞基因表达的差异，
推测细胞外信号与细胞之间的对话关系。

According to the differences in ligand protein abundance and 
receptor cell gene expression between experimental and control groups, 
the crosstalk between extracellular signals and cells was inferred.
"""


# 加载相关数据库文件
# Load the relevant database files
LR_sing, LR_comp, LR_liga, LR_rece = loadDB.loadLRInfo()


def calculateProtLog2FC(data,
                        info,
                        control,
                        exper):
    """
    计算两个分组中蛋白丰度的log2FC值。

    data: 配体蛋白丰度矩阵,
    info: 蛋白样本分组信息,
    control: 对照组名称,
    exper: 实验组名称。

    Calculate the log2FC value of protein abundance in both groups.

    data: ligand protein abundance matrix,
    info: sample grouping information,
    control: control group name,
    exper: experimental group name.
    """

    # 构建字典 {样本名：组别}
    # Build dict {sample: group}
    info_index = info.T.to_dict('records')[0]
    
    # 计算蛋白丰度均值
    # Calculate the mean protein abundance
    data_mean = data.groupby(info_index,axis=1).agg(np.mean)

    # 计算两个分组中蛋白丰度的log2FC值
    # Calculate the log2FC value of protein abundance in both groups
    data_mean['log2FC'] = data_mean[exper] - data_mean[control]
    

    return data_mean


def calculateGeneLog2FC(data,
                        info,
                        control,
                        exper,
                        expr_cutoff):
    """
    计算不同类型细胞中两个分组基因表达的log2FC值。

    data: 受体细胞基因表达矩阵, 
    info: 受体细胞类型分组注释, 
    control: 对照组名称, 
    exper: 实验组名称, 
    expr_cutoff: 0.25, 
        基因在该类型细胞中表达数量大于整体细胞数量*expr_cutoff时，
        才认为该基因在该类型细胞中表达，
        默认数值为 0.25。

    Calculate the log2FC value of gene expression in both groups of 
    different cell types.

    data: receptor cell gene expression matrix,
    info: receptor cell type grouping information,
    control: control group name,
    exper: experimental group name,
    expr_cutoff: 0.25,
        The gene is considered to be expressed in cell types only 
        when the cell numbers of gene expressed is greater than 
        the overall number of cells*expr_cutoff.
        The default value is 0.25.
    """

    # 创建储存结果的标准格式数据框
    # Create a standard format dataframe to store the results
    data_mean = pd.DataFrame(columns=[control,exper,'log2FC','celltype'])

    # 对不同类型细胞进行迭代
    # Circulate different types of cells
    for i in set(info.iloc[:,1].values):

        # 提取该类型细胞注释信息
        # Extract annotation information of the cell type
        cell_info = info[info.iloc[:,1]==i]
        cell_num = len(cell_info)
        
        # 该类型细胞构建字典 {细胞名：组别}
        # Build dict of the cell type {cellname: group}
        cell_index = cell_info.T.to_dict('records')[0]
        
        # 根据expr_cutoff来筛选细胞基因表达矩阵
        # Screen receptor cell gene expression matrix according to 
        # the expr_cutoff value
        cell_filter_expr = data[cell_info.index][data[cell_info.index].\
            astype(bool).sum(axis=1) > (cell_num * expr_cutoff)]
        cell_filter_expr_mean = cell_filter_expr.groupby(cell_index,axis=1).\
            agg(np.mean)
        cell_filter_expr_mean['log2FC'] = ''
        cell_filter_expr_mean['celltype'] = i
        
        data_mean = pd.concat([data_mean,cell_filter_expr_mean])
    
    # 为了避免基因表达量为0时计算过程中所造成的误差，使用data_mean中的非0最小值填充0值
    # To avoid the error caused in the calculation process 
    # when the gene expression level is 0, 
    # use the non-zero minimum value in data_mean to fill the zero value.
    temp_min = data_mean[[control,exper]].\
        apply(lambda x: np.array(x)[np.array(x).nonzero()].min()).min()
    data_mean[[control,exper]] = data_mean[[control,exper]].replace(0,temp_min)

    # 计算两个分组中基因表达的log2FC值
    # Calculate the log2FC value of gene expression in both groups
    data_mean['log2FC'] = np.log2(data_mean[exper] / data_mean[control])
    

    return data_mean


def getLigaMatrix(data,
                  liga=LR_liga):
    """
    根据获得的配体蛋白log2FC矩阵及需要输入的配体蛋白列表，筛选配体蛋白log2FC矩阵。

    data: 蛋白丰度log2FC矩阵, 
    liga: LR_liga, 
        输入的配体蛋白列表，默认从本算法自定义配体-受体关系对中加载。

    Screen the ligand protein log2FC matrix according to the obtained ligand 
    protein log2FC matrix and the ligand protein list that needs to be input.

    data: protein abundance log2FC matrix,
    liga: LR_liga, 
        The list of input ligand proteins, 
        which are loaded from the custom ligand-receptor relationship pairs 
        of the algorithm by default.
    """

    # 获取配体蛋白的交集
    # Get the intersection of ligand proteins
    data_liga = data.loc[data.index.intersection(set(liga))]
    

    return data_liga


def getReceMatrix(data,
                  org,
                  rece_type,
                  rece=LR_rece):
    """
    根据获得的受体细胞基因表达log2FC矩阵及需要输入的受体蛋白基因列表，
    筛选受体细胞基因表达log2FC矩阵。

    data: 细胞基因表达矩阵, 
    org: ['Human],
        研究的物种类型, 
        可选物种类型 Human, Mouse, Pig, 默认为 Human, 
    rece_type: ['Transmembrane','Cytoplasma'], 
        可选受体类型: 胞膜还是胞质, 默认为胞膜和胞质, 
    rece: LR_rece,
        输入的受体蛋白列表，默认从本算法自定义配体-受体关系对中加载。
    
    Screen the gene expression log2FC matrix of recipient cells 
    according to the obtained gene expression log2FC matrix of recipient cells 
    and the receptor protein list that needs to be input. 

    data: cell gene expression matrix,
    org: ['Human],
        The type of species studied,
        Optional species type: Human, Mouse, Pig.
        Default is Human.
    rece_type: ['Transmembrane','Cytoplasma'],
        Optional receptor type: transmembrane or cytoplasm, 
        default is transmembrane and cytoplasm.
    rece: LR_rece,
        The list of input receptor proteins, 
        which are loaded from the custom ligand-receptor relationship pairs 
        of the algorithm by default.    
    """

    # 加载相关数据库文件
    # Load the relevant database files
    transmembrane_prot, cytoplasma_prot = loadDB.loadProtAnno(org)

    # 判断选择的受体类型
    # Determine the receptor type to choose
    if rece_type == 'Transmembrane':
        rece_filter = set(transmembrane_prot['Gene names  (primary )']).\
            intersection(rece)
    elif rece_type == 'Cytoplasma':
        rece_filter = set(cytoplasma_prot['Gene names  (primary )']).\
            intersection(rece)
    else:
        rece_filter = set(list(transmembrane_prot['Gene names  (primary )']) + \
            list(cytoplasma_prot['Gene names  (primary )'])).intersection(rece)
    
    # 获取受体蛋白的交集
    # Get the intersection of receptor proteins
    data_rece = data.loc[data.index.intersection(rece_filter)]
    

    return data_rece


def calculateLRSingScore(liga_matrix,
                         rece_matrix,
                         liga_rece_inter=LR_sing):
    """
    根据自定义配体-受体关系对中的Single结构对, 
    计算配体蛋白与受体细胞之间的相互作用分数。

    liga_matrix: 筛选过后的配体蛋白log2FC矩阵, 
    rece_matrix: 筛选过后的受体细胞log2FC矩阵, 
    liga_rece_inter: 配体-受体关系对,
        本函数中默认为自定义配体-受体关系对中的Single结构对。
    
    Calculate the interaction score between ligand proteins and receptor cells 
    based on the 'Single' pairs in the custom ligand-receptor relationship pairs.

    liga_matrix: the screened ligand protein log2FC matrix,
    rece_matrix: the screened receptor cell log2FC matrix,
    liga_rece_inter: ligand-receptor relationship pair,
        Default value is the 'Single' pairs 
        in the custom ligand-receptor relationship pair.
    """

    # 1.根据liga_matrix中的配体筛选liga_rece_inter
    # 1.Screen liga_rece_inter against ligands in the liga_matrix
    liga_rece_inter_filter1 = liga_rece_inter[liga_rece_inter['ligands'].\
        isin(liga_matrix.index)]
    

    # 2.根据rece_matrix中的受体继续筛选liga_rece_inter
    # 2.Continue to screen liga_rece_inter based on receptors in the rece_matrix
    liga_rece_inter_filter2 = \
        liga_rece_inter_filter1[liga_rece_inter_filter1['receptors'].\
        isin(rece_matrix.index)]
    

    # 3.创建储存结果的标准格式数据框
    # 3.Create a standard format dataframe to store the results
    res_sing = pd.DataFrame(columns=\
        ['ligands','receptors','interaction_name','LR_score','celltype'])
    

    # 4.由于配体蛋白没有细胞类型的差异，构建字典{配体蛋白：log2FC}
    # 4.Build dict {ligand protein: log2FC}
    liga_index = {}
    for row in liga_matrix.itertuples():
        liga_index[row[0]] = row[3]
    

    # 5.对筛选得到的 liga_rece_inter 循环逐行计算LR分数
    # 5.Calculate the LR score line by line for the filtered liga_rece_inter
    for row in liga_rece_inter_filter2.itertuples():

        # 创建临时储存结果的数据框
        # Create a data frame to temporarily store the results
        temp_res = pd.DataFrame(columns=\
            ['ligands','receptors','interaction_name','LR_score','celltype'])
        
        # 判断受体是否存在于多个细胞类型中，
        # 只有一种细胞时需将自动转换的 pd.Series 类型数据转换为 pd.DataFrame
        # To determine whether the receptor exists in multiple cell types, 
        # convert the automatically converted pd.Series type data 
        # to pd.DataFrame when there is only one cell type.
        temp_rece = rece_matrix.loc[row[2]]
        if type(temp_rece) == pd.Series:
            temp_rece = pd.DataFrame([temp_rece.to_dict()])
        
        # 计算LR score, 等于配体log2FC * 受体log2FC
        # Calculate the LR score, 
        # which is equal to ligand log2FC * receptor log2FC
        temp_score = liga_index[row[1]] * temp_rece['log2FC']
        temp_res['LR_score'] = temp_score
        temp_res['ligands'] = row[1]
        temp_res['receptors'] = row[2]
        temp_res['interaction_name'] = '-'.join([row[1],row[2]])
        temp_res['celltype'] = temp_rece.iloc[:,3]
        
        # 存储得到的临时结果
        # Store the temporary result obtained
        res_sing = pd.concat([res_sing,temp_res])
        

    # 6.筛选计算得到的LR score，当LR score > 0时代表配体和受体的变化趋势一致
    # 6.Screen the calculated LR score, when LR score > 0, 
    #   it means that the changing trend of ligand and receptor is consistent.
    res_sing = res_sing[res_sing['LR_score']>0]


    # 7.重建索引
    # 7.Rebuild index
    res_sing.index = res_sing['interaction_name']


    return res_sing


def calculateLRCompScore(liga_matrix,
                         rece_matrix,
                         comp_expr_ratio,
                         liga_rece_inter=LR_comp):
    """
    根据自定义配体-受体关系对中的Complex结构对，
    计算配体蛋白与受体细胞之间的相互作用分数。

    liga_matrix: 筛选过后的配体蛋白log2FC矩阵, 
    rece_matrix: 筛选过后的受体细胞log2FC矩阵, 
    compl_expr_ratio: 0.5, 
        complex结构中基因需要变化的比例, 
        complex结构中需要受体和配体中变化的基因分别超过total genes*complex_expr_ratio, 
        默认为0.5。
    liga_rece_inter: 配体-受体关系对,
        本函数中默认为自定义配体-受体关系对中的Complex结构对。

    Calculate the interaction score between ligand proteins and receptor cells 
    based on the 'Complex' pairs in the custom ligand-receptor relationship pairs.

    liga_matrix: the screened ligand protein log2FC matrix,
    rece_matrix: the screened receptor cell log2FC matrix,
    compl_expr_ratio: 0.5, 
        The ratio of gene that need to be changed in the 'Complex' pairs.
        The genes in the 'Complex' pairs that require changes in receptors and 
        ligands exceed total genes*complex_expr_ratio, respectively.
        The default value is 0.5.
    liga_rece_inter: ligand-receptor relationship pair,
        Default value is the 'Complex' pairs 
        in the custom ligand-receptor relationship pair.
    """
    
    # 1.根据liga_matrix中的配体筛选liga_rece_inter
    # 1.Screen liga_rece_inter against ligands in the liga_matrix

    # 判断liga_rece_inter中的配体是否在liga_matrix中, 
    # 且是否大于配体蛋白总数*complex_expr_ratio
    # Determine whether the ligand in liga_rece_inter is in liga_matrix and 
    # the ligand protein number is greater than 
    # the total number of ligand proteins*complex_expr_ratio.
    liga_rece_inter_liga_exist = []
    for i in range(len(liga_rece_inter)):
        temp1 = [j in liga_matrix.index \
            for j in liga_rece_inter.loc[i,'ligands'].split('_')]
        if temp1.count(True) > (len(temp1) * comp_expr_ratio):
            liga_rece_inter_liga_exist.append(True)
        else:
            liga_rece_inter_liga_exist.append(False)
    
    # 筛选liga_rece_inter
    # Screen liga_rece_inter
    liga_rece_inter_filter1 = liga_rece_inter[liga_rece_inter_liga_exist]

    # 重建索引
    # Rebuild index
    liga_rece_inter_filter1 = liga_rece_inter_filter1.reset_index(drop=True)
    
    
    # 2.根据rece_matrix中的受体继续筛选liga_rece_inter
    # 2.Continue to screen liga_rece_inter based on receptors in the rece_matrix

    # 判断liga_rece_inter中的受体是否在rece_matrix中, 
    # 且是否大于受体蛋白总数*complex_expr_ratio
    # Determine whether the receptor in liga_rece_inter is in rece_matrix and 
    # the receptor protein number is greater than 
    # the total number of receptor proteins*complex_expr_ratio.
    liga_rece_inter_rece_exist = []
    for i in range(len(liga_rece_inter_filter1)):
        temp2 = [j in rece_matrix.index \
            for j in liga_rece_inter_filter1.loc[i,'receptors'].split('_')]
        if temp2.count(True) > (len(temp2) * comp_expr_ratio):
            liga_rece_inter_rece_exist.append(True)
        else:
            liga_rece_inter_rece_exist.append(False)
            
    # 筛选liga_rece_inter
    # Screen liga_rece_inter
    liga_rece_inter_filter2 = liga_rece_inter_filter1[liga_rece_inter_rece_exist]

    # 重建索引
    # Rebuild index
    liga_rece_inter_filter2 = liga_rece_inter_filter2.reset_index(drop=True)
    

    # 3.创建储存结果的标准格式数据框
    # 3.Create a standard format dataframe to store the results
    res_comp = pd.DataFrame(columns=\
        ['ligands','receptors','interaction_name','LR_score','celltype'])
    

    # 4.由于配体蛋白没有细胞类型的差异，构建字典{配体蛋白：log2FC}
    # 4.Build dict {ligand protein: log2FC}

    # 构建字典{单一配体蛋白：log2FC}
    # Build dict {single ligand protein: log2FC}
    liga_index = {}
    for row in liga_matrix.itertuples():
        liga_index[row[0]] = row[3]
    

    def _is_in_liga(x):
        """
        判断Complex结构中的配体是否真的在配体中。

        Determine whether the ligand in the Complex pairs 
        is really in the ligand protein log2FC matrix.
        """

        return x in liga_matrix.index
    
    # 构建字典{Complex结构配体蛋白：log2FC}
    # Build dict {Complex ligand protein: log2FC}
    liga_comp_index = {}
    # Complex结构中配体的log2FC取Complex结构对中所有表达配体的log2FC值之和
    # The log2FC value of the ligand in the Complex pairs 
    # is the sum of log2FC values of all expressed ligands.
    for row in liga_rece_inter_filter2.itertuples():
        temp3 = list(filter(_is_in_liga, row[1].split('_')))
        temp_comp_liga = 0
        for i in temp3:
            temp_comp_liga += liga_index[i]
        liga_comp_index[row[1]] = temp_comp_liga


    # 5.对不同细胞类型迭代
    # 5.Iterate over different cell types
    for celltype in set(rece_matrix['celltype']):
    
        def _is_in_rece(x):
            """
            判断Complex结构中的受体是否真的在该类型细胞的受体中。

            Determine whether the receptor in the Complex pairs 
            is really in the receptor cell log2FC matrix of the cell type.
            """

            return x in rece_matrix[rece_matrix['celltype']==celltype].index
        
        # 对筛选得到的 liga_rece_inter 循环逐行计算LR分数
        # Calculate the LR score line by line for the filtered liga_rece_inter
        for row in liga_rece_inter_filter2.itertuples():
            temp4 = list(filter(_is_in_rece, row[2].split('_')))

            # 判断受体在该细胞类型中是否为空
            # Determine if the receptor is null in the cell type
            if (len(temp4)==0):
                next
            else:
                temp_comp_rece = 0
                # Complex结构中受体的log2FC取Complex结构对中所有表达受体的log2FC值之和
                # The log2FC value of the receptor in the Complex pairs 
                # is the sum of log2FC values of all expressed receptors.
                for i in temp4:
                    temp_comp_rece += \
                        rece_matrix[rece_matrix['celltype']==celltype].\
                        loc[i,'log2FC']
                
                # 计算LR score, 等于配体log2FC * 受体log2FC
                # Calculate the LR score, 
                # which is equal to ligand log2FC * receptor log2FC
                temp_score = liga_comp_index[row[1]] * temp_comp_rece
                
                # 创建临时储存结果的数据框
                # Create a data frame to temporarily store the results
                temp_res = pd.DataFrame({'ligands':row[1],
                                         'receptors':row[2],
                                         'interaction_name':\
                                             '-'.join([row[1],row[2]]),
                                         'LR_score':temp_score,
                                         'celltype':celltype},
                                         index=[0])

                # 存储得到的临时结果
                # Store the temporary result obtained
                res_comp = pd.concat([res_comp,temp_res])
        

    # 6.筛选计算得到的LR score，当LR score > 0时代表配体和受体的变化趋势一致
    # 6.Screen the calculated LR score, when LR score > 0, 
    #   it means that the changing trend of ligand and receptor is consistent.
    res_comp = res_comp[res_comp['LR_score']>0]


    # 7.重建索引
    # 7.Rebuild index
    res_comp.index = res_comp['interaction_name']
    

    return res_comp


def _calculateVirtualLRSingScore(rece_data,
                                 rece_info,
                                 rece_control,
                                 rece_exper,
                                 celltype,
                                 CE_index,
                                 liga_matrix,
                                 rece_list):
    """
    根据细胞类型随机挑选相应分组数量的细胞, 并根据该类型细胞受体筛选相应受体基因表达谱,
    计算virtualLRSingScore。

    rece_data: 受体细胞基因表达矩阵,
    rece_info: 受体细胞类型分组注释,
    rece_control: 受体细胞对照组名称,
    rece_exper: 受体细胞实验组名称, 
    celltype: 需要处理的细胞类型,
    CE_index: 字典{受体细胞: control/exper分组信息},
    liga_matrix: 该细胞类型下所有配体的log2FC值,
    rece_list: 该细胞类型下的所有受体列表。

    The corresponding number cells are randomly selected according to 
    the celltype, and screening the receptor gene expression profiles, 
    followed by calculating the virtualLRSingScore.

    rece_data: receptor cell gene expression matrix,
    rece_info: receptor cell type grouping information,
    rece_control: receptor cell control group name,
    rece_exper: receptor cell experimental group name,
    celltype: the cell type to be processed,
    CE_index: dictionary {receptor cell: control/exper grouping information},
    liga_matrix: log2FC value of all ligands under the cell type,
    rece_list: list of all receptors under the cell type.
    """

    # 该类型细胞总数
    # The total number of cells of the celltype
    temp_cell_num = rece_info[rece_info.iloc[:,1]==celltype].\
        groupby(list(rece_info)[0]).count()
    
    # 对照组细胞数
    # The number of cells in the control group
    temp_control_num = temp_cell_num.loc[rece_control].values

    # 实验组细胞数
    # The number of cells in the experimental group
    temp_exper_num = temp_cell_num.loc[rece_exper].values

    # 随机挑选对照组细胞
    # Randomly selected control group cells
    random_control_cell = rece_info[rece_info.iloc[:,0]==rece_control].\
        sample(temp_control_num)
    
    # 随机挑选实验组细胞
    # Randomly select experimental group cells
    random_exper_cell = rece_info[rece_info.iloc[:,0]==rece_exper].\
        sample(temp_exper_num)
    
    # 整合挑选的两组细胞
    # Integrate the selected two groups of cells
    random_CE_cell = rece_data[pd.concat([random_control_cell,\
        random_exper_cell]).index]


    # 计算两个分组中受体基因表达的log2FC值
    # Calculate the log2FC value of receptor gene expression in both groups
    random_CE_cell_mean = random_CE_cell.loc[rece_list].\
        groupby(CE_index,axis=1).agg(np.mean)
    random_CE_cell_mean['log2FC'] = np.log2(random_CE_cell_mean[rece_exper] / \
        random_CE_cell_mean[rece_control])

    # 计算虚拟LR score, 等于配体log2FC * 虚拟受体log2FC
    # Calculate the virtual LR score, 
    # which is equal to ligand log2FC * virtual receptor log2FC
    virtual_LRSingScore = liga_matrix.values * \
        random_CE_cell_mean['log2FC'].values


    return virtual_LRSingScore


def calculateLRSingPvalue(LR_sing_data,
                          liga_matrix,
                          rece_data,
                          rece_info,
                          rece_control,
                          rece_exper,
                          perm_num,
                          random_seed):
    """
    利用置换检验计算获得的Single配体-受体关系对分数的P值。

    LR_sing_data: calculateLRSingScore函数计算得到的Single配体-受体关系对评分, 
    liga_matrix: calculateProtLog2FC函数计算得到的配体蛋白log2FC矩阵, 
    rece_data: 受体细胞基因表达矩阵, 
    rece_info: 受体细胞类型分组注释，
    rece_control: 受体细胞对照组名称,
    rece_exper: 受体细胞实验组名称, 
    perm_num: 1000,
        置换检验次数, 默认为1000, 
    random_seed: 123, 
        随机种子, 默认为123。

    P-values for scores obtained for 'Single' ligand-receptor relationship pairs 
    were calculated using permutation tests.

    LR_sing_data: the 'Single' ligand-receptor relationship pair score
        calculated by the calculateLRSingScore function,
    liga_matrix: the ligand protein log2FC matrix 
        calculated by the calculateProtLog2FC function,
    rece_data: receptor cell gene expression matrix,
    rece_info: receptor cell type grouping information,
    rece_control: receptor cell control group name,
    rece_exper: receptor cell experimental group name,
    perm_num: 1000,
        number of permutation tests, 
        The default value is 1000.
    random_seed: 123,
        random seed, 
        The default value is 123.
    """

    # 1.构建字典{配体蛋白：log2FC}
    # 1.Build dict {ligand protein: log2FC}
    liga_index = liga_matrix['log2FC'].to_dict()


    # 2.构建字典{受体细胞：control/exper分组信息}，简称CE_index
    # 2.Build dict {receptor cells: control/exper grouping information}, 
    #   referred to as CE_index
    CE_index = rece_info.T.to_dict('records')[0]


    # 3.创建储存结果的标准格式数据框
    # 3.Create a standard format dataframe to store the results
    LR_sing_Pvalue = \
        pd.DataFrame(columns=\
        ['ligands','receptors','interaction_name','LR_score','celltype','Pvalue'])


    # 4.对不同细胞类型迭代
    # 4.Iterate over different cell types
    for celltype in set(rece_info.iloc[:,1]):

        # 提取该细胞类型下的Single关系对
        # Extract the 'Single' relationship pairs under the cell type
        temp_cell = LR_sing_data[LR_sing_data['celltype'] == celltype]

        # 获得该细胞类型下的真实LR score
        # Obtain the true LR score for the cell type
        temp_real_LRScore = temp_cell['LR_score']

        # 获得该细胞类型下的所有配体的log2FC值
        # Obtain log2FC values for all ligands under the cell type
        temp_liga_FC = temp_cell['ligands'].map(liga_index)

        # 获得该细胞类型下的所有受体
        # Obtain all receptors under the cell type
        temp_rece = temp_cell['receptors']
        
        # 创建临时储存结果的数据框
        # Create a data frame to temporarily store the results
        virtual_LRScore = pd.DataFrame(index=temp_cell.index)

        # 设置随机种子
        # Set random seed
        np.random.seed(random_seed)

        # 置换迭代
        # Permutation iteration
        for i in range(int(perm_num-1)):
            temp_virtual_LRScore = \
                pd.DataFrame({'perm_{}'.format(i): \
                _calculateVirtualLRSingScore(rece_data,rece_info,rece_control,\
                rece_exper,celltype,CE_index,temp_liga_FC,temp_rece)}, \
                index=virtual_LRScore.index)
            virtual_LRScore = pd.concat([virtual_LRScore,temp_virtual_LRScore],\
                axis=1)

        # 将真实LR score添加至虚拟LR Score数据框中的最后一列
        # Add real LR score to the last column of the virtual LR Score dataframe
        virtual_LRScore = pd.concat([virtual_LRScore,\
            pd.DataFrame({'real_LRScore':temp_real_LRScore.values}, \
            index=virtual_LRScore.index)], axis=1)

        # 计算P值
        # Calculate the P value
        temp_Pvalue = []
        for row in virtual_LRScore.itertuples():
            temp_Pvalue.\
                append((sum([i > row[-1] for i in row[1:-1]])+1)/(perm_num))
        
        temp_cell = pd.concat([temp_cell,pd.DataFrame({'Pvalue':temp_Pvalue},\
            index=temp_cell.index)],axis=1)
        
        # 存储得到的临时结果
        # Store the temporary result obtained
        LR_sing_Pvalue = pd.concat([LR_sing_Pvalue,temp_cell])
    
            
    return LR_sing_Pvalue


def _calculateVirtualLRCompScore(rece_data,
                                 rece_info,
                                 rece_control,
                                 rece_exper,
                                 celltype,
                                 CE_index,
                                 liga_matrix,
                                 rece_list,
                                 rece_filter_matrix):
    """
    根据细胞类型随机挑选相应分组数量的细胞, 并根据该类型细胞受体筛选相应受体基因表达谱,
    计算virtualLRSingScore, 考虑Complex结构中的配体-受体基因表达情况。

    rece_data: 受体细胞基因表达矩阵,
    rece_info: 受体细胞类型分组注释,
    rece_control: 受体细胞对照组名称,
    rece_exper: 受体细胞实验组名称, 
    celltype: 需要处理的细胞类型,
    CE_index: 字典{受体细胞: control/exper分组信息},
    liga_matrix: 该细胞类型下所有配体的log2FC值,
    rece_list: 该细胞类型下的所有受体列表, 
    rece_filter_matrix: 经getReceMatrix函数筛选后的受体基因log2FC矩阵。

    The corresponding number cells are randomly selected according to 
    the celltype, and screening the receptor gene expression profiles, 
    followed by calculating the virtualLRSingScore. 
    Consider the ligand-receptor gene expression in the 'Complex' pairs.

    rece_data: receptor cell gene expression matrix,
    rece_info: receptor cell type grouping information,
    rece_control: receptor cell control group name,
    rece_exper: receptor cell experimental group name,
    celltype: the cell type to be processed,
    CE_index: dictionary {receptor cell: control/exper grouping information},
    liga_matrix: log2FC value of all ligands under the cell type,
    rece_list: list of all receptors under the cell type,
    rece_filter_matrix: the log2FC matrix of receptor genes 
        filtered by the getReceMatrix function.
    """

    # 该类型细胞总数
    # The total number of cells of the celltype
    temp_cell_num = rece_info[rece_info.iloc[:,1]==celltype].\
        groupby(list(rece_info)[0]).count()
    
    # 对照组细胞数
    # The number of cells in the control group
    temp_control_num = temp_cell_num.loc[rece_control].values

    # 实验组细胞数
    # The number of cells in the experimental group
    temp_exper_num = temp_cell_num.loc[rece_exper].values

    # 随机挑选对照组细胞
    # Randomly selected control group cells
    random_control_cell = rece_info[rece_info.iloc[:,0]==rece_control].\
        sample(temp_control_num)
    
    # 随机挑选实验组细胞
    # Randomly select experimental group cells
    random_exper_cell = rece_info[rece_info.iloc[:,0]==rece_exper].\
        sample(temp_exper_num)
    
    # 整合挑选的两组细胞
    # Integrate the selected two groups of cells
    random_CE_cell = rece_data[pd.concat([random_control_cell,\
        random_exper_cell]).index]
    

    def _is_in_rece(x):
        """
        判断Complex结构中的受体是否真的在该类型细胞的受体中。

        Determine whether the receptor in the Complex pairs 
        is really in the receptor cell log2FC matrix of the cell type.
        """

        temp = x in rece_filter_matrix[rece_filter_matrix['celltype']==\
            celltype].index

        return temp
    
    # 计算Complex结构中每个受体基因在两个分组中的log2FC值之和
    # Calculate the sum of log2FC values for each receptor gene 
    # in the 'Complex' pairs in both groups
    random_CE_cell_FC = []
    for j in [list(filter(_is_in_rece, i.split('_'))) for i in rece_list]:
        temp_rece_FC = 0
        for temp_rece in j:
            temp_rece_mean = \
                random_CE_cell.loc[temp_rece].groupby(CE_index).agg(np.mean)
            temp_rece_FC += np.log2(temp_rece_mean[rece_exper] / \
                temp_rece_mean[rece_control])
        random_CE_cell_FC.append(temp_rece_FC)

    # 计算虚拟LR score, 等于配体log2FC * 虚拟受体log2FC
    # Calculate the virtual LR score, 
    # which is equal to ligand log2FC * virtual receptor log2FC
    virtual_LRCompScore = liga_matrix.values * pd.Series(random_CE_cell_FC)


    return virtual_LRCompScore


def calculateLRCompPvalue(LR_comp_data,
                          liga_matrix,
                          rece_data,
                          rece_info,
                          rece_filter_matrix,
                          rece_control,
                          rece_exper,
                          perm_num,
                          random_seed):
    """
    利用置换检验计算获得的Complex配体-受体关系对分数的P值。

    LR_sing_data: calculateLRSingScore函数计算得到的Complex配体-受体关系对评分, 
    liga_matrix: calculateProtLog2FC函数计算得到的配体蛋白log2FC矩阵, 
    rece_data: 受体细胞基因表达矩阵, 
    rece_info: 受体细胞类型分组注释, 
    rece_filter_matrix: 经getReceMatrix函数筛选后的受体基因log2FC矩阵, 
    rece_control: 受体细胞对照组名称,
    rece_exper: 受体细胞实验组名称, 
    perm_num: 1000,
        置换检验次数, 默认为1000, 
    random_seed: 123, 
        随机种子, 默认为123。

    P-values for scores obtained for 'Complex' ligand-receptor relationship pairs 
    were calculated using permutation tests.

    LR_sing_data: the 'Complex' ligand-receptor relationship pair score
        calculated by the calculateLRSingScore function,
    liga_matrix: the ligand protein log2FC matrix 
        calculated by the calculateProtLog2FC function,
    rece_data: receptor cell gene expression matrix,
    rece_info: receptor cell type grouping information,
    rece_filter_matrix: the log2FC matrix of receptor genes 
        filtered by the getReceMatrix function,
    rece_control: receptor cell control group name,
    rece_exper: receptor cell experimental group name,
    perm_num: 1000,
        number of permutation tests, 
        The default value is 1000.
    random_seed: 123,
        random seed, 
        The default value is 123.
    """

    # 判断Complex结构关系对结果是否为空
    # Determine whether the result of the 'Complex' relationship pair is empty
    if (len(LR_comp_data)==0):
        # 结果为空时输出空数据框
        # Output empty dataframe when result is empty
        LR_comp_Pvalue = \
            pd.DataFrame(
            columns=['ligands','receptors','interaction_name','LR_score',\
            'celltype','Pvalue'])
    
    else:
        # 1.由于配体蛋白没有细胞类型的差异，构建字典{配体蛋白：log2FC}
        # 1.Build dict {ligand protein: log2FC}

        # 构建字典{单一配体蛋白：log2FC}
        # Build dict {single ligand protein: log2FC}
        liga_index = {}
        for row in liga_matrix.itertuples():
            liga_index[row[0]] = row[3]
        

        def _is_in_liga(x):
            """
            判断Complex结构中的配体是否真的在配体中。

            Determine whether the ligand in the Complex pairs 
            is really in the ligand protein log2FC matrix.
            """

            return x in liga_matrix.index
        
        # 构建字典{Complex结构配体蛋白：log2FC}
        # Build dict {Complex ligand protein: log2FC}
        data1_comp_index = {}
        # Complex结构中配体的log2FC取Complex结构对中所有表达配体的log2FC值之和
        # The log2FC value of the ligand in the Complex pairs 
        # is the sum of log2FC values of all expressed ligands.
        for row in LR_comp_data.itertuples():
            temp = list(filter(_is_in_liga, row[1].split('_')))
            temp_comp_liga = 0
            for i in temp:
                temp_comp_liga += liga_index[i]
            data1_comp_index[row[1]] = temp_comp_liga
        

        # 2.构建字典{受体细胞：control/exper分组信息}，简称CE_index
        # 2.Build dict {receptor cells: control/exper grouping information}, 
        #   referred to as CE_index
        CE_index = rece_info.T.to_dict('records')[0]
        

        # 3.创建储存结果的标准格式数据框
        # 3.Create a standard format dataframe to store the results
        LR_comp_Pvalue = \
            pd.DataFrame(
            columns=['ligands','receptors','interaction_name','LR_score',\
            'celltype','Pvalue'])
        

        # 4.对不同细胞类型迭代
        # 4.Iterate over different cell types
        for celltype in set(rece_info.iloc[:,1]):

            # 提取该细胞类型下的Complex关系对
            # Extract the 'Complex' relationship pairs under the cell type
            temp_cell = LR_comp_data[LR_comp_data['celltype'] == celltype]

            # 获得该细胞类型下的真实LR score
            # Obtain the true LR score for the cell type
            temp_real_LRScore = temp_cell['LR_score']

            # 获得该细胞类型下的所有配体的log2FC值
            # Obtain log2FC values for all ligands under the cell type
            temp_liga_FC = temp_cell['ligands'].map(data1_comp_index)

            # 获得该细胞类型下的所有受体
            # Obtain all receptors under the cell type
            temp_rece = temp_cell['receptors']
            
            # 创建临时储存结果的数据框
            # Create a data frame to temporarily store the results
            virtual_LRScore = pd.DataFrame(index=temp_cell.index)

            # 设置随机种子
            # Set random seed
            np.random.seed(random_seed)

            # 置换迭代
            # Permutation iteration
            for i in range(int(perm_num-1)):
                temp_virtual_LRScore = pd.DataFrame({'perm_{}'.format(i): \
                    _calculateVirtualLRCompScore(rece_data,rece_info,\
                    rece_control,rece_exper,celltype,CE_index,temp_liga_FC,\
                    temp_rece,rece_filter_matrix)},index=virtual_LRScore.index)
                virtual_LRScore = pd.concat([virtual_LRScore,\
                    temp_virtual_LRScore],axis=1)

            # 将真实LR score添加至虚拟LR Score数据框中的最后一列
            # Add real LR score to the last column of 
            # the virtual LR Score dataframe
            virtual_LRScore = pd.concat([virtual_LRScore,\
                pd.DataFrame({'real_LRScore':temp_real_LRScore.values},\
                index=virtual_LRScore.index)],axis=1)

            # 计算P值
            # Calculate the P value
            temp_Pvalue = []
            for row in virtual_LRScore.itertuples():
                temp_Pvalue.\
                    append((sum([i > row[-1] for i in row[1:-1]])+1)/(perm_num))
            
            temp_cell = pd.concat([temp_cell,pd.DataFrame({'Pvalue':temp_Pvalue},\
                index=temp_cell.index)],axis=1)
            
            # 存储得到的临时结果
            # Store the temporary result obtained
            LR_comp_Pvalue = pd.concat([LR_comp_Pvalue,temp_cell])
    
            
    return LR_comp_Pvalue


def integrateLRSingCompRes(LRSingRes,
                           LRCompRes):
    """
    整合计算得到的Single和Complex配体-受体关系对结果。

    LRSingRes: Single结构配体-受体关系对结果, 
    LRCompRes: Complex结构配体-受体关系对结果。

    Integrate calculated 'Single' and 'Complex' 
    ligand-receptor relationship pair results.

    LRSingRes: 'Single' ligand-receptor relationship pair results,
    LRCompRes: 'Complex' ligand-receptor relationship pair results.
    """

    # 拼接Single和Complex配体-受体关系对结果
    # Integrate 'Single' and 'Complex' 
    # ligand-receptor relationship pair results.
    LRScore = pd.concat([LRSingRes,LRCompRes])

    # 按LR分数从大到小对最终结果进行排序
    # Sort the final results by LR score from largest to smallest
    LRScore.sort_values(by='LR_score',inplace=True,ascending=False)
    

    return LRScore


def correctPvalue(data):
    """
    利用statsmodels库中的fdrcorrection函数对置换检验获得的P值进行FDR校正。

    data: P值列表。

    P-values obtained from permutation tests were FDR-corrected 
    using the 'fdrcorrection' function in the 'statsmodels' pbrary.

    data: List of p-values.
    """

    data['FDR'] = fdrcorrection(data['Pvalue'])[1]
    

    return data


