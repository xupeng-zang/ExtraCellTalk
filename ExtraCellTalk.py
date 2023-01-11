import warnings
import calculateLRScore
import inferLigaSource
import plotRes

warnings.filterwarnings("ignore")


"""
主函数，程序执行。

Main function, program execution.
"""


def ExtraCellTalk(liga_data,
                  liga_info,
                  liga_control,
                  liga_exper,
                  rece_data,
                  rece_info,
                  rece_control,
                  rece_exper,
                  org=['Human'],
                  rece_type=['Transmembrane','Cytoplasma'],
                  rece_expr_cutoff=0.25,
                  comp_expr_ratio=0.5,
                  pvalue_cutoff=0.05,
                  FDR_cutoff=0.05,
                  perm_num=1000,
                  random_seed=123):
    """
    ExtraCellTalk程序的执行函数。

    liga_data: 配体蛋白丰度矩阵,
    liga_info: 配体蛋白样本分组信息,
    liga_control: 配体蛋白对照组名称,
    liga_exper: 配体蛋白实验组名称, 
    rece_data: 受体细胞基因表达矩阵, 
    rece_info: 受体细胞类型分组注释, 
    rece_control: 受体细胞对照组名称, 
    rece_exper: 受体细胞实验组名称, 
    org: ['Human],
        研究的物种类型, 
        可选物种类型 Human, Mouse, Pig, 默认为 Human, 
    rece_type: ['Transmembrane','Cytoplasma'], 
        可选受体类型: 胞膜还是胞质, 默认为胞膜和胞质,
    rece_expr_cutoff: 0.25, 
        基因在该类型细胞中表达数量大于整体细胞数量*expr_cutoff时，
        才认为该基因在该类型细胞中表达，
        默认数值为 0.25。
    compl_expr_ratio: 0.5, 
        自定义配体-受体关系对中complex结构中基因需要变化的比例, 
        complex结构中需要受体和配体中变化的基因分别超过total genes*complex_expr_ratio, 
        默认数值为0.5。
    pvalue_cutoff: 0.05,
        对LR关系对中未校正的P值过滤,
        默认数值为0.05。
    FDR_cutoff: 0.05,
        对LR关系对中校正后的P值FDR值过滤, FDR<0.05考虑LR关系显著,
        默认数值为0.05。
    perm_num: 1000, 
        置换检验的次数,
        默认数值为1000。
    random_seed: 123,
        随机种子设置,
        默认数值为123。

    The execution function of ExtraCellTalk.

    liga_data: ligand protein abundance matrix,
    liga_info: ligand protein sample grouping information,
    liga_control: ligand protein control group name,
    liga_exper: ligand protein experimental group name,
    rece_data: receptor cell gene expression matrix,
    rece_info: receptor cell type grouping information,
    rece_control: receptor cell control group name,
    rece_exper: receptor cell experimental group name,
    org: ['Human],
        The type of species studied,
        Optional species type: Human, Mouse, Pig.
        Default is Human.
    rece_type: ['Transmembrane','Cytoplasma'],
        Optional receptor type: transmembrane or cytoplasm. 
        Default is transmembrane and cytoplasm,
    rece_expr_cutoff: 0.25,
        The gene is considered to be expressed in cell types only 
        when the cell numbers of gene expressed is greater than 
        the overall number of cells*expr_cutoff.
        The default value is 0.25.
    compl_expr_ratio: 0.5, 
        The ratio of gene in the 'Complex' pairs that need to be changed 
        in the custom ligand-receptor relationship pairs.
        The genes in the 'Complex' pairs that require changes in receptors and 
        ligands exceed total genes*complex_expr_ratio, respectively.
        The default value is 0.5.
    pvalue_cutoff: 0.05,
        Filter the uncorrected P-value in LR relationship pairs.
        The default value is 0.05.
    FDR_cutoff: 0.05,
        Filter the corrected P-value in LR relationship pairs, that is, FDR value. 
        FDR < 0.05 considers the LR relationship pair to be significant.
        The default value is 0.05.
    perm_num: 1000,
        Number of permutation tests.
        The default value is 1000.
    random_seed: 123,
        Random seed setting.
        The default value is 123.
    """

    # 计算两个分组中配体蛋白丰度的log2FC值。
    # Calculate the log2FC value of ligand protein abundance in both groups.
    liga_mean = calculateLRScore.calculateProtLog2FC(liga_data,
                                                     liga_info,
                                                     liga_control,
                                                     liga_exper)

    # 计算不同类型细胞中两个分组受体基因表达的log2FC值。
    # Calculate the log2FC value of receptor gene expression in both groups of 
    # different cell types.
    rece_mean = calculateLRScore.calculateGeneLog2FC(rece_data,
                                                     rece_info,
                                                     rece_control,
                                                     rece_exper,
                                                     rece_expr_cutoff)


    # 根据获得的配体蛋白log2FC矩阵及需要输入的配体蛋白列表，筛选配体蛋白log2FC矩阵。
    # Screen the ligand protein log2FC matrix according to the obtained ligand 
    # protein log2FC matrix and the ligand protein list that needs to be input.
    liga_matrix = calculateLRScore.getLigaMatrix(liga_mean)

    # 根据获得的受体细胞基因表达log2FC矩阵及需要输入的受体蛋白基因列表，
    # 筛选受体细胞基因表达log2FC矩阵。
    # Screen the gene expression log2FC matrix of recipient cells 
    # according to the obtained gene expression log2FC matrix of recipient cells 
    # and the receptor protein list that needs to be input. 
    rece_matrix = calculateLRScore.getReceMatrix(rece_mean,
                                                 org,
                                                 rece_type)


    # 根据自定义配体-受体关系对中的Single结构对, 
    # 计算配体蛋白与受体细胞之间的相互作用分数。
    # Calculate the interaction score between ligand proteins and receptor cells 
    # based on the 'Single' pairs in the custom ligand-receptor relationship pairs.
    LR_sing_score = calculateLRScore.calculateLRSingScore(liga_matrix,
                                                          rece_matrix)

    # 利用置换检验计算获得的Single配体-受体关系对分数的P值。
    # P-values for scores obtained for 'Single' ligand-receptor relationship pairs 
    # were calculated using permutation tests.
    LR_sing_score_pvalue = calculateLRScore.calculateLRSingPvalue(LR_sing_score,
                                                                  liga_mean,
                                                                  rece_data,
                                                                  rece_info,
                                                                  rece_control,
                                                                  rece_exper,
                                                                  perm_num,
                                                                  random_seed)


    # 根据自定义配体-受体关系对中的Complex结构对，
    # 计算配体蛋白与受体细胞之间的相互作用分数。
    # Calculate the interaction score between ligand proteins and receptor cells 
    # based on the 'Complex' pairs in the custom ligand-receptor relationship pairs.
    LR_comp_score = calculateLRScore.calculateLRCompScore(liga_matrix,
                                                          rece_matrix,
                                                          comp_expr_ratio)
    
    # 利用置换检验计算获得的Complex配体-受体关系对分数的P值。
    # P-values for scores obtained for 'Complex' ligand-receptor relationship pairs 
    # were calculated using permutation tests.
    LR_comp_score_pvalue = calculateLRScore.calculateLRCompPvalue(LR_comp_score,
                                                                  liga_mean,
                                                                  rece_data,
                                                                  rece_info,
                                                                  rece_matrix,
                                                                  rece_control,
                                                                  rece_exper,
                                                                  perm_num,
                                                                  random_seed)


    # 整合计算得到的Single和Complex配体-受体关系对结果。
    # Integrate calculated 'Single' and 'Complex' 
    # ligand-receptor relationship pair results.
    LR_score_pvalue = \
        calculateLRScore.integrateLRSingCompRes(LR_sing_score_pvalue,
                                                LR_comp_score_pvalue)
    
    # 利用correctPvalue函数对置换检验获得的P值进行FDR校正。
    # P-values obtained from permutation tests were FDR-corrected 
    # using the 'correctPvalue' function.
    LR_score_pvalue_FDR = calculateLRScore.correctPvalue(LR_score_pvalue)


    # 对P值过滤
    # Filter on p-values
    LR_score_pvalue_FDR = \
        LR_score_pvalue_FDR[LR_score_pvalue_FDR['Pvalue']<pvalue_cutoff]

    # 对FDR值过滤
    # Filter on FDR-values
    LR_score_pvalue_FDR = \
        LR_score_pvalue_FDR[LR_score_pvalue_FDR['FDR']<FDR_cutoff]


    return LR_score_pvalue_FDR


def LigaSource(LR_pair,
               liga_data,
               liga_info,
               liga_control,
               liga_exper,
               liga_source_data,
               liga_source_info,
               liga_source_control,
               liga_source_exper,
               liga_source_expr_cutoff=0.25):
    """
    LigaSource程序的执行函数。

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

    The execution function of LigaSource.

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

    # 利用配体在输入的可能来源细胞中的表达情况推测LR关系对中配体的来源。
    # The source of the ligand in the LR relationship pair was inferred 
    # by the expression of the ligand in the input potential source cells.
    LR_liga_source =  inferLigaSource.inferLigaSource(LR_pair,
                                                      liga_data,
                                                      liga_info,
                                                      liga_control,
                                                      liga_exper,
                                                      liga_source_data,
                                                      liga_source_info,
                                                      liga_source_control,
                                                      liga_source_exper,
                                                      liga_source_expr_cutoff)

    
    return LR_liga_source


def BubblePlot(data,
               top_n=30,
               figsize=(3.5,7),
               dpi=300,
               save='./bubbleplot.pdf',
               show_plot=True
               ):
    """
    BubblePlot程序的执行函数。
    
    data: ExtraCellTalk分析得到的预测结果，
    top_n: 30,
        按LR_score绘制前多少个显著的配体-受体关系对,
        默认数值为30。
    figsize: (3.5,7),
        图片尺寸, 分别为长和宽, 
        默认数值为 3.5 * 7。
    dpi: 300, 
        图片分辨率,
        默认数值为300。
    save: './bubbleplot.pdf',
        图片保存位置,
        默认位置为'./bubbleplot.pdf'。
    show_plot: True,
        是否在Jupyter Notebook上展示图片,
        默认展示。

    The execution function of BubblePlot.

    data: the prediction results obtained by ExtraCellTalk,
    top_n: 30,
        Plot the top n significant ligand-receptor relationship pairs 
        based on LR_score,
        The default value is 30.
    figsize: (3.5,7),
        Picture size, respectively length and width,
        The default value is 3.5 * 7.
    dpi: 300,
        Picture resolution,
        The default value is 300.
    save: './bubbleplot.pdf',
        Picture save location,
        The default location is './bubbleplot.pdf'.
    show_plot: True,
        Whether to display pictures on Jupyter Notebook,
        The default is True.
    """

    # 利用ExtraCellTalk计算得到的结果绘制配体-受体关系对气泡图。
    # Bubble plot of ligand-receptor relationship pairs 
    # based on the results calculated by ExtraCellTalk.
    plotRes.BubblePlot(data,
                       top_n,
                       figsize,
                       dpi,
                       save,
                       show_plot)


def SankeyPlot(data,
               top_n=30,
               figsize=(6,7),
               save='./sankeyplot.pdf'
               ):
    """
    SankeyPlot程序的执行函数。

    data: LigaSource分析得到的预测结果，
    top_n: 30,
        按LR_score绘制前多少个显著的配体来源-配体-受体-受体靶标关系对,
        默认数值为30。
    figsize: (6,7),
        图片尺寸, 分别为长和宽, 
        默认数值为 6 * 7。
    save: './sankeyplot.pdf',
        图片保存位置,
        默认位置为'./sankeyplot.pdf'。

    The execution function of SankeyPlot.

    data: the prediction results obtained by LigaSource,
    top_n: 30,
        Plot the top n significant ligand source-ligand-receptor-receptor target 
        relationship pairs based on LR_score,
        The default value is 30.
    figsize: (6,7),
        Picture size, respectively length and width,
        The default value is 6 * 7.
    save: './sankeyplot.pdf',
        Picture save location,
        The default location is './sankeyplot.pdf'.
    """

    # 利用LigaSource计算得到的结果绘制配体来源-配体-受体-受体靶标关系对桑基图。
    # Sankey plot of ligand source-ligand-receptor-receptor target 
    # relationship pairs based on the results calculated by LigaSource.
    plotRes.SankeyPlot(data,
                       top_n,
                       figsize,
                       save)

