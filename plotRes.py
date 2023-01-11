import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

warnings.filterwarnings("ignore")
pandas2ri.activate()


## 配体-受体关系对气泡图
## Bubble plot of ligand-receptor relationship pairs
def BubblePlot(data,
               top_n,
               figsize,
               dpi,
               save,
               show_plot
               ):
    """
    根据预测结果绘制配体-受体关系对气泡图。
    
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

    Bubble plot of ligand-receptor relationship pairs based on predicted results.

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
    
    # 根据LR score，提取前n个配体-受体关系对进行绘图
    # According to the LR score, 
    # select top n ligand-receptor relationship pairs for plotting
    df = data
    df.sort_values(by='LR_score', inplace=True, ascending=False)
    df = df.head(n=top_n)
    
    # 受体排序
    # Sort by receptor
    df.sort_values(by='receptors', inplace=True)
    
    # 对FDR值标准化
    # Normalize the FDR value
    FDR_value = -np.log10(df['FDR'])
    FDR_normalized_value = [i/np.max(FDR_value)*100 for i in FDR_value.values]

    
    # 绘图尺寸
    # Plot size
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # 绘图
    # Plot
    scatter = ax.scatter(df['celltype'],
                         df['interaction_name'],
                         s=FDR_normalized_value,
                         c=np.log10(df['LR_score']))
    ax.tick_params(axis='x',which='major',rotation=45,pad=0.15)
    ax.tick_params(axis='y',which='major',pad=0.15)
    ax.set_xlabel('')
    ax.set_ylabel('')

    # 设置网格线透明度
    # Set grid lines with some transparency
    ax.grid(alpha=0.4)

    # 确保网格线位于其他对象之后
    # Make sure grid lines are behind other objects
    ax.set_axisbelow(True)

    # 默认垂直限制缩小 0.5
    # Default vertical limits are shrunken by 0.5
    y_shrunk = 0.5
    y_lower, y_upper = ax.get_ylim()
    ax.set_ylim(y_lower + y_shrunk, y_upper - y_shrunk)

    # 点大小图例
    # Legend for dot size
    dot_size = sorted(set([np.min(FDR_value),
              np.percentile(FDR_value, 50),
              np.max(FDR_value)]))
    for label in dot_size:
        ax.scatter([], [], 
                   s = label/np.max(FDR_value)*100,
                   facecolors='black',
                   label = round(label, 2))
    ax.legend(title = '-log10(FDR)', 
              loc='lower left', 
              bbox_to_anchor=(0.98,0),
              frameon = False)

    # 颜色图例
    # Legend for color
    cbar_pos = fig.add_axes([0.96, 0.3, 0.04, 0.5])
    cbar = plt.colorbar(scatter, cax=cbar_pos)
    cbar.set_label(label = 'log10(LR score)')
    
    # 储存图片
    # Save plot
    plt.savefig(save,bbox_inches='tight')
    
    # 判断是否在jupyter notebook中展示
    # Determine whether to show in jupyter notebook
    if show_plot:
        plt.show()


## 配体来源-配体-受体-受体靶标关系对桑基图，调用R语言绘制
## Sankey plot of the ligand source-ligand-receptor-receptor target 
## relationship pairs is drawn by calling R language
def SankeyPlot(data,
               top_n,
               figsize,
               save
               ):
    """
    根据推断的配体来源，绘制配体来源-配体-受体-受体靶标关系对桑基图。

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

    Sankey plot of ligand source-ligand-receptor-receptor target 
    relationship pairs based on predicted results.

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

    # R语言绘图程序
    # R drawing program
    robjects.r('''
    # 导入需要的R包
    # Import the required R packages
    library(ggplot2)
    library(ggalluvial)
    library(RColorBrewer)

    # 桑基图绘图函数
    # Sankey plot function
    sankey_plot <- function(sankey_data, 
                            width,
                            height,
                            save_path){
        plt <- ggplot(sankey_data, 
                      aes(axis1=liga_source, 
                          axis2=interaction_name.1, 
                          axis3=celltype, 
                          y=sankey_count))+
                geom_alluvium(aes(fill=interaction_name.1), 
                              width=1/8, 
                              knot.pos=0.1, 
                              reverse=F) +
                guides(fill=FALSE) +
                geom_stratum(alpha=0.8, 
                             width=1/8, 
                             reverse=F) +
                scale_fill_manual(values=brewer.pal(12,'Paired')) +
                geom_text(stat="stratum", 
                          aes(label=after_stat(stratum)),
                              reverse=FALSE, 
                              size=3) +
                scale_x_discrete(limits=c("Ligand source",
                                          "LR pair",
                                          "Receptor target"), 
                                 expand=c(0, 0))+
                xlab("") + 
                ylab("") +
                theme_bw() +
                theme(panel.grid=element_blank(),
                      panel.border=element_blank(),
                      axis.line=element_blank(),
                      axis.ticks=element_blank(),
                      axis.text.y=element_blank(),
                      axis.text.x=element_text(colour = 'black'))
        
        # 保存绘图结果
        # Save drawing results
        ggsave(filename=save_path, 
               width=width, 
               height=height, 
               plot=plt, 
               device="pdf")
    }
    ''')

    # 根据LR score，提取前n个配体来源-配体-受体-受体靶标关系对进行绘图
    # According to the LR score, 
    # select top n ligand source-ligand-receptor-receptor target 
    # relationship pairs for plotting
    df = data
    df.sort_values(by='LR_score', inplace=True, ascending=False)
    df = df.head(n=top_n)

    # 计算桑基图流的宽
    # Calculate the thickness of Sankey plot flow
    df['sankey_count'] = np.log10(df['liga_source_score'] * df['LR_score'] + 1)

    # 转换python数据格式为R数据格式
    # Convert python data format to R data format
    sankey_data = pandas2ri.py2rpy(df)
    width = robjects.FloatVector([figsize[0]])
    height = robjects.FloatVector([figsize[1]])
    save_path = robjects.StrVector([save])

    # 调用R语言绘图程序
    # Call the R drawing program
    robjects.r.sankey_plot(sankey_data,width,height,save_path)
