import warnings
import pandas as pd

warnings.filterwarnings("ignore")


def loadProtAnno(org):
    """
    加载Uniprot蛋白注释文件，即该蛋白是胞膜蛋白还是胞质蛋白。

    org: ['Human],
        研究的物种类型, 
        可选物种类型 Human, Mouse, Pig, 默认为 Human。

    Load the Uniprot protein annotation files to determine 
    whether the protein is a membrane protein or a cytoplasmic protein.

    org: ['Human],
        The type of species studied,
        Optional species type: Human, Mouse, Pig.
        Default is Human.
    """

    if org == 'Mouse':
        # 细胞膜蛋白 [KW-1003]
        # Cell Membrane Proteins [KW-1003]
        KW1003 = pd.read_excel("./data/database/uniprot-[KW-1003]-Mus+musculus.xlsx")
        KW1003.loc[KW1003['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
            KW1003.loc[KW1003['Gene names  (primary )'].isna(),:].apply(lambda row: \
            row['Entry name'].strip('_PIG'), axis=1)

        # 外周蛋白 [SL-9903]
        # Peripheral Proteins [SL-9903]
        SL9903 = pd.read_excel("./data/database/uniprot-[SL-9903]-Mus+musculus.xlsx")
        SL9903.loc[SL9903['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
            SL9903.loc[SL9903['Gene names  (primary )'].isna(),:].apply(lambda row: \
            row['Entry name'].strip('_PIG'), axis=1)

    elif org == 'Pig':
        # 细胞膜蛋白 [KW-1003]
        # Cell Membrane Proteins [KW-1003]
        KW1003 = pd.read_excel("./data/database/uniprot-[KW-1003]-Sus+scrofa.xlsx")
        KW1003.loc[KW1003['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
            KW1003.loc[KW1003['Gene names  (primary )'].isna(),:].apply(lambda row: \
            row['Entry name'].strip('_PIG'), axis=1)

        # 外周蛋白 [SL-9903]
        # Peripheral Proteins [SL-9903]
        SL9903 = pd.read_excel("./data/database/uniprot-[SL-9903]-Sus+scrofa.xlsx")
        SL9903.loc[SL9903['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
            SL9903.loc[SL9903['Gene names  (primary )'].isna(),:].apply(lambda row: \
            row['Entry name'].strip('_PIG'), axis=1)

    else:
        # 细胞膜蛋白 [KW-1003]
        # Cell Membrane Proteins [KW-1003]
        KW1003 = pd.read_excel("./data/database/uniprot-[KW-1003]-Homo+sapiens.xlsx")
        KW1003.loc[KW1003['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
            KW1003.loc[KW1003['Gene names  (primary )'].isna(),:].apply(lambda row: \
            row['Entry name'].strip('_PIG'), axis=1)

        # 外周蛋白 [SL-9903]
        # Peripheral Proteins [SL-9903]
        SL9903 = pd.read_excel("./data/database/uniprot-[SL-9903]-Homo+sapiens.xlsx")
        SL9903.loc[SL9903['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
            SL9903.loc[SL9903['Gene names  (primary )'].isna(),:].apply(lambda row: \
            row['Entry name'].strip('_PIG'), axis=1)

    # 从细胞膜蛋白中剔除外周蛋白，得到跨膜蛋白
    # Transmembrane proteins are obtained after including peripheral proteins 
    # are deleted from cell membrane proteins.
    transmembrane_prot = KW1003[~KW1003.Entry.isin(SL9903.Entry)]


    if org == 'Mouse':
        # 细胞质蛋白 [KW-0963]
        # Cytoplasm Proteins [KW-0963]
        KW0963 = pd.read_excel("./data/database/uniprot-[KW-0963]-Mus+musculus.xlsx")
        KW0963.loc[KW0963['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
        KW0963.loc[KW0963['Gene names  (primary )'].isna(),:].apply(lambda row: \
        row['Entry name'].strip('_PIG'), axis=1)

    elif org == 'Pig':
        # 细胞质蛋白 [KW-0963]
        # Cytoplasm Proteins [KW-0963]
        KW0963 = pd.read_excel("./data/database/uniprot-[KW-0963]-Sus+scrofa.xlsx")
        KW0963.loc[KW0963['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
        KW0963.loc[KW0963['Gene names  (primary )'].isna(),:].apply(lambda row: \
        row['Entry name'].strip('_PIG'), axis=1)

    else:
        # 细胞质蛋白 [KW-0963]
        # Cytoplasm Proteins [KW-0963]
        KW0963 = pd.read_excel("./data/database/uniprot-[KW-0963]-Homo+sapiens.xlsx")
        KW0963.loc[KW0963['Gene names  (primary )'].isna(),'Gene names  (primary )'] = \
        KW0963.loc[KW0963['Gene names  (primary )'].isna(),:].apply(lambda row: \
        row['Entry name'].strip('_PIG'), axis=1)

    # 细胞质蛋白
    # Cytoplasm Proteins
    cytoplasma_prot = KW0963


    return transmembrane_prot, cytoplasma_prot


def _unPackList(x):
    """
    去除列表中的列表。

    Remove the list from list.
    """

    temp_list = []
    for i in range(len(x)):
        temp_list = temp_list + x[i]
    

    return temp_list


def loadLRInfo():
    """
    加载自定义的配体-受体关系对，分为Single和Complex结构。

    Load custom ligand-receptor relationship pairs, 
    including the Single and Complex structures.
    """

    # 加载自定义的配体-受体关系对
    # Load custom ligand-receptor relationship pairs
    LR_sing = pd.read_excel('./data/database/LR_sing.xlsx')
    LR_comp = pd.read_excel('./data/database/LR_comp.xlsx')

    
    # 获得所有的配体和受体
    # Get all ligands and receptors from 
    # loaded ligand-receptor relationship pairs.
    LR_liga =  set([*list(LR_sing['ligands']), \
        *_unPackList(list(LR_comp['ligands'].str.split('_')))])
    LR_rece =  set([*list(LR_sing['receptors']), \
        *_unPackList(list(LR_comp['receptors'].str.split('_')))])


    return LR_sing, LR_comp, LR_liga, LR_rece


if __name__ == '__main__':
    loadProtAnno()
    loadLRInfo()
