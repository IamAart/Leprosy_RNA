import os
import pandas as pd
from dotenv import load_dotenv

def dataset_with_selection(data, selection=None):
    """
    Filter the dataset based on the selection of biotypes
    :param data: the dataset with the differentially expressed genes from a specific method
    :param selection: the biotypes to filter on, specifically Non-coding, Coding or All genes
    :return: the same dataset filtered on the biotype, thus only specific genes are left
    """
    if selection is None:
        return data
    return data.loc[data["biotype"].isin(selection)]


def sort_on_p_adj_val(ensembles, data, name):
    """
    Sort the ensembl IDs based on the average p.adjust value between the three methods
    :param ensembles: a list of emsembl IDs
    :param data: a dataset of the differentially expressed genes from the three methods
    :param name: a string which is used as name for the column in the final dataframe, of one of the three methods
    :return: a pandas Series object with the ensemble IDs and gene names sorted on the average p.adjust value
    """
    p_adj_list = list()
    genes_list = list()

    for ensemble in ensembles:
        p_adj_value = 0
        added_gene_flag = False
        for i in range(len(data)):
            if ensemble in data[i].index:
                p_adj_value += float(data[i].loc[ensemble]["P.Adjust"])
                if not added_gene_flag:
                    genes_list.append(data[i].loc[ensemble]["gene_name"])
                    added_gene_flag = True
        p_adj_list.append(p_adj_value/len(data))
    ensemble_series = pd.Series([ensemble for _, ensemble in sorted(zip(p_adj_list, ensembles))],name=name).reset_index(drop=True)
    gene_name_list = pd.Series([gene for _, gene in sorted(zip(p_adj_list, genes_list))],name=name).reset_index(drop=True)

    return ensemble_series, gene_name_list


def merge_data(data_deseq, data_limma, data_edger, biotype):
    """
    Obtain ensembl IDs and genes of the intersection of the three methods, and the intersection of each combination of
    two methods and each method alone
    :param data_deseq: differentially expressed genes data from DESeq2
    :param data_limma: differentially expressed genes data from LimmaVoom
    :param data_edger: differentially expressed genes data from EdgeR
    :param biotype: the biotype of the genes, either Non-coding, Coding or All genes
    :return: None
    """
    set_deseq = data_deseq[data_deseq["col"] != "not"].index.tolist()
    set_limma = data_limma[data_limma["col"] != "not"].index.tolist()
    set_edger = data_edger[data_edger["col"] != "not"].index.tolist()

    # get intersects of the 3 libraries
    set_all_3_ens, set_all_3_gene = sort_on_p_adj_val(
        list(set(set_deseq) & set(set_limma) & set(set_edger)),
        [data_deseq, data_limma, data_edger],
        "DESeq2/LimmaVoom/EdgeR"
    )
    set_deseq_edger_ens, set_deseq_edger_gene = sort_on_p_adj_val(
        list(set(set_deseq) & set(set_edger) - set(list(set(set_deseq) & set(set_limma) & set(set_edger)))),
        [data_deseq, data_edger],"DESeq2/EdgeR"
    )
    set_deseq_limma_ens, set_deseq_limma_gene = sort_on_p_adj_val(
        list(set(set_deseq) & set(set_limma) - set(list(set(set_deseq) & set(set_limma) & set(set_edger)))),
        [data_deseq, data_limma], "DESeq2/LimmaVoom"
    )
    set_limma_edger_ens, set_limma_edger_gene = sort_on_p_adj_val(
        list(set(set_limma) & set(set_edger) - set(list(set(set_deseq) & set(set_limma) & set(set_edger)))),
        [data_edger, data_limma], "EdgeR/LimmaVoom"
    )
    set_deseq2_ens, set_deseq2_gene = sort_on_p_adj_val(set(set_deseq) - set(set_limma) - set(set_edger),
                                                        [data_deseq], "DESeq2"
    )
    set_limma2_ens, set_limma2_gene = sort_on_p_adj_val(set(set_limma) - set(set_edger) - set(set_deseq),
                                                        [data_limma], "LimmaVoom"
    )
    set_edger2_ens, set_edger2_gene = sort_on_p_adj_val(set(set_edger) - set(set_limma) - set(set_deseq),
                                                        [data_edger], "EdgeR"
    )

    ensemble_lists = [set_all_3_ens, set_deseq_edger_ens, set_deseq_limma_ens, set_limma_edger_ens, set_deseq2_ens, set_limma2_ens, set_edger2_ens]
    genes_lists = [set_all_3_gene, set_deseq_edger_gene, set_deseq_limma_gene, set_limma_edger_gene, set_deseq2_gene, set_limma2_gene, set_edger2_gene]

    # create combined dataframe
    merged_dataset_ens = pd.concat(ensemble_lists,axis=1)
    merged_dataset_genes = pd.concat(genes_lists,axis=1)

    merged_dataset_ens.to_csv(f"./Venn_Diagrams/Ensemble_data/ENSEMBLE_{biotype}_VENN_DATA.csv", index=False)
    with pd.ExcelWriter(f"./Venn_Diagrams/Ensemble_data/ENSEMBLE_{biotype}_VENN_DATA.xlsx") as writer:
        merged_dataset_ens.to_excel(writer, index=False)

    merged_dataset_genes.to_csv(f"./Venn_Diagrams/Gene_name_data/GENES_{biotype}_VENN_DATA.csv", index=False)
    with pd.ExcelWriter(f"./Venn_Diagrams/Gene_name_data/GENES_{biotype}_VENN_DATA.xlsx") as writer:
        merged_dataset_genes.to_excel(writer, index=False)

load_dotenv()

COLUMN_NAME = "Row.names"
NON_CODING_BIOTYPES = os.getenv("NON_CODING_BIOTYPES").split(",")
CODING_BIOTYPES = os.getenv("CODING_BIOTYPES").split(",")

deseq = pd.read_csv(f"./DESeq2/All_GENES_DESeq2_RNA_SEQ.csv", index_col=1)
limma = pd.read_csv(f"./LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.csv", index_col=1)
edger = pd.read_csv(f"./EdgeR/All_GENES_EdgeR_RNA_SEQ.csv", index_col=1)

for biotype, selection in {"NON_CODING": NON_CODING_BIOTYPES, "CODING": CODING_BIOTYPES, "All_GENES": None}.items():
    deseq_selection = dataset_with_selection(deseq, selection)
    edger_selection = dataset_with_selection(edger, selection)
    limma_selection = dataset_with_selection(limma, selection)
    merge_data(deseq_selection, limma_selection, edger_selection, biotype)
