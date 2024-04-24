import pandas as pd

COLUMN_NAME = "gene_name"  # "Row.names"
GENE_TYPE = "All_GENES"  # "NON_CODING

deseq = pd.read_csv(f"../DESeq2/{GENE_TYPE}_DESeq2_RNA_SEQ.csv")
limma = pd.read_csv(f"../LimmaVoom/{GENE_TYPE}_LimmaVoom_RNA_SEQ.csv")
edger = pd.read_csv(f"../EdgeR/{GENE_TYPE}_EdgeR_RNA_SEQ.csv")


def sort_on_p_adj_val(genes, data, name):
    value_list = list()
    for gene in genes:
        temp_list = list()
        for i in range(len(data)):
            info = data[i][data[i].gene_name == gene]["P.Adjust"].values[0]
            temp_list.append(info)
        value_list.append(sum(temp_list)/len(temp_list))
    print(value_list, genes)
    sorted_value_list = [gene for value, gene in sorted(zip(value_list, genes))]
    print(sorted_value_list)
    return pd.Series(sorted_value_list, name=name)


def merge_data(data_deseq, data_limma, data_edger):
    set_deseq = data_deseq[data_deseq["col"] != "not"][COLUMN_NAME].to_list()
    set_limma = data_limma[data_limma["col"] != "not"][COLUMN_NAME].to_list()
    set_edger = data_edger[data_edger["col"] != "not"][COLUMN_NAME].to_list()

    set_all_3 = sort_on_p_adj_val(
        list(set(set_deseq) & set(set_limma) & set(set_edger)),
        [data_deseq, data_limma, data_edger],
        "DESeq2/LimmaVoom/EdgeR"
    )
    set_deseq_edger = sort_on_p_adj_val(
        list(set(set_deseq) & set(set_edger)), [data_deseq, data_edger], "DESeq2/EdgeR"
    )
    set_deseq_limma = sort_on_p_adj_val(
        list(set(set_deseq) & set(set_limma)), [data_deseq, data_limma], "DESeq2/LimmaVoom"
    )
    set_limma_edger = sort_on_p_adj_val(
        list(set(set_limma) & set(set_edger)), [data_edger, data_limma], "EdgeR/LimmaVoom"
    )
    set_deseq = sort_on_p_adj_val(set_deseq, [data_deseq], "DESeq2")
    set_limma = sort_on_p_adj_val(set_limma, [data_limma], "LimmaVoom")
    set_edger = sort_on_p_adj_val(set_edger, [data_edger], "EdgeR")

    merged_dataset = pd.concat(
        [set_all_3, set_deseq_edger, set_deseq_limma, set_limma_edger, set_deseq, set_limma, set_edger],
        axis=1
    )

    merged_dataset.to_csv(f"{GENE_TYPE}_VENN_DATA.csv", index=False)
    with pd.ExcelWriter(f"{GENE_TYPE}_VENN_DATA.xlsx") as writer:
        merged_dataset.to_excel(writer, index=False)

merge_data(deseq, limma, edger)