import pandas as pd
# TODO: PUSH THIS TO MAIN


def find_gene_name(ensemble, data):
    return data[data["Row.names"] == ensemble].iloc[0]["gene_name"]


def combine_two_algo_data(path_1, path_2, name_1, name_2):
    data_1 = pd.read_excel(path_1)
    print(data_1.columns)
    data_2 = pd.read_excel(path_2)
    # if not ({"gene_type", "feature_selection", "group_size"}.issubset(set(data_1.columns))):
    #     raise Exception("Dataset of Algorithm 1 doesn't have the required columns (gene_type, feature_selection, group_size)")
    # if not ({"gene_type", "feature_selection", "group_size"}.issubset(set(data_2.columns))):
    #     raise Exception("Dataset of Algorithm 2 doesn't have the required columns (gene_type, feature_selection, group_size)")

    data_1 = data_1[~data_1["Library"].isin(["All", "DESeq2", "LimmaVoom", "EdgeR"])]
    data_1["Selection"] = data_1["gene_type"] + "-" + data_1["feature_selection"] + "-" + data_1["group_size"].astype(str) #  + "-" + data_2["libraries"]
    data_2["Selection"] = data_2["gene_type"] + "-" + data_2["feature_selection"] + "-" + data_2["group_size"].astype(str) #  + "-" + data_2["libraries"]

    #  create the list that will gather the data for the excelsheet
    accuracy_list_1 = list()
    accuracy_list_2 = list()
    auc_list_1 = list()
    auc_list_2 = list()
    ensembles_1 = list()
    svm_ensembles = list()
    combined_accuracy_list = list()
    combined_auc_list = list()
    intersect_ensembles = list()
    only_ensembles_1 = list()
    only_ensembles_2 = list()
    intersect_gene_names = list()
    only_gene_name_1 = list()
    only_gene_name_2 = list()

    rna_data = pd.read_excel("./DESeq2/All_GENES_DESeq2_RNA_SEQ.xlsx")
    # collect the data from the rows and put them into the correct lists
    for selection in data_1["Selection"].tolist():
        row_1 = data_1[data_1["Selection"] == selection]
        row_2 = data_2[data_2["Selection"] == selection]
        accuracy_list_1.append(float(row_1["accuracy"].values[0]))
        accuracy_list_2.append(float(row_2["accuracy"].values[0]))
        auc_list_1.append(float(row_1["auc"].values[0]))
        auc_list_2.append(float(row_2["auc"].values[0]))
        ensembles_1.append(row_1["ensembles"].values[0])
        svm_ensembles.append(row_2["ensembles"].values[0])
        combined_accuracy_list.append((float(row_1["accuracy"].values[0])+float(row_2["accuracy"].values[0]))/2)
        combined_auc_list.append((float(row_1["auc"].values[0])+float(row_2["auc"].values[0]))/2)
        ensembles_intersect = set(row_1["ensembles"].values[0].split(", ")) & set(row_2["ensembles"].values[0].split(", "))
        ensembles_only_1 = set(row_1["ensembles"].values[0].split(", ")) - set(row_2["ensembles"].values[0].split(", "))
        ensembles_only_2 = set(row_2["ensembles"].values[0].split(", ")) - set(row_1["ensembles"].values[0].split(", "))
        intersect_ensembles.append(ensembles_intersect)
        only_ensembles_1.append(ensembles_only_1)
        only_ensembles_2.append(ensembles_only_2)

        intersect_genes = list()
        for ensemble in list(ensembles_intersect):
            intersect_genes.append(find_gene_name(ensemble, rna_data))

        only_1_genes = list()
        for ensemble_1 in list(ensembles_only_1):
            only_1_genes.append(find_gene_name(ensemble_1, rna_data))

        only_2_genes = list()
        for ensemble_2 in list(ensembles_only_2):
            only_2_genes.append(find_gene_name(ensemble_2, rna_data))

        intersect_gene_names.append(intersect_genes)
        only_gene_name_1.append(only_1_genes)
        only_gene_name_2.append(only_2_genes)

    try:
        # create a dataset for the combined dataset, this contains all data and no calculation for averages
        dataset = pd.DataFrame(
            zip(data_1["Selection"].tolist(), accuracy_list_1, auc_list_1, ensembles_1, accuracy_list_2, auc_list_2,
                svm_ensembles),
            columns=["selection", f"{name_1}_accuracy", f"{name_1}_auc", f"{name_1}_ensembles", f"{name_2}_accuracy",
                     f"{name_2}_auc", f"{name_2}_ensembles"]
        )
        dataset.to_excel("combined_data.xlsx")
    except Exception as e:
        raise Exception(e)

    try:
        # create a dataset to combine the data of the two algorithms to get intersect, and algo 1 and2 only ensembles
        # and their combined accuracy and auc
        combined_data = pd.DataFrame(
            zip(data_1["Selection"].tolist(), combined_accuracy_list, combined_auc_list, intersect_ensembles,
                only_ensembles_1, only_ensembles_2, intersect_gene_names, only_gene_name_1, only_gene_name_2),
            columns=["selection", "accuracy", "auc", "intersect ensembles", f"{name_1} ensembles",
                     f"{name_2} ensembles", "intersect gene names", f"{name_1} gene names", f"{name_2} gene names"
            ]
        )
        combined_data.to_excel("calculated_combined_data.xlsx")
    except Exception as e:
        raise Exception(e)


def create_excelsheet(input_path, output_path, genes):
    data = pd.read_csv(input_path)
    data[data["Row.names"].isin(genes)].sort_values("P.Adjust").to_excel(output_path)


def get_best_combined_intersects(combined_data, extra_genes=None):
    # prepare for checking each selection for best auc score
    list_selection_no_digits = list()
    digits = list()
    for item in combined_data["selection"].tolist():
        list_selection_no_digits.append(''.join(filter(lambda x: not x.isdigit(), item)))
        digits.append(item.split("-")[2])
    combined_data["selection"] = list_selection_no_digits
    combined_data["amount_of_genes"] = digits

    # check each selection and get the intersect genes of best auc score
    nc_genes = list()
    all_genes = list()
    for selection_option in set(list_selection_no_digits):
        selection_data = combined_data[combined_data["selection"] == selection_option]

        if selection_option.startswith("nc-"):
            print(selection_option)
            print(selection_data.sort_values(by="auc", ascending=False).iloc[0]["amount_of_genes"])
            best_intersect = selection_data.sort_values(by="auc", ascending=False).iloc[0]["intersect ensembles"]
            best_intersect = best_intersect.replace("{", "").replace("}", "").replace("'", "").split(", ")
            print(best_intersect)
            nc_genes.extend(best_intersect)
        elif selection_option.startswith("all-"):
            print(selection_option)
            print(selection_data.sort_values(by="auc", ascending=False).iloc[0]["amount_of_genes"])
            best_intersect = selection_data.sort_values(by="auc", ascending=False).iloc[0]["intersect ensembles"]
            best_intersect = best_intersect.replace("{", "").replace("}", "").replace("'", "").split(", ")
            print(best_intersect)
            all_genes.extend(best_intersect)
        elif selection_option.startswith("coding-"):
            print(selection_option)
            print(selection_data.sort_values(by="auc", ascending=False).iloc[0]["amount_of_genes"])
            best_intersect = selection_data.sort_values(by="auc", ascending=False).iloc[0]["intersect ensembles"]
            best_intersect = best_intersect.replace("{", "").replace("}", "").replace("'", "").split(", ")
            print(best_intersect)
            all_genes.extend(best_intersect)
        else:
            raise Exception("Selection doesn't start with all or nc")

    # if requested add extra genes
    if extra_genes:
        all_genes.extend(extra_genes)

    all_genes = list(set(all_genes))
    all_genes = list(set(all_genes) - set(nc_genes))
    nc_genes = list(set(nc_genes) - set(all_genes))
    print(len(all_genes))
    print(len(nc_genes))
    print(list(set(nc_genes+all_genes)))

    # create excelsheet for each library
    create_excelsheet("./DESeq2/All_GENES_DESeq2_RNA_SEQ.csv", "./deseq2_all_results.xlsx", list(set(nc_genes+all_genes)))
    create_excelsheet("./LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.csv", "./limma_all_results.xlsx", list(set(nc_genes+all_genes)))
    create_excelsheet("./EdgeR/All_GENES_EdgeR_RNA_SEQ.csv", "./edger_all_results.xlsx", list(set(nc_genes+all_genes)))


# combine_two_algo_data(
#     "./Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx",
#     "./Current_Comparison_SVM_RF/svm_analysis_9_22.xlsx",
#     "Random Forest",
#     "Support Vector Machine"
# )

data = pd.read_excel("./Current_Comparison_SVM_RF/calculated_combined_data.xlsx")

#  UBC, TPGS1, C19orf60, MT-ND2
genes_added = ["ENSG00000150991", "ENSG00000141933", "ENSG00000006015", "ENSG00000198763"]
get_best_combined_intersects(data, genes_added)



