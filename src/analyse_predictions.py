import os
import pandas as pd
from dotenv import load_dotenv
import dotenv

load_dotenv()
LIBRARY_COMBINATIONS = eval(os.getenv("LIBRARY_COMBINATIONS").replace('"', ''))
MIN = int(os.getenv("MIN"))

def best_dge_genes(data, dge_path, sorted_by_p_value, model, dge_method):
    """
    This function will return the best genes from the DGE analysis. With best genes, we mean the genes that are in the
    top 10% of the AUC scores. The genes will be sorted by p-adjust-value if sorted_by_p_value is True, otherwise they will be
    sorted by the occurrence of a gene in the top 10% of the AUC scores.
    :param data: the data from the random forest analysis
    :param dge_path: the path to the DGE analysis dataset
    :param sorted_by_p_value: a boolean to indicate if the genes should be sorted by p-adjust-value
    :return: a list of the best genes with their corresponding ensemble IDs, gene names and bio-types
    """
    data = data.sort_values("auc", ascending=False)
    data = data[data["auc"] > 0.96] 
    data = data[data["model"].isin(model)]
    if dge_method: 
        data = data[data["dge_method"] == "Intersection"]
    ensembles_dict = dict()
    for _, row in data.iterrows():
        e_list = row["ensembl"].split(", ")
        for i in range(len(e_list)):
            if e_list[i] not in ensembles_dict.keys(): ensembles_dict[e_list[i]] = 1
            else: ensembles_dict[e_list[i]] += 1

    ensemble_list = sorted(ensembles_dict, key=lambda k: ensembles_dict[k], reverse=True)[0:20]
    dge_df = pd.read_excel(dge_path)
    results = dge_df[dge_df["Row.names"].isin(ensemble_list)]
    if sorted_by_p_value:
        results = results.sort_values("P.Adjust", ascending=True)
    else:
        results = results.set_index("Row.names").reindex(ensemble_list).reset_index()
    results = results[["gene_name", "Row.names", "biotype"]]
    results.rename(columns={"Row.names": "ensemble"}, inplace=True)
    results.to_excel(f"../data/Predictions/best_dge_genes_{model}.xlsx")
    return results["ensemble"].tolist(), results["gene_name"].tolist(), results["biotype"].tolist(), results


def combine_one_run(dictionary, path, min_genes):
    """
    This function will combine the results of one run of the random forest analysis into a pandas DataFrame containing the
    run parameters, use features (gene_names, ensemble_ids) and the average accuracy and AUC score of the amount of
    LOOCV iterations.
    :param dictionary: dictionary to store the data
    :param path: path to the json file containing the results of the random forest analysis
    :param analysis_type: the type of rna that was used for analysis
    :param library: the library/DGE method that was used for analysis
    :param feature_selection: the feature selection method that was used for analysis (RFE, chi2)
    :param min_genes: the minimum amount of genes that were used for analysis
    :return: a dictionary containing the results of the random forest analysis
    """
    df = pd.read_json(path).transpose()
    gene_amount = min_genes
    for _, row in df.iterrows():
        accuracy_score, auc_score, precision_score, recall_score, f1_score, specificity_score = 0, 0, 0, 0, 0, 0
        for _, value in row["results"].items():
            accuracy_score += value["accuracy"]
            auc_score += value["roc_auc"]
            precision_score += value["precision"]
            recall_score += value["recall"]
            f1_score += value["f1"]
            specificity_score += value["specificity"]

        # store every data in a list
        dictionary["analysis_id"].append(row["analysis_id"])
        dictionary["timestamp"].append(row["timestamp"])
        dictionary["model"].append(row["model"])
        dictionary["gene_type"].append(row["gene_type"])
        dictionary["dge_method"].append(row["dge_method"])
        dictionary["feature_selection"].append(row["feature_selection"])
        dictionary["group_size"].append(gene_amount)
        dictionary["accuracy"].append(accuracy_score / len(row["results"].keys()))
        dictionary["precision"].append(precision_score / len(row["results"].keys()))
        dictionary["recall"].append(recall_score / len(row["results"].keys()))
        dictionary["specificity"].append(specificity_score / len(row["results"].keys()))
        dictionary["f1"].append(f1_score / len(row["results"].keys()))
        dictionary["auc"].append(auc_score / len(row["results"].keys()))
        dictionary["ensembl"].append(row["ensembl"])
        dictionary["gene_name"].append(row["gene_names"])
        dictionary["hyperparameters"].append(row["hyperparameters"])

        gene_amount += 1
    return dictionary



def combine_rf_results(subgroup=None):
    """
    This function will combine the results of the random forest analysis into a pandas DataFrame containing the
    run parameters, use features (gene_names, ensemble_ids) and the average accuracy and AUC score of the amount of
    LOOCV iterations.
    :return: a pandas DataFrame containing the results of the random forest analysis
    """
    # prepare for data gathering from jsons
    result_dict = {
        "analysis_id": [],
        "timestamp": [],
        "model": [],
        "gene_type": [],
        "dge_method": [],
        "feature_selection": [],
        "group_size": [],
        "accuracy": [],
        "precision": [],
        "recall": [],
        "specificity": [],
        "f1": [],
        "auc": [],
        "ensembl": [],
        "gene_name": [],
        "hyperparameters": []
    }
    for subdir, dirs, files in os.walk("../data/Predictions"):
        for fname in files:
            if fname == "results.json" and dirs == ["Pictures"]:
                if not subgroup:
                    print(subdir, dirs, fname)  
                    result_dict = combine_one_run(result_dict, f"{subdir}/results.json", MIN)
                else:
                    if subgroup in subdir:
                        print(subdir, dirs, fname)
                        result_dict = combine_one_run(result_dict, f"{subdir}/results.json", MIN)


    return pd.DataFrame(result_dict)

# # Uncomment the following lines to analyse the predictions of the 4 different machine learning methods and save the results in analysis_predictions.xlsx
# result_df = combine_rf_results()
# analysis_ids = pd.read_excel("../data/analysis_ids_machine_learning.xlsx")["analysis_id"].tolist()
# result_df["analysis_id"] = result_df["analysis_id"].astype("category")
# result_df["analysis_id"] = result_df["analysis_id"].cat.set_categories(analysis_ids)
# result_df = result_df.sort_values(["analysis_id", "group_size"])
# result_df.to_csv("../data/Predictions/analysis_predictions.csv")
# result_df.to_excel("../data/Predictions/analysis_predictions.xlsx")

# # Uncomment the following lines to analyse analysis_predictions.xlsx to find the best genes and the intersection between the two paths (See Flowchart)
# result_df = pd.read_csv("../data/Predictions/analysis_predictions.csv")
# ensembles_rf, gene_names_rf, bio_types_rf, results_rf = best_dge_genes(result_df, "../data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx", False, ["RF"], False)
# ensembles_svm, gene_names_svm, bio_types_svm, results_svm = best_dge_genes(result_df, "../data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx", False, ["SVM", "RF"], True)
# ensembles_all = list(set(ensembles_svm) | set(ensembles_rf))
# ensembles_intersection = list(set(ensembles_rf) & set(ensembles_svm))
# for key, dataset in {'union': ensembles_all, 'intersection': ensembles_intersection}.items():
#     final_lists = list()
#     for idx, ensemble in enumerate(dataset):
#         if ensemble in ensembles_rf:
#             final_lists.append(results_rf[results_rf["ensemble"] == ensemble].values[0])
#         else:
#             final_lists.append(results_svm[results_svm["ensemble"] == ensemble].values[0])

#     final_df = pd.DataFrame(final_lists, columns=["gene_name", "ensemble", "biotype"])
#     final_df.to_excel(f"../data/Predictions/best_dge_genes_{key}.xlsx")

# dotenv_file = dotenv.find_dotenv(".env")
# dotenv.set_key(dotenv_file, "SUBSET_ENSEMBLES", ",".join(final_df["ensemble"].tolist()))
# dotenv.set_key(dotenv_file, "SUBSET_GENE_NAMES", ",".join(final_df["gene_name"].tolist()))
# dotenv.set_key(dotenv_file, "SUBSET_MIN", str(3))
# dotenv.set_key(dotenv_file, "SUBSET_MAX", str(len(final_df["ensemble"].tolist())))

# # Uncomment the following lines to analyse the predictions of the subset analysis and save the results in subset_analysis_predictions.xlsx
# result_df = combine_rf_results("Subset")
# result_df = result_df.sort_values(["group_size"])
# result_df.to_csv("../data/Predictions/subset_analysis_predictions.csv")
# result_df.to_excel("../data/Predictions/subset_analysis_predictions.xlsx")