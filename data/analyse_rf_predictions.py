import os

import dotenv
import pandas as pd
from dotenv import load_dotenv

load_dotenv()
LIBRARY_COMBINATIONS = eval(os.getenv("LIBRARY_COMBINATIONS").replace('"', ''))
MIN = int(os.getenv("MIN"))

def best_dge_genes(data, dge_path, sorted_by_p_value):
    """
    This function will return the best genes from the DGE analysis. With best genes, we mean the genes that are in the
    top 10% of the AUC scores. The genes will be sorted by p-adjust-value if sorted_by_p_value is True, otherwise they will be
    sorted by the occurrence of a gene in the top 10% of the AUC scores.
    :param data: the data from the random forest analysis
    :param dge_path: the path to the DGE analysis dataset
    :param sorted_by_p_value: a boolean to indicate if the genes should be sorted by p-adjust-value
    :return: a list of the best genes with their corresponding ensemble IDs, gene names and bio-types
    """
    data = data.sort_values("Average AUC score", ascending=False)
    data = data[data["Average AUC score"] > data["Average AUC score"].quantile(0.9)]
    ensembles_dict = dict()
    for _, row in data.iterrows():
        e_list = row["Ensemble Ids"].split(", ")
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
    results.to_excel("./Predictions_current/best_dge_genes.xlsx")
    return results["Row.names"].tolist(), results["gene_name"].tolist(), results["biotype"].tolist()


def combine_one_run(dictionary, path, analysis_type, library, feature_selection, min_genes):
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
        accuracy_score, auc_score = 0, 0
        for _, value in row["results"].items():
            accuracy_score += value["accuracy"]
            auc_score += value["roc_auc"]
        # store every data in a list
        dictionary["Analysis type"].append(analysis_type)
        dictionary["Gene names"].append(row["genes"])
        dictionary["Ensemble Ids"].append(row["ensembles"])
        dictionary["Average accuracy"].append(accuracy_score / len(row["results"].keys()))
        dictionary["Average AUC score"].append(auc_score / len(row["results"].keys()))
        dictionary["Library"].append(library)
        dictionary["Feature Selection Method"].append(feature_selection)
        dictionary["Number of genes"].append(gene_amount)

        gene_amount += 1
    return dictionary



def combine_rf_results(analysis_type_list: list, used_subset: bool):
    """
    This function will combine the results of the random forest analysis into a pandas DataFrame containing the
    run parameters, use features (gene_names, ensemble_ids) and the average accuracy and AUC score of the amount of
    LOOCV iterations.
    :return: a pandas DataFrame containing the results of the random forest analysis
    """
    # prepare for data gathering from jsons
    result_dict = {
        "Analysis type": [],
        "Library": [],
        "Number of genes": [],
        "Feature Selection Method": [],
        "Average accuracy": [],
        "Average AUC score": [],
        "Ensemble Ids": [],
        "Gene names": []

    }
    for analysis_type in analysis_type_list:
        if used_subset:
            for key in LIBRARY_COMBINATIONS.keys():
                for feature_type in ["RFE", "chi2"]:
                    # create a dataframe from the result json and loop over each score and ensembles /genes
                    json_path = f"./Predictions/{analysis_type}/{key}/{feature_type}/results.json"
                    result_dict = combine_one_run(result_dict, json_path, analysis_type, key, feature_type, MIN)
        else:
            json_path = f"./Predictions/{analysis_type}/results.json"
            result_dict = combine_one_run(result_dict, json_path, analysis_type, "SUBSET", os.getenv("FEATURE_SELECTION_METHOD"), MIN)

    return pd.DataFrame(result_dict)

# save data in csv and excel
result_df = combine_rf_results(["All_GENES", "CODING", "NON_CODING"], True)
result_df.to_csv("./Predictions/analysis_rf_predictions.csv")
result_df.to_excel("./Predictions/analysis_rf_predictions.xlsx")

# call the function to get the best genes
result_df = pd.read_csv("./Predictions/analysis_rf_predictions.csv")
ensembles, gene_names, bio_types = best_dge_genes(result_df, "./LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx", False)

# Set the ensembles and gene names in the .env file
dotenv_file = dotenv.find_dotenv(".env")
dotenv.set_key(dotenv_file, "SUBSET_ENSEMBLES", ",".join(ensembles))
dotenv.set_key(dotenv_file, "SUBSET_GENE_NAMES", ",".join(gene_names))

# IF RF analysis on the best genes is done
if os.getenv("SUBSET_MODEL") in ["RFE", "chi2"]:
    best_result_df = combine_rf_results(["SUBSET"], True)
    best_result_df.to_csv("./Predictions/analysis_rf_predictions_best_genes.csv")
    best_result_df.to_excel("./Predictions/analysis_rf_predictions_best_genes.xlsx")
