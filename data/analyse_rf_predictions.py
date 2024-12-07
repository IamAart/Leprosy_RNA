import os
import pandas as pd
from dotenv import load_dotenv

load_dotenv()
LIBRARY_COMBINATIONS = eval(os.getenv("LIBRARY_COMBINATIONS").replace('"', ''))
MIN = os.getenv("MIN")

# prepare for data gathering from jsons
features = pd.read_csv('./sasc326_features.csv')
analysis_type_list, amount_list, gene_list, ensemble_list, accuracy_average_list, auc_average_list, library_list, feature_type_list = list(), list(), list(), list(), list(), list(), list(), list()
number_of_genes = list()
# loop through each directory to find json
for analysis_type in ["All_GENES", "CODING","NON_CODING"]:
    for key in LIBRARY_COMBINATIONS.keys():
        for feature_type in ["RFE", "chi2"]:
            # create a dataframe from the result json and loop over each score and ensembles /genes
            df = pd.read_json(f"./Predictions_current/{analysis_type}/{key}/{feature_type}/results.json").transpose()
            gene_amount = MIN
            for row_number, row in df.iterrows():
                accuracy_score, auc_score = 0, 0
                # calculate average accuracy and auc score
                for _, value in row["results"].items():
                    accuracy_score += value["accuracy"]
                    auc_score += value["roc_auc"]
                # store every data in a list
                analysis_type_list.append(analysis_type)
                gene_list.append(row["genes"])
                ensemble_list.append(row["ensembles"])
                accuracy_average_list.append(accuracy_score/len(row["results"].keys()))
                auc_average_list.append(auc_score/len(row["results"].keys()))
                library_list.append(key)
                feature_type_list.append(feature_type)

                number_of_genes.append(gene_amount)
                # add one more to the gene count, which resets on new combination
                gene_amount += 1



# create one final dataframe from all data
result_df = pd.DataFrame({
    "Analysis type": analysis_type_list,
    "Library": library_list,
    "Number of genes": number_of_genes,
    "Feature Selection Method": feature_type_list,
    "Average accuracy": accuracy_average_list,
    "Average AUC score": auc_average_list,
    "Ensemble Ids": ensemble_list,
    "Gene names": gene_list
})

# save data in csv and excel
result_df.to_csv("./Predictions_current/analysis_rf_predictions.csv")
result_df.to_excel("./Predictions_current/analysis_rf_predictions.xlsx")
