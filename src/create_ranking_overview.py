import pandas as pd
import os

def create_ranking():
    for subdir, dirs, files in os.walk("../data/Predictions/"):
        output_dict_ranking = {"ensembl": [], "gene_name": []}
        output_dict_weights = {"ensembl": [], "gene_name": []}
        if "kNN" in subdir: continue
        if len(subdir.split("Predictions")[1].split("\\")) == 4:   
            for i in range(3, 21):
                file_name = str(i) + "_weights.csv"
                print(subdir)
                weights_dataframe = pd.read_csv(subdir + "/" + file_name)

                weights_dataframe["ranking"] = range(1, len(weights_dataframe) + 1)
                for _, row in weights_dataframe.iterrows():
                    if row["selector"] not in output_dict_weights["ensembl"] and row["selector"] not in output_dict_ranking["ensembl"]:
                        output_dict_ranking["ensembl"].append(row["selector"])
                        output_dict_ranking["gene_name"].append(row["gene_names"])
                        output_dict_weights["ensembl"].append(row["selector"])
                        output_dict_weights["gene_name"].append(row["gene_names"])

                group_size = i
                output_dict_ranking[group_size] = list()
                output_dict_weights[group_size] = list()
                for ensemble in output_dict_ranking["ensembl"]:
                    if ensemble in weights_dataframe["selector"].values:
                        output_dict_ranking[group_size].append(int(weights_dataframe[weights_dataframe["selector"] == ensemble]["ranking"].values[0]))
                        output_dict_weights[group_size].append(float(weights_dataframe[weights_dataframe["selector"] == ensemble]["weights"].values[0]))
                    else:
                        output_dict_ranking[group_size].append(None)
                        output_dict_weights[group_size].append(None)

                total_length = len(output_dict_ranking["ensembl"])
                for key in output_dict_ranking.keys():
                    if isinstance(key, int):
                        if len(output_dict_ranking[key]) != total_length:
                            output_dict_ranking[key].extend([None] * (total_length - len(output_dict_ranking[key])))
                            output_dict_weights[key].extend([None] * (total_length - len(output_dict_weights[key])))
        
                if i == 20:  
                    pd.DataFrame.from_dict(output_dict_weights).to_csv(subdir + "/weights_ranking.csv", index=False)
                    pd.DataFrame.from_dict(output_dict_ranking).to_csv(subdir + "/ranking.csv", index=False)

create_ranking()