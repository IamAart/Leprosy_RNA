import pandas as pd

analysis = pd.read_csv('data/Predictions/analysis_rf_predictions.csv')
library_types = [["DESeq2", "EdgeR"], ["DESeq2", "LimmaVoom"], ["EdgeR", "LimmaVoom"]]
final_dataframe = pd.DataFrame(columns=["DGE Method", "DGE Method 2", "Average overlap"])
for libraries in library_types:
    subset = analysis[analysis["Library"].isin(libraries)].sort_values(by=["Number of genes", "Feature Selection Method"])
    average_overlap = 0
    for i in range(0, len(subset), 2):
        first = subset.iloc[i]["Ensemble Ids"].split(", ")
        second = subset.iloc[i+1]["Ensemble Ids"].split(", ")
        diff_list = list(set(first) - set(second))
        if diff_list == []:
            new_value = 100
        else:
            new_value = (len(first)-len(diff_list))/len(first)*100
        average_overlap += new_value

    final_dataframe = pd.concat([final_dataframe, pd.DataFrame([[libraries[0], libraries[1], average_overlap/(len(subset)/2), lowest_overlap]], columns=["DGE Method", "DGE Method 2", "Average overlap", "Lowest overlap"])], axis=0)

final_dataframe.to_excel("data/overlap_analysis_dge_method.xlsx", index=False)