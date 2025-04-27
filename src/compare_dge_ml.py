import pandas as pd

analysis = pd.read_csv('../data/Predictions/analysis_predictions.csv')
library_types = [["DESeq2", "EdgeR"], ["DESeq2", "LimmaVoom"], ["EdgeR", "LimmaVoom"]]
final_dataframe = pd.DataFrame(columns=["DGE Method", "DGE Method 2", "Average overlap"])
for libraries in library_types:
    subset = analysis[analysis["dge_method"].isin(libraries)].sort_values(by=["group_size", "feature_selection"])
    average_overlap = 0
    for i in range(0, len(subset), 2):
        first = subset.iloc[i]["ensembl"].split(", ")
        second = subset.iloc[i+1]["ensembl"].split(", ")
        diff_list = list(set(first) - set(second))
        if not diff_list:
            new_value = 100
        else:
            new_value = (len(first)-len(diff_list))/len(first)*100
        average_overlap += new_value

    final_dataframe = pd.concat([final_dataframe, pd.DataFrame([[libraries[0], libraries[1], average_overlap/(len(subset)/2)]], columns=["DGE Method", "DGE Method 2", "Average overlap"])], axis=0)

final_dataframe.to_excel("../data/overlap_analysis_dge_method.xlsx", index=False)