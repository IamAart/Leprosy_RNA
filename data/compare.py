import pandas as pd

data_new = pd.read_excel("./Predictions_current/analysis_rf_predictions.xlsx")
data_old = pd.read_excel("./Predictions_old/analysis_rf_predictions.xlsx")

data_new["Selection"] = data_new["Analysis type"] + "-" + data_new["Feature Selection Method"] + "-" + data_new["Number of genes"].astype(str) + "-" + data_new["Library"].astype(str)
data_old["Selection"] = data_old["Analysis type"] + "-" + data_old["Feature Selection Method"] + "-" + data_old["Number of genes"].astype(str) + "-" + data_old["Library"].astype(str)


combined_data = data_new.merge(data_old, how='left', on='Selection')

intersects = []
news = []
olds = []
intersect_coding = []
new_coding = []
old_coding = []
intersect_non = []
new_non = []
old_non = []

for selection in data_new["Selection"].tolist():
    row_1 = data_new[data_new["Selection"] == selection]
    row_2 = data_old[data_old["Selection"] == selection]
    intersect = set(row_1["Gene names"].values[0].split(", ")) & set(row_2["Gene names"].values[0].split(", "))
    new = set(row_1["Gene names"].values[0].split(", ")) - set(row_2["Gene names"].values[0].split(", "))
    old = set(row_2["Gene names"].values[0].split(", ")) - set(row_1["Gene names"].values[0].split(", "))
    intersects.append(intersect)
    news.append(new)
    olds.append(old)
    
    if "Only combined" in selection:
        if selection.startswith("CODING"):
            intersect_coding.extend(intersect)
            new_coding.extend(row_1["Gene names"].values[0].split(", "))
            old_coding.extend(row_2["Gene names"].values[0].split(", "))
        elif selection.startswith("NON_CODING"):
            intersect_non.extend(intersect)
            new_non.extend(row_1["Gene names"].values[0].split(", "))
            old_non.extend(row_2["Gene names"].values[0].split(", "))

try:
    # create a dataset for the combined dataset, this contains all data and no calculation for averages
    dataset = pd.DataFrame(
        zip(data_new["Selection"].tolist(), intersects, news, olds),
        columns=["selection", "intersect", "new analysis", "old analysis"]
    )
    dataset.to_excel("compare_influence_new_analysis.xlsx")
except Exception as e:
    raise Exception(e)

print("NEW in coding genes", set(new_coding) - set(old_coding), len(set(new_coding) - set(old_coding)))
print("OLD in coding genes",set(old_coding) - set(new_coding), len(set(old_coding) - set(new_coding)))
print("NEW in non-coding genes",set(new_non) - set(old_non), len(set(new_non) - set(old_non)))
print("OLD in non-coding genes",set(old_non) - set(new_non), len(set(old_non) - set(new_non)))
