import pandas as pd

data = pd.read_excel("data/Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx")
data = data.sort_values("auc", ascending=False)
data = data[data["auc"]>=0.96]
gene_names = dict()
ensembles = dict()
print(data.head)

for _, row in data.iterrows():
    e_list = row["ensembles"].split(", ")
    g_list = row["Gene names"].split(", ")
    for i in range(len(e_list)):
        if e_list[i] not in ensembles.keys():
            ensembles[e_list[i]] = 1
            gene_names[g_list[i]] = 1
        else: 
            ensembles[e_list[i]] += 1
            gene_names[g_list[i]] += 1

ensemble_list = sorted(ensembles, key=lambda k: ensembles[k], reverse=True)[0:20]
limmavoom = pd.read_excel("data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx")
results = limmavoom[limmavoom["Row.names"].isin(ensemble_list)]
sorted_values = results.sort_values("P.Adjust", ascending=True)
sorted_values = sorted_values[["gene_name", "Row.names", "biotype"]]
print(sorted_values)