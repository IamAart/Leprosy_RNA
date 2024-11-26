import pandas as pd

data = pd.read_excel("data/Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx")
data = data.sort_values("auc", ascending=False)
gene_names = dict()
ensembles = dict()
print(data.head)

index = 1
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
    
    if index == 25:
        break
    index += 1
# print(ensembles)
print({k: v for k, v in sorted(gene_names.items(), key=lambda item: item[1], reverse=True)})
print(sorted(ensembles, key=lambda k: ensembles[k], reverse=True)[0:22])
print(sorted(gene_names, key=lambda k: gene_names[k], reverse=True)[0:22])