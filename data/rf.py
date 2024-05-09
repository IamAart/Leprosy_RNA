import itertools
import numpy as np
import pandas as pd
import conorm as co
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split


def train_test_total_dataset():
    counts = pd.read_csv("./sasc326_counts.csv", index_col=0)
    samples = pd.read_csv("./sasc326_samples.csv")
    counts_normalized = co.tmm(counts).transpose()

    dict_samples = dict()
    for i, row in samples.iterrows():
        if row["group"] == "HHC":
            dict_samples[row["sample"]] = 0
        elif row["group"] == "First":
            dict_samples[row["sample"]] = 1

    samples_result = pd.Series(dict_samples)
    full_dataset = counts_normalized.merge(samples_result.to_frame(), left_index=True, right_index=True).rename(columns={0: "Result"})
    return full_dataset


if __name__ == "__main__":
    MIN = 10
    MAX = 40
    dataset = train_test_total_dataset()
    dataset = dataset.dropna(subset=["Result"])
    x = dataset.drop(['Result'], axis=1)
    y = dataset['Result'].astype(int)
    genes = pd.read_csv("./Venn_Diagrams/Ensembl_ID_Data/All_GENES_VENN_DATA.csv")["DESeq2/LimmaVoom/EdgeR"].dropna().to_list()
    x = x[genes]
    list_selectors = list()
    total_selectors = list()
    ratio = 0.8
    for i in range(MIN, MAX+1):
        print(i)
        rfe = RFE(estimator=RandomForestClassifier(random_state=42), step=1, n_features_to_select=i)
        selector = rfe.fit(x, y)
        list_selectors.append(selector.get_feature_names_out(genes))
        total_selectors.extend(selector.get_feature_names_out(genes))
    for selector in list_selectors:
        print(selector)
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=(1 - ratio), random_state=42)
        rf = RandomForestClassifier()
        rf.fit(x_train[selector], y_train)
        y_pred = rf.predict(x_test[selector])
        accuracy = accuracy_score(y_test, y_pred)
        print(len(selector), accuracy, selector)

    all_selectors = set(total_selectors)
    Gene_selection_paper = [
        "ENSG00000204387", "ENSG00000198886", "ENSG00000198786", "ENSG00000198763", "ENSG00000283633",
        "ENSG00000198804", "ENSG00000135090", "ENSG00000135597", "ENSG00000198727", "ENSG00000141933",
        "ENSG00000138722", "ENSG00000150991", "ENSG00000248527", "ENSG00000266538", "ENSG00000006015",
        "ENSG00000175602", "ENSG00000225864", "ENSG00000200183", "ENSG00000279227"
    ]
    for i in all_selectors:
        if i not in Gene_selection_paper:
            print(f"{i} not in paper")
        else:
            print(f"{i} in paper")

    # all_combinations = list()
    # for i in range(MIN,MAX+1):
    #     all_combinations.append(set(itertools.combinations(all_selectors, i)))
    #
    # print(all_combinations)
    # for size in all_combinations:
    #     best_accuracy = -np.Inf
    #     best_selector = None
    #     for combination in size:
    #         selector = [i for i in combination]
    #         x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=(1 - ratio), random_state=42)
    #         rf = RandomForestClassifier()
    #
    #         rf.fit(x_train[selector], y_train)
    #         y_pred = rf.predict(x_test[selector])
    #         accuracy = accuracy_score(y_test, y_pred)
    #         if accuracy >= best_accuracy:
    #             best_selector = selector
    #             best_accuracy = accuracy
    #     print(len(best_selector), best_accuracy, best_selector)
