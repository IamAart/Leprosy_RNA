import os

import pandas as pd
import conorm as co
import json
import numpy as np
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, chi2, RFE
from sklearn.metrics import accuracy_score, roc_curve, auc
from sklearn.model_selection import LeaveOneOut, RandomizedSearchCV


# Create Global Variables
MIN = 3
MAX = 20
ratio = 0.8
ITERATIONS_CROSS_VALIDATION = 5
LIBRARY_COMBINATIONS = {
        "Only combined libraries": ["DESeq2/LimmaVoom/EdgeR"],
        "All": ["DESeq2/LimmaVoom/EdgeR", "DESeq2/EdgeR", "DESeq2/LimmaVoom", "DESeq2", "EdgeR/LimmaVoom", "EdgeR",
                "LimmaVoom"],
        # "DESeq2": ["DESeq2/LimmaVoom/EdgeR", "DESeq2/EdgeR", "DESeq2/LimmaVoom", "DESeq2"],
        # "EdgeR": ["DESeq2/LimmaVoom/EdgeR", "DESeq2/EdgeR", "EdgeR/LimmaVoom", "EdgeR"],
        # "LimmaVoom": ["DESeq2/LimmaVoom/EdgeR", "LimmaVoom/EdgeR", "DESeq2/LimmaVoom", "LimmaVoom"]
    }


def save_auc_to_plot(data: dict, path_pic: str):
    image = plt.figure()
    for key, value in data.items():
        plt.plot(value["fpr"], value["tpr"], label=f"Iteration {key+1}: AUC score = {value['auc']}")
    plt.legend(loc='lower right')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    image.savefig(path_pic)
    plt.close(image)


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
    full_dataset = (counts_normalized
                    .merge(samples_result.to_frame(), left_index=True, right_index=True)
                    .rename(columns={0: "Result"}))

    return full_dataset


def get_gene_names(ensembles, features):
    selected_genes = features[features['feature'].isin(ensembles)]
    gene_names = selected_genes["gene_name"].tolist()
    return gene_names


def compare_with_paper(ensembles):
    all_ensembles = set(ensembles)
    Gene_selection_paper = [
        "ENSG00000204387", "ENSG00000198886", "ENSG00000198786", "ENSG00000198763", "ENSG00000283633",
        "ENSG00000198804", "ENSG00000135090", "ENSG00000135597", "ENSG00000198727", "ENSG00000141933",
        "ENSG00000138722", "ENSG00000150991", "ENSG00000248527", "ENSG00000266538", "ENSG00000006015",
        "ENSG00000175602", "ENSG00000225864", "ENSG00000200183", "ENSG00000279227"
    ]
    for gene in all_ensembles:
        if gene not in Gene_selection_paper:
            print(f"{gene} not in paper")
        else:
            print(f"{gene} in paper")


def generate_genes(library_combinations: dict, method):
    venn_data = pd.read_csv(f"./Venn_Diagrams/Ensemble_data/ENSEMBLE_{method}_VENN_DATA.csv")
    # venn_data = pd.read_csv("./Venn_Diagrams/Ensemble_data/ENSEMBLE_NON_CODING_VENN_DATA.csv")
    genes = dict()
    for key, columns in library_combinations.items():
        for id, column in enumerate(venn_data):
            if column in columns:
                if key in genes.keys():
                    addition_list = [x for x in venn_data[column].tolist() if str(x) != 'nan']
                    genes[key].extend(addition_list)
                else:
                    genes[key] = [x for x in venn_data[column].tolist() if str(x) != 'nan']

    return genes


def find_best_selector_genes(x_data, y_data, genes, features, min_feat_amount, max_feat_amount, model):
    selectors_ensemble, selectors_gene_names = list(), list()

    for i in range(min_feat_amount, max_feat_amount + 1):
        # Find the best genes based on Chi squared test for amount of genes

        if model == "chi2":
            best_features_model = SelectKBest(score_func=chi2, k=i)
        elif model == "RFE":
            rf = RandomForestClassifier(
                n_estimators=100,
                max_depth=25,
                min_samples_split=2,
                min_samples_leaf=2,
                max_features="sqrt",
                bootstrap=True
            )
            best_features_model = RFE(rf, n_features_to_select=i)
        else:
            raise KeyError("No model selected")

        best_features_model.fit(x_data[genes], y_data)
        selected_features = best_features_model.get_feature_names_out(genes)

        # Save best genes with ensemble and name in different lists
        selectors_ensemble.append(selected_features)
        print(get_gene_names(selected_features.tolist(), features))
        selectors_gene_names.append(get_gene_names(selected_features.tolist(), features))

    return selectors_ensemble, selectors_gene_names


def train_and_validate_model(x, y, col_names):
    prediction_list, prediction_probability_list, y_test_list = list(), list(), list()

    # Perform Leave One Out Cross Validation
    loo = LeaveOneOut()
    for train_id, test_id in loo.split(x):
        # initiate random forest model
        rf = RandomForestClassifier(
            n_estimators=100,
            max_depth=25,
            min_samples_split=2,
            min_samples_leaf=2,
            max_features="sqrt",
            bootstrap=True
        )

        x_train_cv, x_test_cv = x[col_names].iloc[train_id], x[col_names].iloc[test_id]
        y_train_cv, y_test_cv = y.iloc[train_id], y.iloc[test_id]

        # train random forest model
        rf.fit(x_train_cv, y_train_cv)

        # Save true y value, prediction y value and prediction probability of true value (person gets sick)
        y_test_list.append(y_test_cv)
        prediction_list.append(rf.predict(x_test_cv))
        prediction_probability_list.append(rf.predict_proba(x_test_cv)[:, 1])

    # Calculate Accuracy, FPR, TPR and AUC
    accuracy = accuracy_score(y_test_list, prediction_list)
    fpr_cv, tpr_cv, _ = roc_curve(y_test_list, prediction_probability_list)
    y_roc_auc_cv = auc(fpr_cv, tpr_cv)

    return accuracy, y_roc_auc_cv, fpr_cv, tpr_cv


def find_best_parameters(x, y, col_names):
    rf = RandomForestClassifier()
    # Create the random grid
    random_grid = {'n_estimators': [int(x) for x in np.linspace(start=50, stop=1000, num=10)], # Number of trees in random forest
                   'max_features': ['auto', 'sqrt', None], # Number of features to consider at every split
                   'max_depth': [int(x) for x in np.linspace(5, 100, num=5)].append(None), # Maximum number of levels in tree
                   'min_samples_split': [int(x) for x in np.linspace(start=2, stop=10, num=1)], # Minimum number of samples required to split a node
                   'min_samples_leaf': [1, 2, 4], # Minimum number of samples required at each leaf node
                   'bootstrap': [True, False]} # Method of selecting samples for training each tree
    # Generate a random search cross validation model to find best parameters
    random_search = RandomizedSearchCV(estimator=rf, param_distributions = random_grid, n_iter=100)
    random_search.fit(x[col_names], y)
    print(random_search.best_params_)


def perform_random_forest(selectors, x, y, iterations, gene_names, picture_path):
    results_dict = dict()
    for idx, selector in enumerate(selectors):
        selector_dict = dict()
        plot_data = dict()
        for i in range(iterations):
            accuracy, auc_score, fpr, tpr = train_and_validate_model(x, y, selector)
            plot_data[i] = {"fpr": fpr, "tpr": tpr, "auc": auc_score}

            # Save data in dict
            selector_dict[i] = {"accuracy": accuracy, "roc_auc": auc_score}

        # Save figure with AUC plot
        path = f"{picture_path}/{idx + MIN}.png"
        save_auc_to_plot(plot_data, path)

        results_dict[idx+MIN] = {
                "ensembles": ', '.join(selector),
                "genes": ', '.join(gene_names[idx]),
                "results": selector_dict,
                "plot_path": path
        }
        print(results_dict[idx+MIN])
    return results_dict


if __name__ == "__main__":
    features = pd.read_csv("./sasc326_features.csv")
    # Load and adjust data to correct format
    dataset = train_test_total_dataset()
    bad_samples = dataset.drop(["s103830.003.011", "s103830.004.017"])
    dataset = dataset.dropna(subset=["Result"])

    # Create x and y dataset and get gene sets from file
    x_dataset = dataset.drop(['Result'], axis=1)
    y_dataset = dataset['Result'].astype(int)

    for method in ["All_GENES", "NON_CODING"]:
        all_genes = generate_genes(LIBRARY_COMBINATIONS, method)

        for option, genes in all_genes.items():
            print(option)
            # find_best_parameters(x_dataset, y_dataset, genes) # Only used to find parameters once

            for model in ["RFE", "chi2"]:
                dir_path = f"./Predictions/{method}/{option}/{model}/Pictures"
                if not os.path.exists(dir_path):
                    os.makedirs(dir_path)
                list_selectors, gene_names = find_best_selector_genes(x_dataset, y_dataset, genes, features, MIN, MAX, model)
                result_dict = perform_random_forest(list_selectors, x_dataset, y_dataset, ITERATIONS_CROSS_VALIDATION, gene_names, dir_path)
                with open(f"./Predictions/{method}/{option}/{model}/results.json", 'w') as f:
                    json.dump(result_dict, f)
