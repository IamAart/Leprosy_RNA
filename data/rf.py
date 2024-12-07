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
from dotenv import load_dotenv

# Create Global Variables
load_dotenv()
LIBRARY_COMBINATIONS = eval(os.getenv("LIBRARY_COMBINATIONS").replace('"', ''))
MIN = int(os.getenv("MIN"))
MAX = int(os.getenv("MAX"))
ITERATIONS_CROSS_VALIDATION = int(os.getenv("ITERATIONS_CROSS_VALIDATION"))


def save_auc_to_plot(data: dict, path_pic: str):
    """
    Create an AUC plot from the data of false positive rate and true positive rate and save it to a path
    :param data: data from 5 random forest iteration containing the fpr, tpr and auc score
    :param path_pic: the path to save the plot to
    :return: None
    """
    image = plt.figure()
    for key, value in data.items():
        plt.plot(value["fpr"], value["tpr"], label=f"Iteration {key+1}: AUC score = {value['auc']}")
    plt.legend(loc='lower right')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    image.savefig(path_pic)
    plt.close(image)


def train_test_total_dataset():
    """
    Load the dataset, normalize the data with TMM and create a result column that indicates no disease (0) or disease (1).
    After this create a pandas dataframe which is transposed and has a result column.
    :return: return the pandas dataframe with the normalized data and result column
    """
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
    """
    Get the gene names from the ensembl IDs from the features table
    :param ensembles: list of ensembl IDs
    :param features: pandas dataframe with the features table
    :return: gene names of the ensembl IDs in the same order as ensembles
    """
    gene_names = list()
    for ensembles in ensembles:
        selected_gene = features[features['feature'] == ensembles]['gene_name'].values[0]
        gene_names.append(selected_gene)
    return gene_names


def generate_genes(library_combinations: dict, method):
    """
    create the genes from the venn data based on the library combinations and method (All_GENES, CODING, NON_CODING)
    :param library_combinations: dictionary with the library combinations
    :param method: a string with the method (All_GENES, CODING, NON_CODING)
    :return: a dictionary with the genes for each library combination with the library combination as key
    """
    venn_data = pd.read_csv(f"./Venn_Diagrams/Ensemble_data/ENSEMBLE_{method}_VENN_DATA.csv")
    genes = dict()
    for key, columns in library_combinations.items():
        for _, column in enumerate(venn_data):
            if column in columns:
                if key in genes.keys():
                    addition_list = [x for x in venn_data[column].tolist() if str(x) != 'nan']
                    genes[key].extend(addition_list)
                else:
                    genes[key] = [x for x in venn_data[column].tolist() if str(x) != 'nan']

    return genes


def find_best_selector_genes(x_data, y_data, genes, features, min_feat_amount, max_feat_amount, selection_model):
    """

    :param x_data: training data of the gene counts and samples
    :param y_data: result data of each sample
    :param genes: subset of ensemble IDs to be used for selection
    :param features: a pandas dataframe with the features table
    :param min_feat_amount: an integer with the minimum amount of features to select
    :param max_feat_amount: an integer with the maximum amount of features to select
    :param selection_model: a string with the selection model to use (chi2 or RFE)
    :return:
    """
    selectors_ensemble, selectors_gene_names = list(), list()
    for i in range(min_feat_amount, max_feat_amount + 1):
        if selection_model == "chi2":
            best_features_model = SelectKBest(score_func=chi2, k=i)
        elif selection_model == "RFE":
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
            raise KeyError("No selection_model selected or selection_model not chi2 or RFE")

        best_features_model.fit(x_data[genes], y_data)
        selected_features = best_features_model.get_feature_names_out(genes)

        # Save best genes with ensemble and name in different lists
        selectors_ensemble.append(selected_features)
        selectors_gene_names.append(get_gene_names(selected_features.tolist(), features))

    return selectors_ensemble, selectors_gene_names


def train_and_validate_model(x, y, col_names):
    """
    Train and validate the model with Leave One Out Cross Validation
    :param x: all the data containing the gene counts and samples
    :param y: all the data containing the result of each sample
    :param col_names: the column names of the selected features
    :return: the accuracy, AUC score, false positive rate and true positive rate of the trained model
    """
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
    """
    Find the best parameters for the random forest model
    :param x: all the data containing the gene counts and samples
    :param y: all the data containing the result of each sample
    :param col_names: the column names of the selected features
    :return: None
    """
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


def perform_random_forest(selectors, x, y, iterations, selected_features, picture_path):
    """
    Perform the random forest model with the selected features and save the results in a dictionary and create plots
    :param selectors: gene ensembles to use for the random forest model
    :param x: all the data containing the gene counts and samples
    :param y: all the data containing the result of each sample
    :param iterations: the amount of iterations to perform the random forest model
    :param selected_features: the gene names of the selected features (selectors)
    :param picture_path: the path to save the plots to
    :return: a dictionary with the results of the random forest model
    """
    results_dict = dict()
    for idx, selector in enumerate(selectors):
        print(idx)
        selector_dict = dict()
        plot_data = dict()
        for i in range(iterations):
            accuracy, auc_score, fpr, tpr = train_and_validate_model(x, y, selector)
            plot_data[i] = {"fpr": fpr, "tpr": tpr, "auc": auc_score}

            # Save data in dict
            selector_dict[i] = {"accuracy": accuracy, "roc_auc": auc_score.item()}

        # Save figure with AUC plot
        path = f"{picture_path}/{idx + MIN}.png"
        save_auc_to_plot(plot_data, path)

        results_dict[idx+MIN] = {
                "ensembles": ', '.join(selector),
                "genes": ', '.join(selected_features[idx]),
                "results": selector_dict,
                "plot_path": path
        }
        print(results_dict[idx+MIN])
    return results_dict


if __name__ == "__main__":
    features = pd.read_csv("./sasc326_features.csv")
    # Load and adjust data to correct format
    dataset = train_test_total_dataset()
    dataset = dataset.drop(["s103830.003.011", "s103830.004.017"])
    dataset = dataset.dropna(subset=["Result"])

    # Create x and y dataset and get gene sets from file
    x_dataset = dataset.drop(['Result'], axis=1)
    y_dataset = dataset['Result'].astype(int)

    subset_model = os.getenv("SUBSET_MODEL")
    # if these are not set, find the results for all the genes of all the biotypes of all the library combinations
    if subset_model is None:
        for method in ["All_GENES", "CODING", "NON_CODING"]:
            all_genes = generate_genes(LIBRARY_COMBINATIONS, method)

            # Can switch below for OWN GENES
            for option, genes in all_genes.items():
                print(option)

                if x_dataset is None or y_dataset is None:
                    raise Exception("Training data is not correct")

                # find_best_parameters(x_dataset, y_dataset, genes) # Only used to find parameters once
                for model in ["RFE", "chi2"]:
                        dir_path = f"./Predictions/{method}/{option}/{model}/Pictures"
                        if not os.path.exists(dir_path):
                            os.makedirs(dir_path)
                        list_selectors, gene_names = find_best_selector_genes(x_dataset, y_dataset, genes, features, MIN, MAX, model)
                        result_dict = perform_random_forest(list_selectors, x_dataset, y_dataset, ITERATIONS_CROSS_VALIDATION, gene_names, dir_path)
                        with open(f"./Predictions/{method}/{option}/{model}/results.json", 'w') as f:
                            json.dump(result_dict, f)
    # if these are set, find the results for the subset of genes
    else:
        subset_ensembles = os.getenv("SUBSET_ENSEMBLES").split(",")
        subset_min = int(os.getenv("SUBSET_MIN"))
        subset_max = int(os.getenv("SUBSET_MAX"))
        dir_path = f"./Predictions/Subset/Pictures"
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        list_selectors, gene_names = find_best_selector_genes(x_dataset, y_dataset, subset_ensembles, features, subset_min, subset_max, subset_model)
        result_dict = perform_random_forest(list_selectors, x_dataset, y_dataset, ITERATIONS_CROSS_VALIDATION, gene_names, dir_path)
        with open(f"./Predictions/Subset/results.json", 'w') as f:
            json.dump(result_dict, f)
