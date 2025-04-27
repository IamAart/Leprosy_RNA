import json
import os

import conorm as co
import numpy as np
import pandas as pd
from dotenv import load_dotenv
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap
from scipy.stats import uniform
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, chi2, RFE
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_curve, auc, precision_score, recall_score, f1_score
from sklearn.model_selection import LeaveOneOut, RandomizedSearchCV, GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

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
    counts = pd.read_csv("../data/sasc326_counts.csv", index_col=0)
    samples = pd.read_csv("../data/sasc326_samples.csv")
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


def generate_genes(library_combinations: dict, method: str, model: str):
    """
    create the genes from the venn data based on the library combinations and method (All_GENES, CODING, NON_CODING)
    :param library_combinations: dictionary with the library combinations
    :param method: a string with the method (All_GENES, CODING, NON_CODING)
    :return: a dictionary with the genes for each library combination with the library combination as key
    """
    venn_data = pd.read_csv(f"../data/Venn_Diagrams/Ensemble_data/ENSEMBLE_{method}_VENN_DATA.csv")
    genes = dict()
    if model != "RF":
        library_combinations = {k: v for k, v in library_combinations.items() if k == "Intersection"}
    for key, columns in library_combinations.items():
        for _, column in enumerate(venn_data):
            if column in columns:
                if key in genes.keys():
                    addition_list = [x for x in venn_data[column].tolist() if str(x) != 'nan']
                    genes[key].extend(addition_list)
                else:
                    genes[key] = [x for x in venn_data[column].tolist() if str(x) != 'nan']

    return genes


def select_ml_model(model: str, parameters: dict = None):
    if model == "RF":
        if not parameters:
            parameters = {'n_estimators': 100, 'max_depth': 25, 'min_samples_split': 2, 'min_samples_leaf': 2,
                          'max_features': 'sqrt', 'bootstrap': True}
        return RandomForestClassifier(**parameters)
    elif model == "SVM":
        if not parameters:
            parameters = {'kernel': 'linear', 'shrinking': True, 'C': 3.8, 'gamma': 0.02, 'break_ties': True,
                          'probability': False, 'tol': 0.0004}
        else:
            parameters.update({'kernel': 'linear', 'shrinking': True, 'break_ties': True,
                          'probability': False})
        return SVC(**parameters)
    elif model == "LLR":
        if not parameters:
            return LogisticRegression(max_iter=10000)
        else:
            parameters.update({'max_iter': 10000})
        return LogisticRegression(**parameters)
    elif model == "kNN":
        if not parameters:
            parameters = {'algorithm': 'auto', 'leaf_size': 20, 'metric': 'euclidean', 'n_neighbors': 9, 'p': 1, 'weights': 'uniform'}
        return KNeighborsClassifier(**parameters)
    else:
        raise KeyError("No model selected or model is not equal to RF, SVM, LLR or KNN")


def find_optimal_parameters(X, y, model):
    if model == 'SVM':
        # takes longer than an hour, so run once and then use that
        # clf = SVC(kernel="linear", decision_function_shape='ovr', shrinking=True, max_iter=-1, verbose=0, probability=False)
        # params = {'C': uniform(0.1, 10),
        #           'gamma': uniform(0.001, 0.1),
        #           'tol': uniform(0.0001, 0.0005),
        #           }
        # search = RandomizedSearchCV(
        #     clf,
        #     param_distributions=params,
        #     n_iter=100,
        #     cv=10,
        #     random_state=42,
        #     verbose=0
        # )
        return {'kernel': 'linear', 'shrinking': True, 'C': 3.8, 'gamma': 0.02, 'break_ties': True,
                          'probability': False, 'tol': 0.0004}
    elif model == 'RF':
        clf = RandomForestClassifier()
        params = {
            'n_estimators': [int(n) for n in np.linspace(start=50, stop=200, num=10)],
            'max_features': ['sqrt', 'log2'],
            'max_depth': [int(m) for m in np.linspace(10, 110, num=11)] + [None],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 4],
            'bootstrap': [True, False]
        }
        search = RandomizedSearchCV(
            estimator=clf,
            param_distributions=params,
            n_iter=100, cv=5, verbose=0, random_state=42, n_jobs=-1
        )
    elif model == 'LLR':
        clf = LogisticRegression(max_iter=10000)
        params = {
            'C': [0.01, 0.1, 1, 10, 100],
            'penalty': ['l1', 'l2', 'elasticnet', 'none'],  # Regularization type
            'solver': ['liblinear', 'saga'],  # Solver choice depends on penalty
            'l1_ratio': [0.3, 0.5, 0.7]
        }
        search = GridSearchCV(clf, params, cv=5, scoring='accuracy', verbose=0)
    elif model == 'kNN':
        clf = KNeighborsClassifier()
        params = {
            'n_neighbors': [3, 5, 7, 9, 11],
            'weights': ['uniform', 'distance'],
            'algorithm': ['auto', 'ball_tree', 'kd_tree', 'brute'],
            'leaf_size': [20, 30, 40],
            'p': [1, 2],
            'metric': ['euclidean', 'manhattan', 'minkowski', 'cosine']
        }
        search = GridSearchCV(estimator=clf, param_grid=params, cv=5, scoring='accuracy', verbose=0)
    else:
        raise KeyError('No model selected')

    search.fit(X, y)
    best_params = search.best_params_
    return best_params


def find_best_selector_genes(x_data, y_data, genes, features, min_feat_amount, max_feat_amount, selection_model, model):
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
            print(i, min_feat_amount, max_feat_amount)
            classifier = select_ml_model(model)
            best_features_model = RFE(classifier, n_features_to_select=i)
        else:
            raise KeyError("No selection_model selected or selection_model not chi2 or RFE")

        best_features_model.fit(x_data[genes], y_data)
        selected_features = best_features_model.get_feature_names_out(genes)

        # Save best genes with ensemble and name in different lists
        selectors_ensemble.append(selected_features)
        selectors_gene_names.append(get_gene_names(selected_features.tolist(), features))

    return selectors_ensemble, selectors_gene_names


def train_and_validate_model(x, y, col_names, parameters, model):
    """
    Train and validate the model with Leave One Out Cross Validation
    :param model:
    :param parameters:
    :param x: all the data containing the gene counts and samples
    :param y: all the data containing the result of each sample
    :param col_names: the column names of the selected features
    :return: the accuracy, AUC score, false positive rate and true positive rate of the trained model
    """
    prediction_list, prediction_probability_list, y_test_list, weights, intercepts = list(), list(), list(), list(), list()

    # Perform Leave One Out Cross Validation
    loo = LeaveOneOut()
    for train_id, test_id in loo.split(x):
        # create train and test sets
        x_train_cv, x_test_cv = x[col_names].iloc[train_id], x[col_names].iloc[test_id]
        y_train_cv, y_test_cv = y.iloc[train_id], y.iloc[test_id]

        # initialize machine learning model
        classifier = select_ml_model(model, parameters)
        # train random forest model
        classifier.fit(x_train_cv, y_train_cv)

        # Save true y value, prediction y value and prediction probability of true value (person gets sick)
        y_test_list.append(y_test_cv)
        prediction_list.append(classifier.predict(x_test_cv))
        if model != "SVM":
            prediction_probability_list.append(classifier.predict_proba(x_test_cv)[:, 1])
        else:
            prediction_probability_list.append(classifier.decision_function(x_test_cv))

        # calculate feature importances of the selector
        feature_importance = None
        if model == "RF":
            feature_importance = classifier.feature_importances_
        elif model == "SVM" or model == "LLR":
            feature_importance = classifier.coef_[0]
            if model == "SVM":
                intercepts.append(classifier.intercept_[0])
        weights.append(feature_importance)

    # Calculate Accuracy, FPR, TPR and AUC
    accuracy = accuracy_score(y_test_list, prediction_list)
    fpr_cv, tpr_cv, _ = roc_curve(y_test_list, prediction_probability_list)
    y_roc_auc_cv = auc(fpr_cv, tpr_cv)
    precision = precision_score(y_test_list, prediction_list)
    recall = recall_score(y_test_list, prediction_list)
    f1 = f1_score(y_test_list, prediction_list)
    specificity = recall_score(y_test_list, prediction_list, pos_label=0)
    avg_intercept = None
    avg_weights = None
    if model != "kNN":
        avg_weights = sum(weights) / len(weights)
        if model == "SVM":
            avg_intercept = sum(intercepts) / len(intercepts)

    return accuracy, y_roc_auc_cv, fpr_cv, tpr_cv, precision, recall, f1, specificity, avg_weights, avg_intercept


def get_analysis_id(data, model, gene_type, dge_method, feature_selection):
    result = data[(data["model"]==model) & (data["gene_type"]==gene_type) & (data["DGE method"]==dge_method) & (data["feature_selection"]==feature_selection)]
    return result["analysis_id"].values[0]


def perform_prediction(selectors, x, y, iterations, selected_features, picture_path, model, gene_type, dge_method, feature_selection):
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
    analyses_data = pd.read_excel("../data/analysis_ids_machine_learning.xlsx")
    first_intercept = None
    for idx, selector in enumerate(selectors):
        if gene_type == "Subset":
            analysis_id = "Subset"
        else:
            analysis_id = get_analysis_id(analyses_data, model, gene_type, dge_method, feature_selection)
        print(analysis_id)

        selector_dict = dict()
        plot_data = dict()
        if model == "SVM": iterations = 1
        weights, intercept = None, None

        parameters = find_optimal_parameters(x[selector], y, model)

        for i in range(iterations):
            accuracy, auc_score, fpr, tpr, precision, recall, f1, specificity, weights, intercept = train_and_validate_model(
                x, y, selector, parameters, model
            )
            plot_data[i] = {"fpr": fpr, "tpr": tpr, "auc": auc_score}

            # Save data in dict
            selector_dict[i] = {"accuracy": accuracy, "roc_auc": auc_score.item(), "precision": precision.item(),
                                "recall": recall.item(), "f1": f1.item(), "specificity": specificity.item()}

        if weights is not None:
            gene_importance_df = pd.DataFrame(list(zip(selector, selected_features[idx], weights, [abs(x) for x in weights])),
                                              columns=["selector", "gene_names", "weights", "absolute weights"])
            sorted_gene_importance_df = gene_importance_df.sort_values(by="absolute weights", ascending=False)
            weight_path = picture_path.replace("Pictures", "")
            if not os.path.exists(picture_path):
                os.makedirs(weight_path)
            sorted_gene_importance_df.to_csv(weight_path + f"{idx + MIN}_weights.csv")


        if intercept is not None and idx + MIN == 3 and model == "SVM":
            first_intercept = intercept
        # Save figure with AUC plot
        path = f"{picture_path}/{idx + MIN}.png"
        save_auc_to_plot(plot_data, path)

        results_dict[idx+MIN] = {
            "analysis_id": analysis_id,
            "timestamp": pd.Timestamp.now().strftime("%d-%m-%Y %H:%M:%S"),
            "model": model,
            "gene_type": gene_type,
            "dge_method": dge_method,
            "feature_selection": feature_selection,
            "group_size": idx + MIN,
            "results": selector_dict,
            "ensembl": ', '.join(selector),
            "gene_names": ', '.join(selected_features[idx]),
            "hyperparameters": parameters,
            "plot_path": path
        }
    return results_dict, first_intercept


def hyperplane_creation(X, y, intercept, method, feature_selection, path):
    feature_importance = pd.read_csv(path + "/3_weights.csv")
    gene1 = feature_importance.iloc[0]
    gene2 = feature_importance.iloc[1]
    gene3 = feature_importance.iloc[2]

    xx, yy = np.meshgrid(np.linspace(X[gene1["selector"]].min(), X[gene1["selector"]].max(), 50),
                         np.linspace(X[gene2["selector"]].min(), X[gene2["selector"]].max(), 50))
    zz = (- gene1["weights"] * xx - gene2["weights"] * yy - intercept) / gene3["weights"]

    cmap = ListedColormap(['#882255', '#44aa99'])

    # Initialize the figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot and decision boundary
    ax.scatter(X[gene1["selector"]], X[gene2["selector"]], X[gene3["selector"]], c=y, cmap=cmap, s=30)
    ax.plot_surface(xx, yy, zz, alpha=0.2)

    # Animation function
    def update(frame):
        ax.view_init(elev=0, azim=frame)
        ax.set_xlabel(gene1["gene_names"], labelpad=15)
        ax.set_ylabel(gene2["gene_names"], labelpad=15)
        ax.set_zlabel(gene3["gene_names"], labelpad=15)

        xticks = ax.get_xticks()
        yticks = ax.get_yticks()
        zticks = ax.get_zticks()

        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_zticks(zticks)

        if 45 < frame < 135 or 225 < frame < 315:
            ax.set_yticklabels([""] * len(yticks))
        else:
            ax.set_yticklabels([f"{tick:.0f}" for tick in yticks])

        # Hide X-axis labels when X-axis is facing away
        if 135 < frame < 225 or 315 < frame < 360 or 0 < frame < 45:
            ax.set_xticklabels([""] * len(xticks))
        else:
            ax.set_xticklabels([f"{tick:.0f}" for tick in xticks])

        # Always keep Z-axis labels visible
        ax.set_zticklabels([f"{tick:.0f}" for tick in zticks])

    # Create animation
    anim = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), interval=50, blit=False)

    # Save as a GIF
    output_path = f"{path}/hyperplane.gif"
    anim.save(output_path, writer='pillow', fps=5, dpi=720)
    plt.close(fig)


if __name__ == "__main__":
    features = pd.read_csv("../data/sasc326_features.csv")
    # Load and adjust data to correct format
    dataset = train_test_total_dataset()
    dataset = dataset.drop(["s103830.003.011", "s103830.004.017"])
    dataset = dataset.dropna(subset=["Result"])

    # Create x and y dataset and get gene sets from file
    x_dataset = dataset.drop(['Result'], axis=1)
    y_dataset = dataset['Result'].astype(int)
    print("Step 1: loading data")
    subset_model = os.getenv("SUBSET_MODEL")
    # if these are not set, find the results for all the genes of all the biotypes of all the library combinations
    if subset_model is None:
        for ml_model in ["RF", "SVM", "LLR", "kNN"]:
            for method in ["ALL", "NC", "CODING"]:
                all_genes = generate_genes(LIBRARY_COMBINATIONS, method, ml_model) # INTERSECTION of alle dingen
                if ml_model != "RF":
                    all_genes = {k: v for k, v in all_genes.items() if k == "Intersection"}
                # Can switch below for OWN GENES
                for option, genes in all_genes.items():
                    if x_dataset is None or y_dataset is None:
                        raise Exception("Training data is not correct")

                    # find_best_parameters(x_dataset, y_dataset, genes) # Only used to find parameters once
                    for model in ["RFE", "chi2"]:
                        if ml_model == "kNN" and model == "RFE": continue
                        
                        dir_path = f"../data/Predictions/{ml_model}/{method}/{option}/{model}"
                        if os.path.exists(dir_path):
                            dir_path = f"../data/Predictions_new/{ml_model}/{method}/{option}/{model}"
                        os.makedirs(f"{dir_path}/Pictures")
                        list_selectors, gene_names = find_best_selector_genes(x_dataset, y_dataset, genes, features, MIN, MAX, model, ml_model)
                        result_dict, intercept = perform_prediction(list_selectors, x_dataset, y_dataset, ITERATIONS_CROSS_VALIDATION, gene_names, f"{dir_path}/Pictures", ml_model, method, option, model)
                        
                        with open(f"{dir_path}/results.json", 'w') as f:
                            print("Dumped the json at ", f"{dir_path}/results.json")
                            json.dump(result_dict, f)
                        if intercept:
                            hyperplane_creation(x_dataset, y_dataset, intercept, method, model, dir_path)
    # if these are set, find the results for the subset of genes
    else:
        subset_ensembles = os.getenv("SUBSET_ENSEMBLES").split(",")
        subset_min = int(os.getenv("SUBSET_MIN"))
        subset_max = int(os.getenv("SUBSET_MAX"))
        dir_path = f"../data/Predictions/Subset"
        os.makedirs(f"{dir_path}/Pictures")
        list_selectors, gene_names = find_best_selector_genes(x_dataset, y_dataset, subset_ensembles, features, subset_min, subset_max, "RFE", subset_model)
        result_dict, intercept = perform_prediction(list_selectors, x_dataset, y_dataset, ITERATIONS_CROSS_VALIDATION, gene_names, f"{dir_path}/Pictures", "SVM", "Subset", "Subset", "RFE")
        with open(f"{dir_path}/results.json", 'w') as f:
            print("Dumped the json at ", f"{dir_path}/results.json")
            json.dump(result_dict, f)
        if intercept:
            hyperplane_creation(x_dataset, y_dataset, intercept, "Subset", "RFE", dir_path)
