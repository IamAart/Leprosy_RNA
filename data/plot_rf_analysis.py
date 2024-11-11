import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def adjust_data(data, nc_adjust):
    if nc_adjust:
        nc_data = data[data["Analysis type"] == "NON_CODING"]
        nc_data_10 = nc_data[nc_data["Number of genes"].astype(int) <= 10]
        all_data = data[data["Analysis type"] == "All_GENES"]
        data = pd.concat([all_data, nc_data_10])
        return data
    return data


def box_plot_rfe_vs_chi2(data):
    plot = plt.figure()
    sns.boxplot(
        data=data,
        x="Number of genes",
        y="Average AUC score",
        hue="Feature Selection Method",
        width=0.8,
        showmeans=True,
        fliersize=1,
    )
    plot.savefig("./RF_analysis_plots/boxplot_rfe_vs_chi2_with_number_of_genes.png")
    plt.close(plot)
    plot2 = plt.figure()
    sns.boxplot(
        data=data,
        x="Feature Selection Method",
        y="Average AUC score",
        width=0.4,
        showmeans=True,
        fliersize=1,
    )
    plot2.savefig("./RF_analysis_plots/boxplot_rfe_vs_chi2.png")
    plt.close(plot2)


def box_plot_type_comparison(data):
    data = data[data["Feature Selection Method"] == "RFE"]
    plot2 = plt.figure()
    sns.boxplot(
        data=data,
        x="Analysis type",
        y="Average AUC score",
        width=0.4,
        showmeans=True,
        fliersize=1,
    )
    plot2.savefig("./RF_analysis_plots/boxplot_non_coding_vs_coding.png")
    plt.close(plot2)


def box_plot_libraries(data):
    data = data[data["Feature Selection Method"] == "RFE"]
    data1 = data[data["Analysis type"] == "All_GENES"]
    plot = plt.figure()
    sns.boxplot(
        x=data1["Library"],
        y=data1["Average AUC score"],
        width=0.8,
        showmeans=True,
        fliersize=1,
        palette="tab10"
    )
    plot.savefig("./RF_analysis_plots/boxplot_libraries_all_genes.png")
    plt.close(plot)
    data2 = data[data["Analysis type"] == "NON_CODING"]
    plot2 = plt.figure()
    sns.boxplot(
        x=data2["Library"],
        y=data2["Average AUC score"],
        width=0.8,
        showmeans=True,
        fliersize=1,
        palette="tab10"
    )
    plot2.savefig("./RF_analysis_plots/boxplot_libraries_non_coding.png")
    plt.close(plot2)
    data3 = data[data["Analysis type"] == "CODING"]
    plot3 = plt.figure()
    sns.boxplot(
        x=data3["Library"],
        y=data3["Average AUC score"],
        width=0.8,
        showmeans=True,
        fliersize=1,
        palette="tab10"
    )
    plot3.savefig("./RF_analysis_plots/boxplot_libraries_coding.png")
    plt.close(plot3)


rf_results = pd.read_excel('./Predictions_current/analysis_rf_predictions.xlsx')
box_plot_rfe_vs_chi2(rf_results)
box_plot_type_comparison(rf_results)
box_plot_libraries(rf_results)