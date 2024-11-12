library("readxl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("viridis")

rfe_vs_chi2_plot <- function(data) {
    data$feature_selection <- as.factor(data$feature_selection)
    plot <- ggplot(data) +
        geom_boxplot(aes(x=feature_selection, y=auc, color=feature_selection)) +
        labs(title = "Difference in auc score between feature selection methods; RFE and Chi2", x = "Feature Selection Method", y = "Auc Score", color="Feature Selection Method") + 
        scale_color_viridis(
            discrete = TRUE,
            labels = c(rfe = "Recursive Feature Elimination.", chi2 = "Chi-squared Feature Elimination")
        ) +
        theme(legend.position = "right")
    
    ggsave("./data/RF_analysis_plots/boxplot_rfe_vs_chi2.png", plot, width=7, height=7, dpi=300)
}

type_comparison_plot <- function(data) {
    data <- data[data$feature_selection == "rfe", ]
    data$gene_type[data$gene_type == "nc"] <- "non_coding" 
    plot <- ggplot(data) +
        geom_boxplot(aes(x=gene_type, y=auc, color=gene_type)) +
        labs(title = "Difference in auc score between different biotype", x = "Biotype", y = "Auc Score", color="Biotype Information") + 
        scale_color_viridis(
            discrete = TRUE,
            labels = c(all = "All Genes", coding = "Coding genes", non_coding = "Non-coding genes")
        ) +
        theme(legend.position = "right")
    
    ggsave("./data//RF_analysis_plots/boxplot_non_coding_vs_coding.png", plot, width=7, height=7, dpi=300) 
}

box_plot_libraries <- function(data, biotype) {
    data <- data[data$feature_selection == "rfe", ]
    data <- data[data["gene_type"] == biotype, ]
    if (biotype == "nc") {
        data$gene_type[data$gene_type == "nc"] <- "non_coding" 
    }
    plot <- ggplot(data) +
        geom_boxplot(aes(x=Library, y=auc, color=Library)) +
        labs(title = "Difference in auc score between RNA sequencing methods for non-coding genes" , x = "RNA Sequencing Method", y = "Auc Score", color="RNA Sequencing Method") + 
        scale_color_viridis(
            discrete = TRUE,
            labels = c("All" = "All unique genes of each library added together", "DESeq2" = "DESeq2 genes", "EdgeR" = "EdgeR genes", "LimmaVoom" = "LimmaVoom genes", "Only combined libraries" = "Intersect genes of the three libraries")
        ) +
        theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
    
    ggsave("./data/RF_analysis_plots/boxplot_libraries_non_coding_genes.png", plot, width=12, height=6, dpi=300) 
}

data <- read_excel("./data/Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx")
rfe_vs_chi2_plot(data)
type_comparison_plot(data)
# box_plot_libraries(data, "all")
# box_plot_libraries(data, "coding")
box_plot_libraries(data, "nc")
