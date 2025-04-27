library("readxl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("edgeR")
library("ggpubr")
library("reshape2")
library("stringr")
library("dotenv")

# .libPaths("./R_Packages/")
dotenv::load_dot_env(file = ".env")
NON_CODING_BIOTYPES <- as.vector(strsplit(Sys.getenv("NON_CODING_BIOTYPES"), ","))[[1]]
CODING_BIOTYPES <- as.vector(strsplit(Sys.getenv("CODING_BIOTYPES"), ","))[[1]]
CUT_OFF <- as.numeric(Sys.getenv("CUT_OFF"))
LOG2_OFF <- log2(as.numeric(Sys.getenv("LOG2_OFF")))

colors <- c("#88ccee", "#44aa99", "#117733", "#332288", "#ddcc77", "#999933", "#cc6677", "#882255", "#aa4499", "#dddddd")

rfe_vs_chi2_plot <- function(data, plot_path) {
    data <- data[data$model == "RF", ]
    data$`feature_selection` <- as.factor(data$`feature_selection`)
    print(unique(data$`feature_selection`))
    means <- data %>% group_by("feature_selection") %>% summarize(m=mean(`auc`))
    print(means)
    plot <- ggplot(data, aes(x=feature_selection, y=`auc`, fill=feature_selection)) +
        geom_boxplot() +
        theme_minimal() +
        labs(x = "feature_selection", y = "auc", color="feature_selection") +
        scale_fill_manual(
            values = c(RFE = colors[7], chi2 = colors[2]),
            labels = c(RFE = "RFE", chi2 = expression(chi^2))
        ) +
        scale_x_discrete(
            labels = c("RFE" = "RFE", "chi2" = expression(chi^2))
        ) +
        scale_y_continuous(
            breaks = c(0.85, 0.9, 0.95, 1.0),
            labels = c("0.85", "0.90", "0.95", "1.00")
        ) +
        theme(
            legend.position = "none",
            # Axis titles font size
            axis.title.x = element_text(size = 19),    # x-axis label size
            axis.title.y = element_text(size = 19),    # y-axis label size

            # Axis tick labels font size
            axis.text.x = element_text(size = 17),     # x-axis ticks (labels) size
            axis.text.y = element_text(size = 17),     # y-axis ticks (labels) size

            # Legend title and text font size
            legend.title = element_text(size = 17),    # Legend title size
            legend.text = element_text(size = 16)      # Legend labels size
        ) +
        stat_compare_means(comparisons = list(c("chi2", "RFE")), size=6) +
        stat_summary(
            fun = "mean", 
            geom = "point",
            shape = 4,
            size = 2,
            stroke = 1.5
        )

    ggsave(plot_path, plot, width=15, height=10, dpi=720)
}

type_comparison_plot <- function(data, plot_path) {
    # data$`gene_type`[data$`gene_type` == "NON_CODING"] <- "NC"
    # data$`gene_type`[data$`gene_type` == "CODING"] <- "CODING"
    # data$`gene_type`[data$`gene_type` == "All_GENES"] <- "ALL"
    data <- data[data$model == "RF", ]
    means <- data %>% group_by(`gene_type`) %>% summarize(m=mean(`auc`))
    print(means)
    plot <- ggplot(data, aes(x=`gene_type`, y=`auc`, fill=`gene_type`)) +
        geom_boxplot() +
        labs(x = "gene_type", y = "auc", fill="gene_type") +
        scale_fill_manual(
            values = c(ALL = colors[1], NC = colors[1], CODING = colors[1]),
        ) +
        scale_x_discrete(
            limits = c("ALL", "NC", "CODING")
        ) +
        scale_y_continuous(
            breaks = c(0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
            labels = c("0.75", "0.80", "0.85", "0.90", "0.95", "1.00")
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            # Axis titles font size
            axis.title.x = element_text(size = 19),    # x-axis label size
            axis.title.y = element_text(size = 19),    # y-axis label size

            # Axis tick labels font size
            axis.text.x = element_text(size = 17),     # x-axis ticks (labels) size
            axis.text.y = element_text(size = 17),     # y-axis ticks (labels) size

            # Legend title and text font size
            legend.title = element_text(size = 17),    # Legend title size
            legend.text = element_text(size = 16)      # Legend labels size
        ) +
        stat_compare_means(
            comparisons = list(c(1,2), c(1,3), c(2,3)), 
            size=6
        ) +
        stat_summary(
            fun = "mean", 
            geom = "point",
            shape = 4,
            size = 2,
            stroke = 1.5
        )

    ggsave(plot_path, plot, width=15, height=10, dpi=720)
}

box_plot_libraries <- function(data, biotype, plot_path) {
    # data$`gene_type`[data$`gene_type` == "NON_CODING"] <- "non-coding"
    # data$`gene_type`[data$`gene_type` == "CODING"] <- "coding"
    # data$`gene_type`[data$`gene_type` == "All_GENES"] <- "all"
    data <- data[data$model == "RF", ]
    data <- data[data["gene_type"] == biotype, ]
    means <- data %>% group_by("dge_method") %>% summarize(m=mean(`auc`))
    print(means)
    plot <- ggplot(data, aes(x=dge_method, y=auc, fill=`dge_method`)) +
        geom_boxplot() +
        theme_minimal() +
        labs(x = "DGE Analysis Method", y = "auc", fill="RNA Sequencing Method") +
        scale_fill_manual(
            values = c("Union" = colors[9], "DESeq2" = colors[9], "EdgeR" = colors[9], "LimmaVoom" = colors[9], "Intersection" = colors[9]),
        ) +
        scale_x_discrete(
            labels = c("Union" = "Union", "DESeq2" = "DESeq2", "EdgeR" = "EdgeR", "LimmaVoom" = "LimmaVoom", "Intersection" = "Intersect")
        ) +
        scale_y_continuous(
            breaks = c(0.85, 0.9, 0.95, 1.0),
            labels = c("0.85", "0.90", "0.95", "1.00")
        ) +
        theme(
            legend.position = "none",
            # Axis titles font size
            axis.title.x = element_text(size = 19),    # x-axis label size
            axis.title.y = element_text(size = 19),    # y-axis label size

            # Axis tick labels font size
            axis.text.x = element_text(size = 17),     # x-axis ticks (labels) size
            axis.text.y = element_text(size = 17),     # y-axis ticks (labels) size

            # Legend title and text font size
            legend.title = element_text(size = 17),    # Legend title size
            legend.text = element_text(size = 16)      # Legend labels size
        ) +
        stat_compare_means(
            comparisons = list(
                c("Union", "DESeq2"), c("Union", "EdgeR"), c("Union", "LimmaVoom"), c("Union", "Intersection"),
                c("DESeq2", "EdgeR"), c("DESeq2", "LimmaVoom"), c("DESeq2", "Intersection"),
                c("EdgeR", "LimmaVoom"), c("EdgeR", "Intersection"),
                c("LimmaVoom", "Intersection")
            ), size=6
        ) +
        stat_summary(
            fun = "mean",
            geom = "point",
            shape = 4,
            size = 2,
            stroke = 1.5
        )

    ggsave(plot_path, plot, width=15, height=10, dpi=720)
}

make_boxplot_gene_expression <- function(ensembles, genes, plot_path) {
    counts <- read.csv("./data/sasc326_counts.csv")
    row.names(counts) <- counts[, 1]
    counts[, 1] <- NULL
    samples <- read.csv("./data/sasc326_samples.csv")
    dge <- DGEList(counts = counts, samples = samples)
    normfactors <-  calcNormFactors(dge, method= "TMM")
    normalized_counts <- cpm(normfactors, log=TRUE)

    hhc_t1 <- samples[samples$group %in% c("HHC", "First"), ][, c("sample", "group")]
    data_hhc_t1 <- normalized_counts[ensembles, ]
    data_hhc_t1 <- subset(data_hhc_t1, select=hhc_t1$sample)
    data_hhc_t1 <- data.frame(GeneName = genes, data_hhc_t1)
    data_hhc_t1 <- data_hhc_t1 %>% select(-c("s103830.003.011", "s103830.004.017"))

    melted_data <- reshape2::melt(data_hhc_t1)
    colnames(melted_data) <- c("GeneName", "SampleID", "CPM")

    hhc <- hhc_t1[hhc_t1$group == "HHC", ]$sample
    t1 <- hhc_t1[hhc_t1$group == "First", ]$sample
    melted_data <- melted_data %>%
        mutate(Condition=case_when(
            SampleID %in% hhc ~ "HHC",
            SampleID %in% t1 ~ "First",
            TRUE ~ "Second"
        ))


    plot <- ggplot(melted_data, aes(x = GeneName, y = CPM)) +
            geom_boxplot(aes(fill=factor(Condition, levels=c("HHC", "First"))), position = position_dodge(0.9)) +
            theme_minimal() +
            scale_fill_manual(
                values = c("HHC" = colors[7], "First" = colors[2]),
                labels = c("HHC" = "HHC", "First" = "First")
            ) +
            scale_y_continuous(
                breaks = c(-4, -2, 0, 2, 4, 8, 16),
                labels = c("-2", "-1", "0", "1", "2", "3", "4")
            ) +
            labs(y = "CPM", x="Gene names") +
            guides(fill=guide_legend("Condition")) +
            theme(
                legend.position = "right",
                # Axis titles font size
                axis.title.x = element_text(size = 19),    # x-axis label size
                axis.title.y = element_text(size = 19),    # y-axis label size

                # Axis tick labels font size
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 17),     # x-axis ticks (labels) size
                axis.text.y = element_text(size = 17),     # y-axis ticks (labels) size

                # Legend title and text font size
                legend.title = element_text(size = 17),    # Legend title size
                legend.text = element_text(size = 16)      # Legend labels size
            )


    ggsave(plot_path, plot, width=18, height=10, dpi=720)
}

mutate_plot_data_for_color <- function(data, log2_off, cut_off, coding, non_coding, genes) {
    if (is.null(genes)) {
        d <- data %>%
            mutate( Colorcode = case_when(
                Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "up_c" ,
                Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "down_c",
                Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "up_nc" ,
                Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "down_nc",
                TRUE ~ "not"
            ) )
    }
    else {
        d <- data %>%
            mutate( Colorcode = case_when(
                Row.names %in% genes & Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "up_c" ,
                Row.names %in% genes & Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "down_c",
                Row.names %in% genes & Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "up_nc" ,
                Row.names %in% genes & Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "down_nc",
                TRUE ~ "not"
            ) )
    }
    return(d)
}


VolcanoPlot_DGE <- function(data, log2_off, cut_off, coding, non_coding, genes, top, plot_path) {
    d <- mutate_plot_data_for_color(data, log2_off, cut_off, coding, non_coding, genes)
    text_repel_1 <- d[d$Colorcode != "not", ]
    if (!is.null(genes)) {
      text_repel_1 <- d[d$Row.names %in% genes, ]
    }
    if (!is.null(top)) {
      text_repel_1 <- head(text_repel_1, top)
    }
    p <- ggplot(data = d, aes( x = Log2FoldChange, y = -log10( P.Adjust ), color = Colorcode, label = gene_name) ) +
      geom_hline(yintercept = -log10(cut_off), linetype = "dashed", colour="darkgray", alpha=0.75) +
      geom_vline(xintercept = c(-log2_off, log2_off), linetype = "dashed", colour="darkgray", alpha=0.75) +
      theme_minimal() +
      geom_point() +
      geom_text_repel(data = text_repel_1, size=4.4) +
      labs(x = "Log2 Fold Change", y="-Log10(P-Adjust value)") +
      scale_color_manual(
        values = c( up_c = colors[8], down_c = colors[4], up_nc = colors[3], down_nc = colors[6], not = colors[10]),
        labels = c(up_c = "UP - CODING", down_c = "DOWN - CODING", up_nc = "UP - NC", down_nc = "DOWN - NC", not = "Stable state" )
      ) +
      scale_y_continuous(breaks = c(0, 1.00, -log10(cut_off), 2.00, 3.00, 4.00),
                         labels = c(0, 1.00, round(-log10(cut_off), digits = 2), 2, 3, 4)) +
      scale_x_continuous(breaks = c(-3, -2, -1, -log2_off, 0, log2_off, 1, 2, 3),
                         limits = c(-3, 3),
                         labels = c(-3, -2, -1, "-Log2(1.5)", 0 , "Log2(1.5)", 1, 2, 3)) +
      guides(color=guide_legend("Regulatory State")) +
      theme(
        # Axis titles font size
        axis.title.x = element_text(size = 19),    # x-axis label size
        axis.title.y = element_text(size = 19),    # y-axis label size

        # Axis tick labels font size
        axis.text.x = element_text(size = 17),     # x-axis ticks (labels) size
        axis.text.y = element_text(size = 17),     # y-axis ticks (labels) size

        # Legend title and text font size
        legend.title = element_text(size = 17),    # Legend title size
        legend.text = element_text(size = 16)      # Legend labels size
      )
    ggsave(plot_path, p, width=18, height=10, dpi=720)

}

# Run comparison plots
print("RUNNING COMPARISON PLOTS")
data <- readxl::read_excel("./data/Predictions/analysis_predictions.xlsx")
rfe_vs_chi2_plot(data, "./data/Plots/boxplot_rfe_vs_chi2.png")
type_comparison_plot(data, "./data/Plots/boxplot_rna_type.png")
box_plot_libraries(data, "ALL", "./data/Plots/boxplot_libraries_ALL_genes.png")
box_plot_libraries(data, "CODING", "./data/Plots/boxplot_libraries_coding_genes.png")
box_plot_libraries(data, "NC", "./data/Plots/boxplot_libraries_non_coding_genes.png")

# RUN gene expression top 30 genes union
dataset_union <- readxl::read_excel("./data/Predictions/best_dge_genes_union.xlsx")
subset_ensembles_union <- dataset_union$ensemble[1:30]
subset_genes_union <- dataset_union$gene_name[1:30]
make_boxplot_gene_expression(subset_ensembles_union, subset_genes_union, "./data/Plots/Union_Top_30_gene_expression_box_plot.png")

# RUN gene expression intersection
dataset_intersection <- readxl::read_excel("./data/Predictions/best_dge_genes_intersection.xlsx")
subset_ensembles_intersection <- dataset_intersection$ensemble
subset_genes_intersection <- dataset_intersection$gene_name
make_boxplot_gene_expression(subset_ensembles_intersection, subset_genes_intersection, "./data/Plots/Intersection_gene_expression_box_plot.png")

# # run Volcanoplot for results from every DGE analysis
print("RUNNING DGE VOLCANO PLOTS")
VolcanoPlot_DGE(
  readxl::read_excel( "./data/DESeq2/All_GENES_DESeq2_RNA_SEQ.xlsx"),
  LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, NULL, NULL,
  "./data/Plots/DESeq2_volcano_plot.png"
)
VolcanoPlot_DGE(
  readxl::read_excel( "./data/EdgeR/All_GENES_EdgeR_RNA_SEQ.xlsx"),
  LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, NULL, NULL,
  "./data/Plots/EdgeR_volcano_plot.png"
)
VolcanoPlot_DGE(
  readxl::read_excel( "./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"),
  LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, NULL, NULL,
  "./data/Plots/LimmaVoom_volcano_plot.png"
)

# run Volcanoplot for results from every DGE analysis with subset of genes found by Path A
dataset_rf <- readxl::read_excel("./data/Predictions/best_dge_genes_['RF'].xlsx")
subset_ensembles_rf <- dataset_rf$ensemble
VolcanoPlot_DGE(
    readxl::read_excel("./data/DESeq2/All_GENES_DESeq2_RNA_SEQ.xlsx"),
    LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles_rf, NULL,
    "./data/Plots/Subset_PathA_DESeq2_volcano_plot.png"
)
VolcanoPlot_DGE(
    readxl::read_excel("./data/EdgeR/All_GENES_EdgeR_RNA_SEQ.xlsx"),
    LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles_rf, NULL,
    "./data/Plots/Subset_PathA_EdgeR_volcano_plot.png"
)
VolcanoPlot_DGE(
    readxl::read_excel("./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"),
    LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles_rf, NULL,
    "./data/Plots/Subset_PathA_LimmaVoom_volcano_plot.png"
)

# run Volcanoplot for results from every DGE analysis with subset of genes found by Path B
dataset <- readxl::read_excel("./data/Predictions/best_dge_genes_['SVM', 'RF'].xlsx")
subset_ensembles_rf_svm <- dataset$ensemble
VolcanoPlot_DGE(
    readxl::read_excel("./data/DESeq2/All_GENES_DESeq2_RNA_SEQ.xlsx"),
    LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles_rf_svm, NULL,
    "./data/Plots/Subset_PathB_DESeq2_volcano_plot.png"
)
VolcanoPlot_DGE(
    readxl::read_excel("./data/EdgeR/All_GENES_EdgeR_RNA_SEQ.xlsx"),
    LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles_rf_svm, NULL,
    "./data/Plots/Subset_PathB_EdgeR_volcano_plot.png"
)
VolcanoPlot_DGE(
    readxl::read_excel("./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"),
    LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles_rf_svm, NULL,
    "./data/Plots/Subset_PathB_LimmaVoom_volcano_plot.png"
)

# Run final Volconaplot for comparison of all analytical approaches
data <- readxl::read_excel("./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx")
d <- data %>%
    mutate( Colorcode = case_when(
        Row.names %in% subset_ensembles_intersection ~ "Intersection",
        Row.names %in% subset_ensembles_rf ~ "A",
        Row.names %in% subset_ensembles_rf_svm ~ "B",
        TRUE ~ "not"
    ) )
text_repel_1 <- d[d$Row.names %in% unique(c(subset_ensembles_rf, subset_ensembles_rf_svm)), ]
p <- ggplot(data = d, aes( x = Log2FoldChange, y = -log10( P.Adjust ), color = factor(Colorcode, levels=c("A", "B", "Intersection", "not")), label = gene_name) ) +
    geom_hline(yintercept = -log10(CUT_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
    geom_vline(xintercept = c(-LOG2_OFF, LOG2_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
    theme_minimal() +
    geom_point() +
    geom_text_repel(data = text_repel_1, size=5.5) +
    labs(x = "Log2 Fold Change", y="-Log10(P-Adjust value)") +
    scale_color_manual(
        values = c(A = colors[2], B = colors[7], Intersection=colors[4], not = colors[10]),
        labels = c(A = "Path A", B = "Path B", Intersection="Intersection", not = "Other" )
    ) +
    scale_y_continuous(breaks = c(0, 1.00, -log10(CUT_OFF), 2.00, 3.00, 4.00),
                        labels = c(0, 1.00, round(-log10(CUT_OFF), digits = 2), 2, 3, 4)) +
    scale_x_continuous(breaks = c(-3, -2, -1, -LOG2_OFF, 0, LOG2_OFF, 1, 2, 3),
                        limits = c(-3, 3),
                        labels = c(-3, -2, -1, "-Log2(1.5)", 0 , "Log2(1.5)", 1, 2, 3)) +
    guides(color=guide_legend("Analytical approaches")) +
    theme(
    # Axis titles font size
    axis.title.x = element_text(size = 19),    # x-axis label size
    axis.title.y = element_text(size = 19),    # y-axis label size

    # Axis tick labels font size
    axis.text.x = element_text(size = 17),     # x-axis ticks (labels) size
    axis.text.y = element_text(size = 17),     # y-axis ticks (labels) size

    # Legend title and text font size
    legend.title = element_text(size = 17),    # Legend title size
    legend.text = element_text(size = 16)      # Legend labels size
    )
ggsave("./data/Plots/Comparison_analytical_approaches_Volcanoplot.png", p, width=18, height=10, dpi=720)
