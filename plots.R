library("readxl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("edgeR")
library("ggpubr")
library("reshape2")
library("stringr")
library("dotenv")

.libPaths("./R_Packages/")
dotenv::load_dot_env(file = ".env")
NON_CODING_BIOTYPES <- as.vector(strsplit(Sys.getenv("NON_CODING_BIOTYPES"), ","))[[1]]
CODING_BIOTYPES <- as.vector(strsplit(Sys.getenv("CODING_BIOTYPES"), ","))[[1]]
CUT_OFF <- as.numeric(Sys.getenv("CUT_OFF"))
LOG2_OFF <- log2(as.numeric(Sys.getenv("LOG2_OFF")))

colors <- c("#88ccee", "#44aa99", "#117733", "#332288", "#ddcc77", "#999933", "#cc6677", "#882255", "#aa4499", "#dddddd")

rfe_vs_chi2_plot <- function(data, plot_path) {
    data$`Feature Selection Method` <- as.factor(data$`Feature Selection Method`)
    means <- data %>% group_by(`Feature Selection Method`) %>% summarize(m=mean(`Average AUC score`))
    print(means)
    plot <- ggplot(data, aes(x=`Feature Selection Method`, y=`Average AUC score`, fill=`Feature Selection Method`)) +
        geom_boxplot() +
        theme_minimal() +
        labs(x = "Feature Selection Method", y = "Average AUC score", color="Feature Selection Method") +
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
    data$`Analysis type`[data$`Analysis type` == "NON_CODING"] <- "NC"
    data$`Analysis type`[data$`Analysis type` == "CODING"] <- "CODING"
    data$`Analysis type`[data$`Analysis type` == "All_GENES"] <- "ALL"
    means <- data %>% group_by(`Analysis type`) %>% summarize(m=mean(`Average AUC score`))
    print(means)
    plot <- ggplot(data, aes(x=`Analysis type`, y=`Average AUC score`, fill=`Analysis type`)) +
        geom_boxplot() +
        labs(x = "Type of RNA", y = "Average AUC score", fill="Type of RNA") +
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
    data$`Analysis type`[data$`Analysis type` == "NON_CODING"] <- "non-coding"
    data$`Analysis type`[data$`Analysis type` == "CODING"] <- "coding"
    data$`Analysis type`[data$`Analysis type` == "All_GENES"] <- "all"
    data <- data[data["Analysis type"] == biotype, ]
    means <- data %>% group_by(`Library`) %>% summarize(m=mean(`Average AUC score`))
    print(means)
    plot <- ggplot(data, aes(x=`Library`, y=`Average AUC score`, fill=`Library`)) +
        geom_boxplot() +
        theme_minimal() +
        labs(x = "DGE Analysis Method", y = "Average AUC score", fill="RNA Sequencing Method") +
        scale_fill_manual(
            values = c("All" = colors[9], "DESeq2" = colors[9], "EdgeR" = colors[9], "LimmaVoom" = colors[9], "Only combined libraries" = colors[9]),
        ) +
        scale_x_discrete(
            labels = c("All" = "Union", "DESeq2" = "DESeq2", "EdgeR" = "EdgeR", "LimmaVoom" = "LimmaVoom", "Only combined libraries" = "Intersect")
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
                c("All", "DESeq2"), c("All", "EdgeR"), c("All", "LimmaVoom"), c("All", "Only combined libraries"),
                c("DESeq2", "EdgeR"), c("DESeq2", "LimmaVoom"), c("DESeq2", "Only combined libraries"),
                c("EdgeR", "LimmaVoom"), c("EdgeR", "Only combined libraries"),
                c("LimmaVoom", "Only combined libraries")
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
# print("RUNNING COMPARISON PLOTS")
data <- readxl::read_excel("./data/Predictions/analysis_rf_predictions.xlsx")
rfe_vs_chi2_plot(data, "./data/Plots/boxplot_rfe_vs_chi2.png")
type_comparison_plot(data, "./data/Plots/boxplot_rna_type.png")
box_plot_libraries(data, "all", "./data/Plots/boxplot_libraries_all_genes.png")
box_plot_libraries(data, "coding", "./data/Plots/boxplot_libraries_coding_genes.png")
box_plot_libraries(data, "non-coding", "./data/Plots/boxplot_libraries_non_coding_genes.png")

# subset_ensembles <- "ENSG00000215908,ENSG00000253683,ENSG00000179085,ENSG00000258920,ENSG00000234741,ENSG00000225864,ENSG00000222020,ENSG00000188290,ENSG00000272677,ENSG00000232229,ENSG00000167700,ENSG00000215414,ENSG00000145337,ENSG00000274012,ENSG00000278771,ENSG00000272906,ENSG00000273599,ENSG00000236552,ENSG00000228205,ENSG00000170889,ENSG00000204387,ENSG00000203875,ENSG00000269893,ENSG00000172803,ENSG00000184986,ENSG00000226287,ENSG00000233493,ENSG00000141933,ENSG00000226085,ENSG00000255559"
# subset_ensembles <- as.vector(strsplit(subset_ensembles, ","))[[1]]
# subset_genes <- "CROCCP2,CTB-79E8.3,DPM3,FOXN3-AS1,GAS5,HCG4P11,HDAC4-AS1,HES4,HNRNPD-DT,LINC00865,MFSD3,PSMA6P1,PYURF,RN7SL2,RN7SL3,RP11-533E19.7,RP11-59C5.3,RPL13AP5,RPS3P3,RPS9,SNHG32,SNHG5,SNHG8,SNX32,TMEM121,TMEM191A,TMEM238,TPGS1,UQCRFS1P1,ZNF252P-AS1"
# subset_genes <- as.vector(strsplit(subset_genes, ","))[[1]]
# make_boxplot_gene_expression(subset_ensembles, subset_genes, "./data/Plots/Top_30_gene_expression_box_plot.png")

# # run Volcanoplot for results from every DGE analysis
# print("RUNNING DGE VOLCANO PLOTS")
# VolcanoPlot_DGE(
#   readxl::read_excel( "./data/DESeq2/All_GENES_DESeq2_RNA_SEQ.xlsx"),
#   LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, NULL, NULL,
#   "./data/Plots/DESeq2_volcano_plot.png"
# )
# VolcanoPlot_DGE(
#   readxl::read_excel( "./data/EdgeR/All_GENES_EdgeR_RNA_SEQ.xlsx"),
#   LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, NULL, NULL,
#   "./data/Plots/EdgeR_volcano_plot.png"
# )
# VolcanoPlot_DGE(
#   readxl::read_excel( "./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"),
#   LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, NULL, NULL,
#   "./data/Plots/LimmaVoom_volcano_plot.png"
# )

# # run Volcanoplot for results from every DGE analysis with subset of genes found by Random Forest
# if (Sys.getenv("SUBSET_ENSEMBLES") != "") {
#   print("RUNNING SUBSET DGE VOLCANO PLOTS")
#   subset_ensembles <- as.vector(strsplit(Sys.getenv("SUBSET_ENSEMBLES"), ","))[[1]]
#   VolcanoPlot_DGE(
#     readxl::read_excel("./data/DESeq2/All_GENES_DESeq2_RNA_SEQ.xlsx"),
#     LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles, NULL,
#     "./data/Plots/Subset_DESeq2_volcano_plot.png"
#   )
#   VolcanoPlot_DGE(
#     readxl::read_excel("./data/EdgeR/All_GENES_EdgeR_RNA_SEQ.xlsx"),
#     LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles, NULL,
#     "./data/Plots/Subset_EdgeR_volcano_plot.png"
#   )
#   VolcanoPlot_DGE(
#     readxl::read_excel("./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"),
#     LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES, subset_ensembles, NULL,
#     "./data/Plots/Subset_LimmaVoom_volcano_plot.png"
#   )

#   if (Sys.getenv("SUBSET_GENE_NAMES") != "") {
#     print("RUNNING SUBSET GENE EXPRESSION PLOT")
#     subset_genes <- as.vector(strsplit(Sys.getenv("SUBSET_GENE_NAMES"), ","))[[1]]
#     make_boxplot_gene_expression(subset_ensembles, subset_genes, "./data/Plots/Subset_gene_expression_box_plot.png")
#   }
# }

# TODO: MAKE THIS AS A FUNCTION / INTEGRATED FUNCTION
# data <- readxl::read_excel("./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx")
# subset_ensembles <- as.vector(strsplit(Sys.getenv("SUBSET_ENSEMBLES"), ","))[[1]]
# subset_ensembles_2 <- as.vector(strsplit(Sys.getenv("SUBSET_ENSEMBLES_2"), ","))[[1]]
# intersection <- intersect(subset_ensembles, subset_ensembles_2)
# path_a <- setdiff(subset_ensembles, subset_ensembles_2)
# path_b <- setdiff(subset_ensembles_2, subset_ensembles)
# d <- data %>%
#     mutate( Colorcode = case_when(
#         Row.names %in% path_a ~ "DGE",
#         Row.names %in% path_b ~ "ML",
#         Row.names %in% intersection ~ "Intersection",
#         TRUE ~ "not"
#     ) )
# text_repel_1 <- d[d$Row.names %in% c(subset_ensembles, subset_ensembles_2), ]
# p <- ggplot(data = d, aes( x = Log2FoldChange, y = -log10( P.Adjust ), color = factor(Colorcode, levels=c("DGE", "ML", "Intersection", "not")), label = gene_name) ) +
#     geom_hline(yintercept = -log10(CUT_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
#     geom_vline(xintercept = c(-LOG2_OFF, LOG2_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
#     theme_minimal() +
#     geom_point() +
#     geom_text_repel(data = text_repel_1, size=5.5) +
#     labs(x = "Log2 Fold Change", y="-Log10(P-Adjust value)") +
#     scale_color_manual(
#         values = c(DGE = colors[2], ML = colors[7], Intersection=colors[4], not = colors[10]),
#         labels = c(DGE = "Path A", ML = "Path B", Intersection="Intersection", not = "Other" )
#     ) +
#     scale_y_continuous(breaks = c(0, 1.00, -log10(CUT_OFF), 2.00, 3.00, 4.00),
#                         labels = c(0, 1.00, round(-log10(CUT_OFF), digits = 2), 2, 3, 4)) +
#     scale_x_continuous(breaks = c(-3, -2, -1, -LOG2_OFF, 0, LOG2_OFF, 1, 2, 3),
#                         limits = c(-3, 3),
#                         labels = c(-3, -2, -1, "-Log2(1.5)", 0 , "Log2(1.5)", 1, 2, 3)) +
#     guides(color=guide_legend("Analytical approaches")) +
#     theme(
#     # Axis titles font size
#     axis.title.x = element_text(size = 19),    # x-axis label size
#     axis.title.y = element_text(size = 19),    # y-axis label size

#     # Axis tick labels font size
#     axis.text.x = element_text(size = 17),     # x-axis ticks (labels) size
#     axis.text.y = element_text(size = 17),     # y-axis ticks (labels) size

#     # Legend title and text font size
#     legend.title = element_text(size = 17),    # Legend title size
#     legend.text = element_text(size = 16)      # Legend labels size
#     )
# ggsave("./data/Plots/Comparison_analytical_approaches_Volcanoplot.png", p, width=18, height=10, dpi=720)
