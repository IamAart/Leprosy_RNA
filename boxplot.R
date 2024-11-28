library("readxl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("viridis")
library("edgeR")
library("ggpubr")

colors = c("#88ccee", "#44aa99", "#117733", "#332288", "#ddcc77", "#999933", "#cc6677", "#882255", "#aa4499", "#dddddd")

rfe_vs_chi2_plot <- function(data) {
    data$feature_selection <- as.factor(data$feature_selection)
    plot <- ggplot(data, aes(x=feature_selection, y=auc, fill=feature_selection)) +
        geom_boxplot() +
        theme_minimal() +
        labs(x = "Feature Selector", y = "Average Auc score", color="Feature selector") + 
        scale_fill_manual(
            values = c(rfe = colors[7], chi2 = colors[2]),
            labels = c(rfe = "RFE", chi2 = expression(chi^2))
        ) +
        scale_x_discrete(
            labels = c("rfe" = "RFE", "chi2" = expression(chi^2))
        ) +
        scale_y_continuous(
            breaks = c(0.85, 0.9, 0.95, 1.0),
            limits = c(0.85, 1.00),
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
        stat_compare_means(comparisons = list(c("chi2", "rfe")))
    
    ggsave("./data/RF_analysis_plots/boxplot_rfe_vs_chi2_mann_whitney_9.png", plot, width=15, height=10, dpi=720)
}

type_comparison_plot <- function(data) {
    # data <- data[data$feature_selection == "rfe", ]
    data$gene_type[data$gene_type == "nc"] <- "non_coding" 
    plot <- ggplot(data, aes(x=gene_type, y=auc, fill=gene_type)) +
        geom_boxplot() +
        
        labs(x = "Type of RNA", y = "Average Auc score", fill="Type of RNA") + 
        scale_fill_manual(
            values = c(all = colors[1], non_coding = colors[1], coding = colors[1]),
        ) +
        scale_x_discrete(
            limits = c("all", "non_coding", "coding"),
            labels = c("all" = "ALL", "non_coding" = "NC", "coding" = "CODING")
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
        stat_compare_means(comparisons = list(c("all", "coding"), c("coding", "non_coding"), c("all", "non_coding")))
    
    ggsave("./data//RF_analysis_plots/boxplot_non_coding_vs_coding_mann_whitney.png", plot, width=15, height=10, dpi=720) 
}

box_plot_libraries <- function(data, biotype) {
    # data <- data[data$feature_selection == "rfe", ]
    data <- data[data["gene_type"] == biotype, ]
    if (biotype == "nc") {
        data$gene_type[data$gene_type == "nc"] <- "non_coding" 
    }
    plot <- ggplot(data, aes(x=Library, y=auc, fill=Library)) +
        geom_boxplot() +
        theme_minimal() +
        labs(x = "DGE Analysis Method", y = "Average Auc score", fill="RNA Sequencing Method") + 
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
            )
            # method = "t.test"
        )
    
    ggsave("./data/RF_analysis_plots/boxplot_libraries_non_coding_genes_wilcox_mann_whitney.png", plot, width=12, height=6, dpi=720) 
}

make_boxplot_gene_expression <- function() {
    counts <- read.csv("./data/sasc326_counts.csv")
    row.names(counts) <- counts[, 1]
    counts[, 1] <- NULL
    samples <- read.csv("./data/sasc326_samples.csv")
    dge <- DGEList(counts = counts, samples = samples)
    normfactors <-  calcNormFactors(dge, method= "TMM")
    normalized_counts <- cpm(normfactors, log=TRUE)

    genes <- c('ENSG00000204387', 'ENSG00000233493', 'ENSG00000179085', 'ENSG00000170889', 'ENSG00000236552', 'ENSG00000225864', 'ENSG00000141933', 'ENSG00000272906', 'ENSG00000215908', 'ENSG00000269893', 'ENSG00000226287', 'ENSG00000203875', 'ENSG00000274012', 'ENSG00000278771', 'ENSG00000145337', 'ENSG00000255559', 'ENSG00000258920', 'ENSG00000215414', 'ENSG00000234741', 'ENSG00000253683')
    gene_names <- c('SNHG32 or C6orf48', 'TMEM238', 'DPM3', 'RPS9', 'RPL13AP5', 'HCG4P11', 'TPGS1', 'RP11-533E19.7', 'CROCCP2', 'SNHG8', 'TMEM191A', 'SNHG5', 'RN7SL2', 'Metazoa_SRP', 'PYURF', 'ZNF252P-AS1', 'FOXN3-AS1', 'PSMA6P1', 'GAS5', 'CTB-79E8.3')

    hhc_t1 <- samples[samples$group %in% c("HHC", "First"), ][, c("sample", "group")]
    hhc <- hhc_t1[hhc_t1$group == "HHC", ]$sample
    t1 <- hhc_t1[hhc_t1$group == "First", ]$sample
    print(head(hhc))

    data_hhc_t1 <- normalized_counts[genes, ]
    data_hhc_t1 <- subset(data_hhc_t1, select=hhc_t1$sample)
    data_hhc_t1 <- data.frame(GeneName = gene_names, data_hhc_t1)
    data_hhc_t1 <- data_hhc_t1 %>% select(-c("s103830.003.011", "s103830.004.017"))

    melted_data <- reshape2::melt(data_hhc_t1)
    colnames(melted_data) <- c("GeneName", "SampleID", "CPM")

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
                labels = c("HHC" = "Household Contacts", "First" = "Progressors without symptoms")
            ) +
            scale_y_continuous(
                breaks = c(-4, -2, 0, 2, 4, 8, 16),
                labels = c("-2", "-1", "0", "1", "2", "3", "4")
            ) + 
            labs(y = "CPM", x="Gene Names") +
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
            

    ggsave("./data/RF_analysis_plots/boxplot_final_result.png", plot, width=18, height=10, dpi=720) 
}

data <- read_excel("./data/Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx")


# rfe_vs_chi2_plot(data)
# type_comparison_plot(data)
# box_plot_libraries(data, "all")
# box_plot_libraries(data, "coding")
# box_plot_libraries(data, "nc")
# make_boxplot_gene_expression()
