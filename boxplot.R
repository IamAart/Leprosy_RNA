library("readxl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("viridis")
library("edgeR")

rfe_vs_chi2_plot <- function(data) {
    data$feature_selection <- as.factor(data$feature_selection)
    plot <- ggplot(data) +
        geom_boxplot(aes(x=feature_selection, y=auc, color=feature_selection)) +
        labs(title = "Difference in auc score between feature selection methods; RFE and Chi2", x = "Feature Selection Method", y = "Auc Score", color="Feature Selection Method") + 
        scale_color_viridis(
            discrete = TRUE,
            labels = c(rfe = "Recursive Feature Elimination", chi2 = "Chi-squared Feature Elimination")
        ) +
        theme(legend.position = "right")
    
    ggsave("./data/RF_analysis_plots/boxplot_rfe_vs_chi2.png", plot, width=7, height=7, dpi=490)
}

type_comparison_plot <- function(data) {
    data <- data[data$feature_selection == "rfe", ]
    data$gene_type[data$gene_type == "nc"] <- "non_coding" 
    plot <- ggplot(data) +
        geom_boxplot(aes(x=gene_type, y=auc, color=gene_type)) +
        labs(title = "Difference in auc score between different biotypes", x = "Biotype", y = "Auc Score", color="Biotype Information") + 
        scale_color_viridis(
            discrete = TRUE,
            labels = c(all = "All Genes", coding = "Coding genes", non_coding = "Non-coding genes")
        ) +
        theme(legend.position = "right")
    
    ggsave("./data//RF_analysis_plots/boxplot_non_coding_vs_coding.png", plot, width=7, height=7, dpi=490) 
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
            labels = c("All" = "genes from the union of all libraries", "DESeq2" = "DESeq2 genes", "EdgeR" = "EdgeR genes", "LimmaVoom" = "LimmaVoom genes", "Only combined libraries" = "genes from the intersect of all libraries")
        ) +
        theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
    
    ggsave("./data/RF_analysis_plots/boxplot_libraries_non_coding_genes.png", plot, width=12, height=6, dpi=720) 
}

make_boxplot_gene_expression <- function() {
    counts <- read.csv("./data/sasc326_counts.csv")
    row.names(counts) <- counts[, 1]
    counts[, 1] <- NULL
    samples <- read.csv("./data/sasc326_samples.csv")
    dge <- DGEList(counts = counts, samples = samples)
    normfactors <-  calcNormFactors(dge, method= "TMM")
    normalized_counts <- cpm(normfactors, log=TRUE)

    genes <- c("ENSG00000204387", "ENSG00000233493", "ENSG00000179085", "ENSG00000236552", "ENSG00000225864", "ENSG00000272906", "ENSG00000269893", "ENSG00000226287", "ENSG00000203875", "ENSG00000274012", "ENSG00000278771", "ENSG00000255559", "ENSG00000258920", "ENSG00000215414", "ENSG00000177469")
    gene_names <- c("SNHG32 or C6orf48", "TMEM238", "DPM3", "RPL13AP5", "HCG4P11", "RP11-533E19.7", "SNHG8", "TMEM191A", "SNHG5", "RN7SL2", "Metazoa SRP", "ZNF252P-AS1", "FOXN3-AS1", "PSMA6P1", "PTRF")

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
            geom_boxplot(aes(fill=factor(Condition, levels=c("HHC", "First"), labels = c("HHC" = "Household Contacts", "First" = "Early Progressors"))), position = position_dodge(0.9)) +
            
            labs(y = "Counts Per Million (Log2(x))", x="Gene Names") +
            guides(fill=guide_legend("Condition")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            

    ggsave("./data/RF_analysis_plots/boxplot_final_result.png", plot, width=12, height=6, dpi=720) 
}

# data <- read_excel("./data/Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx")
# rfe_vs_chi2_plot(data)
# type_comparison_plot(data)
# # box_plot_libraries(data, "all")
# # box_plot_libraries(data, "coding")
# box_plot_libraries(data, "nc")
make_boxplot_gene_expression()
