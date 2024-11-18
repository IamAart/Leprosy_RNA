library("readxl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("reshape2")
library("stringr")



VolcanoPlot <- function( d, top_20 ) {
    text_repel <- d[d$Colorcode %in% c("up_c", "down_c", "up_nc", "down_nc"), ]
    if (top_20 == TRUE) {
      text_repel <- head(text_repel, 20)
    }
    p <- ggplot(data = d, aes( x = Log2FoldChange, y = -log10( P.Adjust ), color = Colorcode, label = gene_name) ) +
      geom_hline(yintercept = -log10(CUT_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
      geom_vline(xintercept = c(-LOG2_OFF, LOG2_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
      geom_point() +
      geom_text_repel(data = text_repel, size=3) +
      labs(title = "Volcano plot different expressed genes", x = "Log2 (Fold Change)", y="Adjusted P-Value") +
      scale_color_manual(
        values = c( up_c = "#ff7f0e", down_c = "blue", up_nc = "#e31a1c", down_nc = "darkgreen", not = "gray" ),
        labels = c(up_c = "Up-regulated coding genes", down_c = "Down-regulated coding genes", up_nc = "Up-regulated non-coding genes", down_nc = "Down-regulated non-coding genes", not = "Not differentially expressed genes" )
      ) +
      scale_y_continuous(breaks = c(0, -log10(CUT_OFF), 2.00, 3.00, 4.00),
                         labels = c(1, 0.05, 0.01, 0.001, 0.0001)) +
      scale_x_continuous(breaks = c(-3, -2, -1, -LOG2_OFF, 0, LOG2_OFF, 1, 2, 3),
                         limits = c(-3, 3),
                         labels = c(-3, -2, -1, "-Log2(1.5)", 0 , "Log2(1.5)", 1, 2, 3)) +
      guides(color=guide_legend("Gene Information")) +
      theme(legend.position = "right")
    return(p)
}

VolcanoPlot_Final_Results <- function( data, log2_off, cut_off, coding, non_coding ) {
    first <- c("ENSG00000204387", "ENSG00000233493", "ENSG00000203875", "ENSG00000269893", "ENSG00000274012")
    second <- c("ENSG00000177469", "ENSG00000179085", "ENSG00000215414", "ENSG00000225864", "ENSG00000226287", "ENSG00000236552", "ENSG00000255559", "ENSG00000258920", "ENSG00000272906", "ENSG00000278771")
    
    d <- data %>%
          mutate( Colorcode = case_when(
              Row.names %in% c(first, second) & Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "up_c" ,
              Row.names %in% c(first, second) & Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "down_c",
              Row.names %in% c(first, second) & Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "up_nc" ,
              Row.names %in% c(first, second) & Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "down_nc",
              TRUE ~ "not"
          ) )
    text_repel_1 <- d[d$Row.names %in% c(first), ]
    text_repel_2 <- d[d$Row.names %in% c(second), ]
    p <- ggplot(data = d, aes( x = Log2FoldChange, y = -log10( P.Adjust ), color = Colorcode, label = gene_name) ) +
      geom_hline(yintercept = -log10(CUT_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
      geom_vline(xintercept = c(-LOG2_OFF, LOG2_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
      geom_point() +
      geom_text_repel(data = text_repel_1, size=3, colour="purple") +
      geom_text_repel(data = text_repel_2, size=3, colour="darkgreen") +
      labs(x = "Log2 (Fold Change)", y="Adjusted P-Value") +
      scale_color_manual(
        values = c( up_c = "#ff7f0e", down_c = "blue", up_nc = "#e31a1c", down_nc = "darkgreen", not = "gray" ),
        labels = c(up_c = "Up-regulated coding genes", down_c = "Down-regulated coding genes", up_nc = "Up-regulated non-coding genes", down_nc = "Down-regulated non-coding genes", not = "Not differentially expressed genes" )
      ) +
      scale_y_continuous(breaks = c(0, -log10(CUT_OFF), 2.00, 3.00, 4.00),
                         labels = c(1, 0.05, 0.01, 0.001, 0.0001)) +
      scale_x_continuous(breaks = c(-3, -2, -1, -LOG2_OFF, 0, LOG2_OFF, 1, 2, 3),
                         limits = c(-3, 3),
                         labels = c(-3, -2, -1, "-Log2(1.5)", 0 , "Log2(1.5)", 1, 2, 3)) +
      guides(color=guide_legend("Gene Information")) +
      theme(legend.position = "right")
    return(p)
}

get_subset_genes_all <- function(path) {
  data <- read_excel(path)["ensembles"]
  total_list = c()
  for(i in 1:nrow(data)) {
    vector = str_split(data[i, ], ", ")
    total_list <- append(total_list, vector[[1]])
  }
  return(unlist(unique(total_list), use.names = FALSE))
}

get_subset_genes_best_of_three <- function(path) {
  total_list <- c()
  data <- read_excel(path)

  all_data <- data[data$gene_type == "all", ]
  vector_all <- str_split(all_data[which.max(all_data$auc), "ensembles"], ", ")
  total_list <- append(total_list, vector_all[[1]])
  
  c_data <- data[data$gene_type == "coding", ]
  vector_c <- str_split(c_data[which.max(c_data$auc), "ensembles"], ", ")
  total_list <- append(total_list, vector_c[[1]])

  nc_data <- data[data$gene_type == "nc", ]
  vector_nc <- str_split(nc_data[which.max(nc_data$auc), "ensembles"], ", ")
  total_list <- append(total_list, vector_nc[[1]])
  
  return(unlist(unique(total_list), use.names = FALSE))
}

get_subset_genes_best <- function(path) {
  total_list <- c()
  data <- read_excel(path)
  vector <- str_split(data[which.max(data$auc), "ensembles"], ", ")
  total_list <- append(total_list, vector[[1]])
  
  return(unlist(unique(total_list), use.names = FALSE))
}



run_volcanoplot <- function(data_path_list, cut_off, log2_off, non_coding, coding, plot_names, subset) {
  for (i in seq_along(data_path_list)) {
    # load dge analysis data
    data <- read_excel(data_path_list[[i]])
    if (is.null(subset)) { # if no subset, create data without filtering on gene_names
      result_data <- data %>%
          mutate( Colorcode = case_when(
              Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "up_c" ,
              Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "down_c",
              Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "up_nc" ,
              Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "down_nc",
              TRUE ~ "not"
          ) )
    } else {
      if (subset == "All") { # if subset, filter on row.names with the subset
        subset_genes <- get_subset_genes_all("./data/Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx")
      } else if (subset == "Best_of_Three") {
        subset_genes <- get_subset_genes_best_of_three("./data/Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx")
        print(subset_genes)
      } else if (subset == "Best") {
        subset_genes <- get_subset_genes_best("./data/Current_Comparison_SVM_RF/rf_analysis_9_22.xlsx")
      }
      result_data <- data %>%
            mutate( Colorcode = case_when(
                Row.names %in% subset_genes & Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "up_c" ,
                Row.names %in% subset_genes & Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "down_c",
                Row.names %in% subset_genes & Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "up_nc" ,
                Row.names %in% subset_genes & Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "down_nc",
                TRUE ~ "not"
            ) )
    }
      
    ggplot_build(VolcanoPlot(result_data, FALSE))
    ggsave(paste0("./data/VolcanoPlots/", plot_names[i], ".png"), VolcanoPlot(result_data, FALSE), width=15, height=10, dpi=300)
  }
}

CUT_OFF <- 0.05
LOG2_OFF <- log2(1.5)
NON_CODING_BIOTYPES <- c("unitary_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_unprocessed_pseudogene", "antisense", "transcribed_unitary_pseudogene", "polymorphic_pseudogene", "lincRNA", "sense_intronic", "transcribed_processed_pseudogene", "sense_overlapping", "IG_V_pseudogene", "pseudogene", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "snRNA", "miRNA", "misc_RNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "TR_V_pseudogene", "TR_J_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "scRNA", "scaRNA", "vaultRNA", "sRNA", "macro_lncRNA", "non_coding", "IG_pseudogene", "processed_transcript", "ribozyme")
CODING_BIOTYPES <- c("IG_D_gene", "protein_coding", "TR_V_gene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene", "TR_D_gene", "TEC")

paths <- list(
  "./data/DESeq2/All_GENES_DESeq2_RNA_SEQ.xlsx",
  "./data/EdgeR/All_GENES_EdgeR_RNA_SEQ.xlsx",
  "./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"
)

plot_names <- list(
  "DESeq2_DGE_analysis_volcano",
  "EdgeR_DGE_analysis_volcano",
  "LimmaVoom_DGE_analysis_volcano"
)

subset_genes <- NULL # "Best"  # "All" NULL "Best_of_Three"
# run_volcanoplot(paths, CUT_OFF, LOG2_OFF, NON_CODING_BIOTYPES, CODING_BIOTYPES, plot_names, subset_genes)


ggplot_build(VolcanoPlot_Final_Results(read_excel("./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"), LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES))
ggsave(paste0("./data/VolcanoPlots/Final_RF.png"), VolcanoPlot_Final_Results(read_excel("./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"), LOG2_OFF, CUT_OFF, CODING_BIOTYPES, NON_CODING_BIOTYPES), width=15, height=10, dpi=300)