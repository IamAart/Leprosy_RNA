library("readxl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("reshape2")

VolcanoPlot <- function( d ) {
    p <- ggplot(data = d, aes( x = Log2FoldChange, y = -log10( P.Adjust ), color = Colorcode, label = gene_name) ) +
      geom_hline(yintercept = -log10(CUT_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
      geom_vline(xintercept = c(-LOG2_OFF, LOG2_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
      geom_point() +
      geom_text_repel(data = head(d[d$Colorcode %in% c("up_c", "down_c", "up_nc", "down_nc"), ], 20), size=3) +
      labs(title = "Volcano plot different expressed genes", x = "Log2 (Fold Change)", y="Adjusted P-Value") +
      scale_color_manual(
        values = c( up_c = "red", down_c = "blue", up_nc = "darkorange", down_nc = "darkgreen", not = "gray" ),
        labels = c(up_c = "Up-regulated coding genes", down_c = "Down-regulated coding genes", up_nc = "Up-regulated non-coding genes", down_nc = "Down-regulated non-coding genes", not = "Not differentially expressed genes" )
      ) +
      scale_y_continuous(breaks = c(0, -log10(CUT_OFF), 2.00, 3.00, 4.00),
                         labels = c(1, "<0.05", 0.01, 0.001, ">0.0001")) +
      scale_x_continuous(breaks = c(seq(-3, 3, 1)), limits = c(-3, 3)) +
      guides(color=guide_legend("Gene Information")) +
      theme(legend.position = "right")
    return(p)
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
    } else { # if subset, filter on row.names with the subset
      result_data <- data %>%
          mutate( Colorcode = case_when(
              Row.names %in% subset & Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "up_c" ,
              Row.names %in% subset & Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% coding ~ "down_c",
              Row.names %in% subset & Log2FoldChange >= log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "up_nc" ,
              Row.names %in% subset & Log2FoldChange <= -log2_off & P.Adjust <= cut_off & biotype %in% non_coding ~ "down_nc",
              TRUE ~ "not"
          ) )
    }
    ggplot_build(VolcanoPlot(result_data))
    ggsave(paste0("./data/VolcanoPlots/", plot_names[i], ".png"), VolcanoPlot(result_data))
  }
}

CUT_OFF <- 0.001
LOG2_OFF <- log2(1.5)
NON_CODING_BIOTYPES <- c("unitary_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_unprocessed_pseudogene", "antisense", "transcribed_unitary_pseudogene", "polymorphic_pseudogene", "lincRNA", "sense_intronic", "transcribed_processed_pseudogene", "sense_overlapping", "IG_V_pseudogene", "pseudogene", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "snRNA", "miRNA", "misc_RNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "TR_V_pseudogene", "TR_J_pseudogene", "IG_D_gene", "IG_C_pseudogene", "IG_J_pseudogene", "scRNA", "scaRNA", "vaultRNA", "sRNA", "macro_lncRNA", "non_coding", "IG_pseudogene")
CODING_BIOTYPES <- c("protein_coding", "processed_transcript", "TR_V_gene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene", "ribozyme", "TR_D_gene", "TEC")

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

run_volcanoplot(paths, CUT_OFF, LOG2_OFF, NON_CODING_BIOTYPES, CODING_BIOTYPES, plot_names, NULL)
