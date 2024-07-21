library("readxl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("reshape2")

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
    ggplot_build(VolcanoPlot(result_data, FALSE))
    ggsave(paste0("./data/VolcanoPlots/", plot_names[i], ".png"), VolcanoPlot(result_data, FALSE))
  }
}

CUT_OFF <- 0.05
LOG2_OFF <- log2(1.5)
NON_CODING_BIOTYPES <- c("unitary_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_unprocessed_pseudogene", "antisense", "transcribed_unitary_pseudogene", "polymorphic_pseudogene", "lincRNA", "sense_intronic", "transcribed_processed_pseudogene", "sense_overlapping", "IG_V_pseudogene", "pseudogene", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "snRNA", "miRNA", "misc_RNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "TR_V_pseudogene", "TR_J_pseudogene", "IG_D_gene", "IG_C_pseudogene", "IG_J_pseudogene", "scRNA", "scaRNA", "vaultRNA", "sRNA", "macro_lncRNA", "non_coding", "IG_pseudogene")
CODING_BIOTYPES <- c("protein_coding", "processed_transcript", "TR_V_gene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene", "ribozyme", "TR_D_gene", "TEC")

paths <- list(
  "./data/DESeq2/All_GENES_DESeq2_RNA_SEQ.xlsx",
  "./data/EdgeR/All_GENES_EdgeR_RNA_SEQ.xlsx",
  "./data/LimmaVoom/All_GENES_LimmaVoom_RNA_SEQ.xlsx"
)
plot_names <- list(
  "test_DESeq2_DGE_analysis_volcano",
  "test_EdgeR_DGE_analysis_volcano",
  "test_LimmaVoom_DGE_analysis_volcano"
)
subset_genes <- list('ENSG00000272906', 'ENSG00000269893', 'ENSG00000175602', 'ENSG00000161179', 'ENSG00000272677', 'ENSG00000278771', 'ENSG00000266538', 'ENSG00000150991', 'ENSG00000265393', 'ENSG00000215030', 'ENSG00000198763', 'ENSG00000273599', 'ENSG00000111678', 'ENSG00000188290', 'ENSG00000213442', 'ENSG00000152082', 'ENSG00000277383', 'ENSG00000226085', 'ENSG00000225864', 'ENSG00000276791', 'ENSG00000006015', 'ENSG00000237842', 'ENSG00000240877', 'ENSG00000204387', 'ENSG00000179085', 'ENSG00000228205', 'ENSG00000233493', 'ENSG00000225178', 'ENSG00000198804', 'ENSG00000171858', 'ENSG00000232229', 'ENSG00000222020', 'ENSG00000141933', 'ENSG00000167700', 'ENSG00000236552', 'ENSG00000130748', 'ENSG00000272256', 'ENSG00000103254', 'ENSG00000203875', 'ENSG00000155428', 'ENSG00000231351')
run_volcanoplot(paths, CUT_OFF, LOG2_OFF, NON_CODING_BIOTYPES, CODING_BIOTYPES, plot_names, subset_genes)
