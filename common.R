library("writexl")
library("dplyr")
library("ggplot2")
library("ggrepel")

BIOTYPES <- c("unitary_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_unprocessed_pseudogene", "antisense", "transcribed_unitary_pseudogene", "polymorphic_pseudogene", "lincRNA", "sense_intronic", "transcribed_processed_pseudogene", "sense_overlapping", "IG_V_pseudogene", "pseudogene", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "snRNA", "miRNA", "misc_RNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "TR_V_pseudogene", "TR_J_pseudogene", "IG_D_gene", "IG_C_pseudogene", "IG_J_pseudogene", "scRNA", "scaRNA", "vaultRNA", "sRNA", "macro_lncRNA", "non_coding", "IG_pseudogene")
CUT_OFF <- 0.05
LOG2_OFF <- 0.6
FILEPATH <- "./data/sasc326.rda"
COUNTS_NAME <- "countsMatrix"
FEATURES_NAME <- "featuresData"
SAMPLES_NAME <- "samplesData"
VENN_NAMES <- "gene_name" # "gene_name" or "Row.names"

load_one_rdata_object <- function(file_path) {
    # Gets one object from rdata (in our case "gte") and returns it
    # Benefit: variable name can be defined and used easily
    res <- local(
        {
            load(file_path)
            return(get(ls()))
        }
    )
    return(res)
}

column_to_rownames <- function(df, column_number) {
    # Function to swap the first column to the rownames
    rownames(df) <- df[,column_number]
    df[,column_number] <- NULL
    df
}

check_order_rownames <- function (data1, data2) {
    # Check whether the rownames of the samples are in the same order as the columns of the counts.
    # If not adjust so that counts columns are adjusted to the standard of the rownames in the samples
    # If now all rownames of the counts exist in the colnames of the counts return false which triggers an error
    if (!all(rownames(data1) %in% colnames(data2))) {
        return(FALSE)
    }
    if (!all(rownames(data1) == colnames(data2))) {
        data2 <- data2[, rownames(data1)]
    }
    return(data2)
}

prepare_data <- function(path, counts_name, features_name, samples_name, biotypes) {
    # Prepare data so that it easy to use during normalization
    # It adjusts the rownames to the column containing the ensemble ids or the samples to the rownames. It also selects all non-coding genes
    # Returns a list containing the 3 data structures counts, samples and features
    data <- load_one_rdata_object(path)

    counts <- column_to_rownames(data[[counts_name]], 1)
    features <- column_to_rownames(data[[features_name]], 1)
    counts <- counts[merge(counts, features, by = 'row.names')$biotype %in% biotypes,]

    samples <- data[[samples_name]]
    samples <- tibble::column_to_rownames(samples, "sample")
    return(list(counts = counts, features = features, samples = samples))
}

normalization_edgeR <- function(counts, samples, design) {
    # Normalize data from of the counts and sample and put them into a DGE List before going into the DGE Analysis
    normalized_data <- DGEList(counts = counts, samples = samples)

    keep <- filterByExpr(y = normalized_data, design)
    normalized_data <- normalized_data[keep, , keep.lib.sizes = FALSE]

    normalized_data <- calcNormFactors(normalized_data, method = "TMM")
    return(normalized_data)
}

make_save_table <- function(result, features, filepath_name) {
    # Saves the result of all three libraries to an excel file and a csv file
    data <- merge(as.data.frame(result), features, by="row.names")
    names(data)[names(data) %in% c("padj", "FDR", "adj.P.Val")] <- "P.Adjust" # DESeq2, EdgeR, LimmaVoom
    names(data)[names(data) %in% c("log2FoldChange", "logFC")] <- "Log2FoldChange" # DESeq2, EdgeR & LimmaVoom
    data <- data[order(data$P.Adjust), ]
    data <- data %>%
        mutate( col = case_when(
            P.Adjust <= CUT_OFF & Log2FoldChange >= LOG2_OFF ~ "up",
            P.Adjust <= CUT_OFF & Log2FoldChange < -LOG2_OFF ~ "down",
            TRUE ~ "not"
        ) )
    # data <- data[which(data[["P.Adjust"]] <= CUT_OFF),]

    write.csv(data, file = paste0(filepath_name, ".csv"))
    writexl::write_xlsx(data, paste0(filepath_name, ".xlsx"))
    return(data)
}

VolcanoPlot <- function( d ) {
    # text_label <-
    p <- ggplot( d ) +
        aes( x = Log2FoldChange, y = -log10( P.Adjust ), color = col, label = gene_name) +
        geom_hline(yintercept = -log10(CUT_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
        geom_vline(xintercept = c(-LOG2_OFF, LOG2_OFF), linetype = "dashed", colour="darkgray", alpha=0.75) +
        geom_point() +
        geom_text_repel(data = head(d[d$col %in% c("up", "down"), ], 20), size=3) +
        # geom_text_repel(data = tail(d, n=length(d)-21) ) +
        # geom_text_repel( data = d %>% filter( col != "not" ) ) +
        # geom_text_repel( data = d %>% filter( col == "not" ), size = 2 ) +
        xlab("Log2 (Fold Change)") +
        ylab("Adjusted P-Value") +
        scale_y_continuous(breaks = c(0, -log10(CUT_OFF), 2.00, 3.00, 4.00),
                           labels = c(1, "<0.05", 0.01, 0.001, ">0.0001")) +
        scale_x_continuous(breaks = c(seq(-7, 7, 2)), limits = c(-7, 7)) +
        theme_classic() +
        scale_color_manual( values = c( up = "red", down = "blue", not = "gray" ), guide = "none" )
    return(p)
}