library("writexl")
library("dplyr")
library("ggplot2")
library("ggrepel")

BIOTYPES <- c("unitary_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_unprocessed_pseudogene", "antisense", "transcribed_unitary_pseudogene", "polymorphic_pseudogene", "lincRNA", "sense_intronic", "transcribed_processed_pseudogene", "sense_overlapping", "IG_V_pseudogene", "pseudogene", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "snRNA", "miRNA", "misc_RNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "TR_V_pseudogene", "TR_J_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "scRNA", "scaRNA", "vaultRNA", "sRNA", "macro_lncRNA", "non_coding", "IG_pseudogene", "processed_transcript", "ribozyme")
CODING_BIOTYPES <- c("IG_D_gene", "protein_coding", "TR_V_gene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene", "TR_D_gene", "TEC")

CUT_OFF <- 0.05
LOG2_OFF <- log2(1.5)

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
    counts <- merge(counts, features, by = 'row.names')

    samples <- data[[samples_name]]
    samples <- tibble::column_to_rownames(samples, "sample")
    return(list(counts = counts, features = features, samples = samples))
}

normalization_edgeR <- function(counts, samples, design) {
    # Normalize data from of the counts and sample and put them into a DGE List before going into the DGE Analysis
    normalized_data <- DGEList(counts = counts, samples = samples)
    # TODO this reduces the amount of rows to 17.000
    keep <- filterByExpr(y = normalized_data, design)
    normalized_data <- normalized_data[keep, keep.lib.sizes = FALSE]

    normalized_data <- calcNormFactors(normalized_data, method = "TMM")
    # plotMDS(normalized_data)
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

    write.csv(data, file = paste0(filepath_name, ".csv"))
    writexl::write_xlsx(data, paste0(filepath_name, ".xlsx"))
    return(data)
}