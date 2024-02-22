library("writexl")

BIOTYPES <- c("non_coding", "lincRNA", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "macro_lncRNA")
CUT_OFF <- 0.05
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

make_save_table <- function(result, cuf_off_col_name, cuf_off_val, order_col_name, features, filepath_name) {
    # Saves the result of all three libraries to an excel file and a csv file
    result <- result[which(result[[cuf_off_col_name]] <= cuf_off_val),]
    resOrdered <- result[order(result[[order_col_name]]),]

    data <- merge(as.data.frame(resOrdered), features, by="row.names")

    write.csv(data, file = paste0(filepath_name, ".csv"))
    writexl::write_xlsx(data, paste0(filepath_name, ".xlsx"))
}