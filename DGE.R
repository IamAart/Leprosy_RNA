# Loading R packages
library(DESeq2)
library(dplyr)

source("common.R")

# load in data with rownames
data <- prepare_data(FILEPATH, COUNTS_NAME, FEATURES_NAME, SAMPLES_NAME, BIOTYPES)
counts <- data[["counts"]]
features <- data[["features"]]
samples <- data[["samples"]]

# make sure participation names are at the same spot for counts and sample data
counts <- check_order_rownames(samples, counts)
if (is.logical(counts)) {
	stop("Error when adjusting colnames and rownames to be similar")
}

# perform DESEQ2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~ group)
deseq <- DESeq(dds)
res <- results(deseq, c("group", "First", "HHC"), pAdjustMethod = "BH")

# Save data in table with cut off added to cut the data at that p value
res <- make_save_table(res, features, "./data/DESeq2/NON_CODING_DESeq2_RNA_SEQ")