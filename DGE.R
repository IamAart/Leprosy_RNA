# Loading R packages
library(DESeq2)
library(Glimma)
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

glimmaMA(deseq)

make_save_table(res, "padj", CUT_OFF, "pvalue", features, "./data/DESeq2/all_nc_data-First-HHC")
