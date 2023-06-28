source("./dge_analysis.R")

COUNT_DATA_FILENAME <- "./data/sasc326.rda"
COUNTS_VARIABLE_NAME <- "countsMatrix"
SAMPLES_VARIABLE_NAME <- "samplesData"
FEATURES_VARIABLE_NAME <- "featuresData"
DGE_TYPE <- "LimmaVoom" # Can be "DESeq2", "LimmaVoom" or "EdgeR"
ADJUST_COMMON <- TRUE
INCLUDE_FEATURES <- TRUE
LOG2_LIMIT <- 2.0
PVAL_LIMIT <- 0.05

COMPARISON <- "First-HHC"
DESIGN <- ~ group

count_data <- load_one_rdata_object(COUNT_DATA_FILENAME)
dge_result <- perform_DGE_analysis(
    data = count_data,
    dge_type = DGE_TYPE,
    adjust_common = ADJUST_COMMON,
    include_features = INCLUDE_FEATURES,
    design = DESIGN,
    c_name = COUNTS_VARIABLE_NAME,
    f_name = FEATURES_VARIABLE_NAME,
    s_name = SAMPLES_VARIABLE_NAME
)

if (DGE_TYPE == "DESeq2") {
	data <- as.data.frame(dge_result$result)
	genes <- dge_result$genes
} else if (DGE_TYPE == "LimmaVoom") {
	data <- topTable(dge_result$result, n=Inf)
	genes <- dge_result$genes
} else if (DGE_TYPE == "EdgeR") {
	data <-dge_result$result$table
	genes <- dge_result$genes
}

save_and_show(
	data = as.data.frame(data),
	genes = genes,
	dge_type = DGE_TYPE,
	comparison = COMPARISON,
	p_val = PVAL_LIMIT,
	log_val = LOG2_LIMIT,
	labels = TRUE
)