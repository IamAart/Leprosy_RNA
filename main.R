library(DESeq2)
library(ggplot2)
library(edgeR)

# GLOBAL VARIABLES
COUNT_DATA_FILENAME <- "data/sasc326.rda"
COUNTS_VARIABLE_NAME <- "countsMatrix"
SAMPLES_VARIABLE_NAME <- "samplesData"
FEATURES_VARIABLE_NAME <- "featuresData"
# DESIGN <- ~ type
DESIGN <- ~ group + YOB + gender

load_one_rdata_object <- function(file_path) {
    # Gets one object from rdata (in our case "gte") and returns it
    # Benefit: variable name can be defined and used easily
    res <- local({
        load(file_path)
        return(get(ls()))
    })
    return(res)
}

create_deseq_result <- function(count_data, design) {
    # prepare the dataset with the features column as rownames in order to match it with features dataset
    rownames(count_data[[COUNTS_VARIABLE_NAME]]) <- count_data[[COUNTS_VARIABLE_NAME]][,1]
    count_data[[COUNTS_VARIABLE_NAME]][,1] <- NULL

    dds <- DESeqDataSetFromMatrix(
        countData = count_data[[COUNTS_VARIABLE_NAME]],
        colData = count_data[[SAMPLES_VARIABLE_NAME]],
        design = design
    )
    # TODO: Prepare Feature table to be same row amount as Counts Table
    # Currently feature data is not complete
    # rownames(count_data[[FEATURES_VARIABLE_NAME]]) <- count_data[[FEATURES_VARIABLE_NAME]][,1]
    # filtered_feature_data <- count_data[[FEATURES_VARIABLE_NAME]][rownames(count_data[[COUNTS_VARIABLE_NAME]]),]
    # apply feature data in the DESeqDataSet
    # mcols(dds) <- DataFrame(mcols(dds), count_data[[FEATURES_VARIABLE_NAME]])
    deseq <- DESeq(dds)
    deseq
}

create_edgeR_result <- function(count_data) {
    # put gene name into rownames
    rownames(count_data[[COUNTS_VARIABLE_NAME]]) <- count_data[[COUNTS_VARIABLE_NAME]][,1]
    count_data[[COUNTS_VARIABLE_NAME]][,1] <- NULL

    dge_list <- DGEList(
        counts = count_data[[COUNTS_VARIABLE_NAME]],
        samples = count_data[[SAMPLES_VARIABLE_NAME]],
        # genes = count_data[[FEATURES_VARIABLE_NAME]]
    )
    # filter by genes with low counts
    keep <- filterByExpr(y = dge_list)
    dge <- dge_list[keep, , keep.lib.sizes = FALSE]
    #normalization
    dge <- calcNormFactors(dge)
    # estimating dispersions with qCML
    disp <- estimateDisp(dge)

    # Check DGE
    et <- exactTest(disp)
    et
    # top_dge <- topTags(et, n="Inf")
    # print(summary(decideTests(object = et, lfc = 1)))

}

count_data <- load_one_rdata_object(COUNT_DATA_FILENAME)

# edgeR
# plotMD(create_loom_result(count_data))

# DESeq2
# DESeq2::plotMA(results(create_deseq_result(count_data, DESIGN)))
