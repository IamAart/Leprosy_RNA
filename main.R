library(DESeq2)
library(ggplot2)

# GLOBAL VARIABLES
COUNT_DATA_FILENAME <- "data/sasc326.rda"
COUNTS_VARIABLE_NAME <- "countsMatrix"
SAMPLES_VARIABLE_NAME <- "samplesData"
FEATURES_VARIABLE_NAME <- "featuresData"
DESIGN <- ~ type

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
    result <- results(DESeq(dds))
}

count_data <- load_one_rdata_object(COUNT_DATA_FILENAME)
test <- create_deseq_result(count_data, DESIGN)

# showcase results
summary(test)
plotMA(test, ylim=c(-2,2))