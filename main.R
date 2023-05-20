library(DESeq2)
library(ggplot2)
library(edgeR)

# GLOBAL VARIABLES
COUNT_DATA_FILENAME <- "data/sasc326.rda"
COUNTS_VARIABLE_NAME <- "countsMatrix"
SAMPLES_VARIABLE_NAME <- "samplesData"
FEATURES_VARIABLE_NAME <- "featuresData"
# DESIGN <- ~ type
DESIGN <- ~group + YOB + gender

load_one_rdata_object <- function(file_path) {
    # Gets one object from rdata (in our case "gte") and returns it
    # Benefit: variable name can be defined and used easily
    res <- local({
        load(file_path)
        return(get(ls()))
    })
    return(res)
}

create_deseq_result <- function(count_data, design, common_bool, feat_bool) {
    # prepare the dataset with the features column as rownames in order to match it with features dataset
    counts <- column_to_rownames(count_data[[COUNTS_VARIABLE_NAME]], 1)
    features <- column_to_rownames(count_data[[FEATURES_VARIABLE_NAME]], 1)

    # find the genes with the same name in both feature and count data
    if (common_bool) {
        common <- find_common_rownames(list(counts, features))
        counts <- counts[common, ]
        features <- features[common, ]
    }

    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = count_data[[SAMPLES_VARIABLE_NAME]],
        design = design
    )

    if (feat_bool) {
        # mcols(dds) cannot have columns named "seqnames", "ranges", "strand", "start", "end", "width", "element"
        names(features)[names(features) == 'start'] <- 'start_loc'
        names(features)[names(features) == 'end'] <- 'end_loc'
        mcols(dds) <- cbind(mcols(dds), features)
    }

    deseq <- DESeq(dds)
    deseq
}

prepare_dge_list <- function(count_data, common_bool, feat_bool) {
    # prepare the dataset with the features column as rownames in order to match it with features dataset
    counts <- column_to_rownames(count_data[[COUNTS_VARIABLE_NAME]], 1)
    features <- column_to_rownames(count_data[[FEATURES_VARIABLE_NAME]], 1)

    # find the genes with the same name in both feature and count data
    if (common_bool) {
        common <- find_common_rownames(list(counts, features))
        counts <- counts[common, ]
        features <- features[common, ]
    }

    if (feat_bool) {
        dge_list <- DGEList(
            counts = counts,
            samples = count_data[[SAMPLES_VARIABLE_NAME]],
            genes = features
        )
    }
    else {
        dge_list <- DGEList(
            counts = counts,
            samples = count_data[[SAMPLES_VARIABLE_NAME]],
        )
    }

    dge <- calcNormFactors(dge_list)

    keep <- filterByExpr(y = dge)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge
}

create_edgeR_result <- function(count_data, common, feat) {
    dge <- prepare_dge_list(count_data, common, feat)
    #normalization
    dge <- calcNormFactors(dge)
    # estimating dispersions with qCML
    disp <- estimateDisp(dge)

    # Check DGE
    et <- exactTest(disp)
    et

}

create_limma_voom_result <- function(count_data, common, feat) {
    dge <- prepare_dge_list(count_data, common, feat)

    design <- model.matrix(DESIGN, count_data[[SAMPLES_VARIABLE_NAME]])
    v <- voom(dge, design, plot=TRUE)
    vfit <- lmFit(v, design)
    tfit <- treat(vfit, lfc=1)
    tfit

}

find_common_rownames <- function(list_of_df) {
    Reduce(intersect, lapply(list_of_df, row.names))
}

column_to_rownames <- function(df, column_number) {
    rownames(df) <- df[,column_number]
    df[,column_number] <- NULL
    df
}

count_data <- load_one_rdata_object(COUNT_DATA_FILENAME)

# edgeR
# plotMD(create_loom_result(count_data))

# Limma Voom
# plotMD(create_limma_voom_result(count_data))

#DESeq2
DESeq2::plotMA(results(create_deseq_result(count_data, DESIGN, TRUE, TRUE)))
