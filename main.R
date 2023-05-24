library(DESeq2)
library(ggplot2)
library(edgeR)
library(Glimma)

# GLOBAL VARIABLES
COUNT_DATA_FILENAME <- "data/sasc326.rda"
COUNTS_VARIABLE_NAME <- "countsMatrix"
SAMPLES_VARIABLE_NAME <- "samplesData"
FEATURES_VARIABLE_NAME <- "featuresData"
# DESIGN <- ~ type
DESIGN <- ~ group

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

    if (feat_bool) {
        # rowData cannot have columns named "seqnames", "ranges", "strand", "start", "end", "width", "element"
        names(features)[names(features) == 'start'] <- 'start_loc'
        names(features)[names(features) == 'end'] <- 'end_loc'
        dds <- DESeqDataSetFromMatrix(
            countData = counts,
            colData = count_data[[SAMPLES_VARIABLE_NAME]],
            rowData = features,
            design = design
        )
    } else {
        dds <- DESeqDataSetFromMatrix(
            countData = counts,
            colData = count_data[[SAMPLES_VARIABLE_NAME]],
            design = design
        )
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
    print(et)
    et
}

create_limma_voom_result <- function(count_data, common, feat) {
    dge <- prepare_dge_list(count_data, common, feat)

    design <- model.matrix(DESIGN, count_data[[SAMPLES_VARIABLE_NAME]])
    v <- voom(dge, design, plot=TRUE)
    fit <- lmFit(v, design)
    # contr <- contrasts(group)
    # tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(fit)
    # tfit <- treat(vfit, lfc=1)
    # tfit
    tmp

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
dge_result <- create_edgeR_result(count_data, TRUE, TRUE)

# Limma Voom
# dge_result <- create_limma_voom_result(count_data, TRUE, TRUE)

#DESeq2
# dge_result <- create_deseq_result(count_data, DESIGN, TRUE, TRUE)

# TODO Showcase results in the same way (thinking Heatmap and Histogram)

# showcase results in a MA plot and Table
# glimmaMA(dge_result)

# showcase results in a Volcano plot
glimmaVolcano(dge_result)

#TODO optimize Design