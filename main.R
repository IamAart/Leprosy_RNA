library(DESeq2)
library(ggplot2)
library(edgeR)
library(Glimma)

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
    rownames(df) <- df[,column_number]
    df[,column_number] <- NULL
    df
}

equalize_feature_counts <- function(counts, features, adjust_common) {
    # find the genes with the same name in both feature and count data and remove the rest
    if (adjust_common) {
        common <- Reduce(intersect, lapply((list(counts, features)), row.names))
        return(list(counts = counts[common, ], features = features[common, ]))
    }
    return(c(counts = counts, features = features))
}

remove_rowData_dup_colnames <- function(data) {
    # rowData cannot have columns named "seqnames", "ranges", "strand", "start", "end", "width", "element"
    forbidden_colnames <- c("seqnames", "ranges", "strand", "start", "end", "width", "element")
    replace_colnames <- c("sequence_names", "seq_ranges", "seq_strand", "start_loc", "end_loc", "width_size", "seq_element")
    for (i in seq_along(forbidden_colnames)) {
        names(data)[names(data) == forbidden_colnames[[i]]] <- replace_colnames[[i]]
    }
    return(data)
}

prepare_DGE_data <- function(data, dge_type, adjust_common, include_features, design, c_name, f_name, s_name) {
    # prepare the dataset with the features column as rownames in order to match it with features dataset
    counts <- column_to_rownames(data[[c_name]], 1)
    features <- column_to_rownames(data[[f_name]], 1)

    # find the genes with the same name in both feature and count data
    c_f <- equalize_feature_counts(counts, features, adjust_common)
    counts <- c_f[["counts"]]
    features <- c_f[["features"]]
    if (dge_type == "Limma-Voom" || dge_type == "EdgeR") {
        if (include_features) {
            return(DGEList(counts = counts, samples = data[[s_name]], genes = features))
        } else {
            return(DGEList(counts = counts, samples = data[[s_name]]))
        }
    } else if (dge_type == "DESeq2") {
        if (include_features) {
            features <- remove_rowData_dup_colnames(features)
            return(DESeqDataSetFromMatrix(countData = counts, colData = data[[s_name]], rowData = features, design = design))
            return(dds)
        } else {
            return(DESeqDataSetFromMatrix(countData = counts, colData = data[[SAMPLES_VARIABLE_NAME]], design = design))

        }
    } else {
        stop("dge_type parameter should be either 'Limma-Voom', 'EdgeR' or 'DESeq2'")
    }
}

perform_DGE_analysis <- function(data, dge_type, adjust_common, include_features, design, c_name, f_name, s_name) {
    dge_data <- prepare_DGE_data(data, dge_type, adjust_common, include_features, design, c_name, f_name, s_name)
    # TODO: return data for table, plot and/or more + save table data in file

    if (dge_type == "DESeq2") {
        contr <- c("group", "HHC", "First") # This is used to showcase the result of the difference in gene expression between First and HHC
        deseq <- DESeq(dge_data)
        result <- results(deseq, contrast = contr, alpha = 0.05)
        return(result)
    } else {
        group <- factor(dge_data$group)
        design <- model.matrix(~ 0 + group, data[[s_name]])
        colnames(design) <- levels(data[[s_name]]$group)
        # design <- model.matrix(design, data[[s_name]])
        # colnames(design) <- c(levels(data[[s_name]]$group), "YOB", "gender")
        dge_data <- calcNormFactors(dge_data)
        # keep <- filterByExpr(y = dge_data, design)
        # dge_data <- dge_data[keep, , keep.lib.sizes = FALSE]
        if (dge_type == "Limma-Voom") {
            v <- voom(dge_data, design, plot=TRUE)
            fit <- lmFit(v, design)
            contrast <- makeContrasts(First_vs_HHC = HHC - First, levels = v$design)
            fit <- contrasts.fit(fit, contrasts = contrast)
            result <- eBayes(fit)
            return(result)
        } else if (dge_type == "EdgeR") {
            # TODO fix: contrast
            disp <- estimateDisp(dge_data, design)
            fit <- glmQLFit(disp, design)
            contrast <- makeContrasts(First_vs_HHC = HHC - First, levels = design)
            result <- glmQLFTest(fit, contrast = contrast[, "First_vs_HHC"])
            # et <- exactTest(disp)
            return(result)
        } else {
            stop("dge_type parameter should be either 'Limma-Voom', 'EdgeR' or 'DESeq2'")
        }
    }
}

COUNT_DATA_FILENAME <- "./data/sasc326.rda"
COUNTS_VARIABLE_NAME <- "countsMatrix"
SAMPLES_VARIABLE_NAME <- "samplesData"
FEATURES_VARIABLE_NAME <- "featuresData"
DGE_TYPE <- "EdgeR" # Can be "DESeq2", "Limma-Voom" or "EdgeR"
ADJUST_COMMON <- TRUE
INCLUDE_FEATURES <- TRUE

DESIGN <- ~ group + YOB + gender

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
# TODO: save data as csv
# TODO: showcase in different kind of figures: MA plot, Volcano Plot, Heatmap, p-value Histogram

# glimmaMA(dge_result)
# glimmaVolcano(dge_result)

# TODO: Optimize Design
# TODO: Improve normalization and preprocessing