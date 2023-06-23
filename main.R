library(DESeq2)
library(ggplot2)
library(edgeR)
library(Glimma)

# GLOBAL VARIABLES
COUNT_DATA_FILENAME <- "./data/sasc326.rda"
COUNTS_VARIABLE_NAME <- "countsMatrix"
SAMPLES_VARIABLE_NAME <- "samplesData"
FEATURES_VARIABLE_NAME <- "featuresData"
DESIGN <- ~ groupFirst + groupHHC + type + YOB + gender


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
    return(list(counts = counts, features = features))
}

remove_rowData_dup_colnames <- function(data) {
    # rowData cannot have columns named "seqnames", "ranges", "strand", "start", "end", "width", "element"
    forbidden_colnames <- list("seqnames", "ranges", "strand", "start", "end", "width", "element")
    replace_colnames <- list("sequence_names", "seq_ranges", "seq_strand", "start_loc", "end_loc", "width_size", "seq_element")
    for (i in len(forbidden_colnames)) {
        names(data)[names(data) == forbidden_colnames[i]] <- replace_colnames[i]
    }
    return(data)
}

prepare_DGE_data <- function(data, dge_type, adjust_common, include_features, c_name, f_name, s_name) {
    # prepare the dataset with the features column as rownames in order to match it with features dataset
    counts <- column_to_rownames(data[[c_name]], 1)
    features <- column_to_rownames(data[[f_name]], 1)

    # find the genes with the same name in both feature and count data
    list[counts, features] <- equalize_feature_counts(counts, features, adjust_common)
    if (dge_type == "Limma-Voom" || dge_type == "EdgeR") {
        if (include_features) {
            return(DGEList(counts = counts, samples = data[[s_name]], genes = features))
        } else {
            return(DGEList(counts = counts, samples = data[[s_name]]))
        }
    } else if (dge_type == "DESeq2") {
        if (include_features) {
            features <- remove_rowData_dup_colnames(features)
            return(DESeqDataSetFromMatrix(countData = counts, colData = data[[s_name]], rowData = features, design = DESIGN))
            return(dds)
        } else {
            return(DESeqDataSetFromMatrix(countData = counts, colData = data[[SAMPLES_VARIABLE_NAME]], design = DESIGN))

        }
    } else {
        stop("dge_type parameter should be either 'Limma-Voom', 'EdgeR' or 'DESeq2'")
    }
}

perform_DGE_analysis <- function(data, dge_type, adjust_common, include_features, c_name, f_name, s_name) {
    dge_data <- prepare_DGE_data(data, dge_type, adjust_common, include_features, c_name, f_name, s_name)
    # TODO: return data for table, plot and/or more + save table data in file
    # TODO: Adjust contrast within each function, maybe add the contrast by variable
    # TODO: Add this or update this in the Limma-Voom and EdgeR procedure. It should come between dge_data and the first sentence in the if statement
    # dge <- calcNormFactors(dge_list)
    #
    # keep <- filterByExpr(y = dge)
    # dge <- dge[keep, , keep.lib.sizes = FALSE]
    # dge
    if (dge_type == "DESeq2") {
        contr <- c("HHC vs First", "groupHHC", "groupFirst") # TODO: Fix the contrast
        deseq <- DESeq(dge_data)
        deseq_results <- results(deseq, contrast = contr, alpha = 0.05)
        return(deseq_results)
    } else if (dge_type == "Limma-Voom") {
        design <- model.matrix(DESIGN, data[[s_name]])
        v <- voom(dge_data, design, plot=TRUE)
        fit <- lmFit(v, design)
        # Make contrast here
        tmp <- eBayes(fit, robust = TRUE)
        return(tmp)
    } else if (dge_type == "EdgeR") {
        dge <- calcNormFactors(dge_data)
        design <- model.matrix(DESIGN, count_data[[s_name]])
        rownames(design) <- colnames(dge)
        disp <- estimateDisp(dge, design)
        et <- exactTest(disp)
        return(et)
    } else {
        stop("dge_type parameter should be either 'Limma-Voom', 'EdgeR' or 'DESeq2'")
    }
}

count_data <- load_one_rdata_object(COUNT_DATA_FILENAME)
count_data[[SAMPLES_VARIABLE_NAME]]['groupFirst'] <- ifelse(count_data[[SAMPLES_VARIABLE_NAME]]$group=="First",1,0)
count_data[[SAMPLES_VARIABLE_NAME]]['groupSecond'] <- ifelse(count_data[[SAMPLES_VARIABLE_NAME]]$group=="Second",1,0)
count_data[[SAMPLES_VARIABLE_NAME]]['groupHHC'] <- ifelse(count_data[[SAMPLES_VARIABLE_NAME]]$group=="HHC",1,0)
count_data[[SAMPLES_VARIABLE_NAME]]['typeMB'] <- ifelse(count_data[[SAMPLES_VARIABLE_NAME]]$type=="MB",1,0)
count_data[[SAMPLES_VARIABLE_NAME]]['typePB'] <- ifelse(count_data[[SAMPLES_VARIABLE_NAME]]$type=="PB",1,0)
count_data[[SAMPLES_VARIABLE_NAME]]['typeHHC'] <- ifelse(count_data[[SAMPLES_VARIABLE_NAME]]$type=="HHC",1,0)

# edgeR
# dge_result <- create_edgeR_result(count_data, FALSE, FALSE)

# Limma Voom
# dge_result <- create_limma_voom_result(count_data, TRUE, TRUE)

#DESeq2
# dge_result <- create_deseq_result(count_data, TRUE, TRUE)

# TODO: Showcase results in the same way (thinking Heatmap and Histogram)
# TODO: Pick column DGE you would like to see plotted
# glimmaMA(dge_result)

# showcase results in a Volcano plot
# glimmaVolcano(dge_result)

# TODO: Optimize Design
# TODO: Improve normalization and preprocessing