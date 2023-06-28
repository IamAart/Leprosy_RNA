library(DESeq2)
library(ggplot2)
library(edgeR)
library(tibble)
library(ggplot2)
library(ggrepel)

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
    # Make feature contain the same genes as in the counts table, if the gene does not exist make it an NA row
    if (adjust_common) {
        counts_rows <- rownames(counts)[!(rownames(counts) %in%  rownames(features))]
        features <- features[rownames(counts) %in% rownames(features),]
        features[counts_rows, ] <- NA
        return(list(counts = counts, features = features))
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
    samples <- data[[s_name]]
    if (!all(rownames(samples) %in% colnames(counts))) {
        samples <- tibble::column_to_rownames(samples, "sample")
    }
    if (!all(rownames(samples) == colnames(counts))) {
        counts <- counts[, rownames(samples)]
    }

    # find the genes with the same name in both feature and count data
    c_f <- equalize_feature_counts(counts, features, adjust_common)
    counts <- c_f[["counts"]]
    features <- c_f[["features"]]
    if (!all(rownames(features) == rownames(counts)) && adjust_common) {
        features <- features[rownames(counts), ]
    }
    if (dge_type == "LimmaVoom" || dge_type == "EdgeR") {
        if (include_features) {
            return(DGEList(counts = counts, samples = samples, genes = features))
        } else {
            return(dge = DGEList(counts = counts, samples = samples))
        }
    } else if (dge_type == "DESeq2") {
	    dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = design)
        if (include_features) {
            features <- remove_rowData_dup_colnames(features)
            mcols(dds) <- DataFrame(mcols(dds), features)
        }
        return(list(dge = dds, genes = features))
    } else {
        stop("dge_type parameter should be either 'LimmaVoom', 'EdgeR' or 'DESeq2'")
    }
}

perform_DGE_analysis <- function(data, dge_type, adjust_common, include_features, design, c_name, f_name, s_name) {
    dge_data <- prepare_DGE_data(data, dge_type, adjust_common, include_features, design, c_name, f_name, s_name)

    if (dge_type == "DESeq2") {
	    if (include_features) {
		    deseq <- DESeq(dge_data$dge)
		    result <- results(deseq, contrast = c("group", "First", "HHC"), alpha = 0.05)
	        return(list(result = result, genes = dge_data$genes))
	    } else {
		    deseq <- DESeq(dge_data)
		    result <- results(deseq, contrast = c("group", "First", "HHC"), alpha = 0.05)
	        return(result)
	    }
    } else {
        group <- factor(data[[s_name]]$group, levels = c("HHC", "First", "Second"))
        design <- model.matrix(~ group + YOB + gender, data[[s_name]])
        colnames(design) <- c(levels(data[[s_name]]$group), "YOB", "gender")
        dge_data <- calcNormFactors(dge_data)

        keep <- filterByExpr(y = dge_data, design)
        dge_data <- dge_data[keep, , keep.lib.sizes = FALSE]

        if (dge_type == "LimmaVoom") {
            v <- voom(dge_data, design, plot=TRUE)
            fit <- lmFit(v, design)
            contrast <- makeContrasts(First_vs_HHC = First - HHC, levels = design)
            fit <- contrasts.fit(fit, contrasts = contrast)
            result <- eBayes(fit)
            return(list(result = result, genes = result$genes))
        } else if (dge_type == "EdgeR") {
            disp <- estimateDisp(dge_data, design)
            result <- exactTest(disp, pair= c("HHC", "First"))
            return(list(result = result, genes = result$genes))
        } else {
            stop("dge_type parameter should be either 'LimmaVoom', 'EdgeR' or 'DESeq2'")
        }
    }
}

volcano_dge <- function(data, gene_data, dge_type, p_lim, log_lim, comparison, labels) {
    if (dge_type == "DESeq2") {
        log2foldchange_name <- "log2FoldChange"
        pvalue_name <- "pvalue"
    } else if (dge_type == "LimmaVoom") {
        log2foldchange_name <- "logFC"
        pvalue_name <- "P.Value"
    } else if (dge_type == "EdgeR") {
        log2foldchange_name <- "logFC"
        pvalue_name <- "PValue"
    } else {
        stop("dge_type parameter should be either 'LimmaVoom', 'EdgeR' or 'DESeq2'")
    }

    data$diffexpressed <- "NO"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
    data$diffexpressed[data[[log2foldchange_name]] > log_lim & data[[pvalue_name]] < p_lim] <- "UP"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    data$diffexpressed[data[[log2foldchange_name]] < -log_lim & data[[pvalue_name]]] <- "DOWN"


    data$delabel <- NA
	if (labels) {
		data$delabel[which(data$diffexpressed != "NO")] <- gene_data$gene_name[which(data$diffexpressed != "NO")]
	}

    plot <- ggplot(data=data, aes(x=data[[log2foldchange_name]], y=-log10(data[[pvalue_name]]), col=diffexpressed, label=delabel)) +
        geom_point() +
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-log_lim, log_lim), col="red") +
        geom_hline(yintercept=-log10(p_lim), col="red")

    ggsave(
        filename = paste0("plot-", comparison, ".png"),
        plot = plot,
        path = paste0("./data/", dge_type, "/")
    )
}

save_and_show <- function(data, genes, dge_type, comparison, p_val, log_val, labels) {
    filepath_data <- paste0("./data/", dge_type, "/data-", comparison, ".csv")
	write.csv(data, file = filepath_data)
	filepath_genes <- paste0("./data/", dge_type, "/genes-", comparison, ".csv")
    write.csv(genes, file = filepath_genes)

	volcano_dge(
	    data,
	    genes,
	    dge_type,
	    p_val,
	    log_val,
	    comparison,
	    labels
	)
}
