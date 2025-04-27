library("limma")
library("edgeR")
library("DESeq2")
library("dplyr")
source("common.R")

.libPaths("./R_Packages/")

# These global variables are gathered from the .env file in the common.R file
data <- prepare_data(FILEPATH, COUNTS_NAME, FEATURES_NAME, SAMPLES_NAME)
counts <- data[["counts"]]
features <- data[["features"]]
samples <- data[["samples"]]

# Make sure the order of the rownames and colnames is the same
counts <- check_order_rownames(samples, counts)
if (is.logical(counts)) {
  stop("Error when adjusting colnames and rownames to be similar")
}

if (DGE_METHOD == "DESeq2") {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~ group)
  # Create MDS plot if selected
  if (BOOL_MDS_PLOT) {
    sf <- estimateSizeFactors(dds)
    normalized_data <- DESeq2::counts(sf, normalized=TRUE)
    normalized_data <- log2(normalized_data)
    make_mds_plot(normalized_data, sprintf("./data/Plots/%s_MDS_plot.png", DGE_METHOD))
  }
  # Perform DESEQ2
  deseq <- DESeq(dds)
  res <- results(deseq, c("group", "First", "HHC"), pAdjustMethod = "BH")
} else {
  # Create constrast parameter (HHC vs First)
  design <- model.matrix(~0+factor(samples$group, c("HHC", "First", "Second")))
  colnames(design) <- c("HHC", "First", "Second")
  contrast <- makeContrasts(First_vs_HHC = First - HHC, levels = design)

  # Normalize data
  normalized_data <- normalization_edgeR(counts, samples, design, features)
  # Create MDS plot if selected
  if (BOOL_MDS_PLOT) { make_mds_plot(normalized_data, sprintf("./data/Plots/%s_MDS_plot.png", DGE_METHOD)) }
  if (DGE_METHOD == "LimmaVoom") {
    # Prepare data for limma voom
    voom <- voom(normalized_data, design, plot=FALSE, normalize="quantile")
    fit <- lmFit(voom, design)
    fit <- contrasts.fit(fit, contrasts = contrast)
    # Perform DGE Analysis
    ebayes <- eBayes(fit)
    res <-  topTable(ebayes, n=Inf, adjust.method = "BH")
  }
  else if (DGE_METHOD == "EdgeR") {
    # Prepare data for edgeR
    disp <- estimateDisp(normalized_data, design)
    # disp <- estimateGLMCommonDisp(normalized_data,design)
    # disp <- estimateGLMTrendedDisp(disp, design, method="power")
    # disp <- estimateGLMTagwiseDisp(disp, design)
    fit <- glmFit(disp, design)
    # Perform DGE Analysis
    lrt <- glmLRT(fit, contrast=contrast)
    res<- as.data.frame(topTags(lrt, n=Inf, adjust.method = "BH"))
  }
}

# Save results in table
res <- make_save_table(res, features, sprintf("./data/%s/TEST_%s_RNA_SEQ", DGE_METHOD, DGE_METHOD))
