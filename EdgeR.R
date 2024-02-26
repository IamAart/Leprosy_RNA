# Loading R packages
library(edgeR)

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

# Create constrast parameter (HHC vs First)
design <- model.matrix(~0+factor(samples$group, c("HHC", "First", "Second")))
colnames(design) <- c("HHC", "First", "Second")
contrast <- makeContrasts(First_vs_HHC = First - HHC, levels = design)

# Normalize data
normalized_data <- normalization_edgeR(counts, samples, design)

# Perform EdgeR Analysis
disp <- estimateGLMCommonDisp(normalized_data,design)
disp <- estimateGLMTrendedDisp(disp, design, method="power")
disp <- estimateGLMTagwiseDisp(disp, design)
fit <- glmFit(disp, design)
lrt <- glmLRT(fit, contrast=contrast)
res<- as.data.frame(topTags(lrt, n=Inf, adjust.method = "BH"))

# Save data in table with cut off added to cut the data at that p value
res <- make_save_table(res, features, "./data/EdgeR/NON_CODING_EdgeR_RNA_SEQ")

# Showcase Result in Volcano Plot
plot <- VolcanoPlot(res)
ggplot_build(plot)