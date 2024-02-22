# Loading R packages
library(limma)
library(edgeR)
library(Glimma)
library(dplyr)

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

# Create Contrast (HHC vs First)
design <- model.matrix(~0+factor(samples$group, c("HHC", "First", "Second")))
colnames(design) <- c("HHC", "First", "Second")
contrast <- makeContrasts(First_vs_HHC = First - HHC, levels = design)

# Normalization
normalized_data <- normalization_edgeR(counts, samples, design)

# Perform Limma Voom Analysis
voom <- voom(normalized_data, design, plot=TRUE, normalize="quantile")
fit <- lmFit(voom, design)
fit <- contrasts.fit(fit, contrasts = contrast)
ebayes <- eBayes(fit)
res <-  topTable(ebayes, n=Inf, adjust.method = "BH")

# Showcase Results with Glimma
glimmaMA(ebayes)

# Save data in table with cut off added to cut the data at that p value
make_save_table(res, "adj.P.Val", CUT_OFF, "P.Value", features, "./data/LimmaVoom/all_nc_data-First-HHC")