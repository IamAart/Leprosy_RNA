library("writexl")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("tibble")
library("stringr")
library("dotenv")

dotenv::load_dot_env(file = ".env")
NON_CODING_BIOTYPES <- as.vector(strsplit(Sys.getenv("NON_CODING_BIOTYPES"), ","))
CODING_BIOTYPES <- as.vector(strsplit(Sys.getenv("CODING_BIOTYPES"), ","))
CUT_OFF <- as.numeric(Sys.getenv("CUT_OFF"))
LOG2_OFF <- log2(as.numeric(Sys.getenv("LOG2_OFF")))

FILEPATH <- Sys.getenv("FILEPATH")
COUNTS_NAME <- Sys.getenv("COUNTS_NAME")
FEATURES_NAME <- Sys.getenv("FEATURES_NAME")
SAMPLES_NAME <- Sys.getenv("SAMPLES_NAME")

DGE_METHOD <- Sys.getenv("DGE_METHOD")
BOOL_MDS_PLOT <- as.logical(Sys.getenv("BOOL_MDS_PLOT"))

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

column_to_rownames <- function(df, column_name) {
    # Function to swap the first column to the rownames
    rownames(df) <- df[,column_name]
    df[,column_name] <- NULL
    df
}

check_order_rownames <- function (data1, data2) {
    # Check whether the rownames of the samples are in the same order as the columns of the counts.
    # If not adjust so that counts columns are adjusted to the standard of the rownames in the samples
    # If now all rownames of the counts exist in the colnames of the counts return false which triggers an error
    if (!all(rownames(data1) %in% colnames(data2))) {
        return(FALSE)
    }
    if (!all(rownames(data1) == colnames(data2))) {
        data2 <- data2[, rownames(data1)]
    }
    return(data2)
}

prepare_data <- function(path, counts_name, features_name, samples_name) {
    # Prepare data so that it easy to use during normalization
    # It adjusts the rownames to the column containing the ensemble ids or the samples to the rownames. It also selects all non-coding genes
    # Returns a list containing the 3 data structures counts, samples and features
    data <- load_one_rdata_object(path)
    counts <- column_to_rownames(data[[counts_name]], "feature")
    features <- column_to_rownames(data[[features_name]], "feature")
    counts <- merge(counts, features, by = 'row.names')
    counts <- column_to_rownames(counts, "Row.names")

    samples <- data[[samples_name]]
    samples <- tibble::column_to_rownames(samples, "sample")
    return(list(counts = counts, features = features, samples = samples))
}

make_mds_plot <- function(norm_data, plot_path) {
    colors <- c("#88ccee", "#44aa99", "#117733", "#332288", "#ddcc77", "#999933", "#cc6677", "#882255", "#aa4499", "#dddddd")

    plot <- plotMDS(norm_data)
    plot <- as.data.frame(plot)
    plot_data <- plot[, c("x", "y")]
    plot_data["Row.names"] <- rownames(plot_data)
    samples["Row.names"] <- rownames(samples)

    hhc <- samples[samples$group == "HHC", ]$Row.names
    t1 <- samples[samples$group == "First", ]$Row.names
    data <- plot_data %>%
        mutate( Colorcode = case_when(
            Row.names %in% hhc ~ "HHC",
            Row.names %in% t1 ~ "First",
            TRUE ~ "Second"
        ) )
    data$Colorcode <- factor(data$Colorcode, levels = c("HHC", "First", "Second"))

    text_repel <- data
    # text_repel <- data[data$x <= -1.25, ]
    # text_repel <- data[data$Row.names %in% c("s103830.003.011", "s103830.004.017", "s103830.003.003"), ]
    p <- ggplot(data, aes(x=x, y=y, label = Row.names, color=Colorcode)) +
            geom_point(size=3) +
            geom_text_repel(data = text_repel, size=7) +
            theme_minimal() +
            labs(x = "Log2 Fold Change Dimension 1", y="Log 2 Fold Change Dimension 2") +
            scale_color_manual(
                values = c(HHC = colors[1], First = colors[4], Second = colors[7]),
                labels = c(HHC = "Healthy", First = "Progressors t=1", Second = "Progressors t=2")
            ) +
            guides(color=guide_legend("Patient Information")) +
            theme(
                # Axis titles font size
                axis.title.x = element_text(size = 19),    # x-axis label size
                axis.title.y = element_text(size = 19),    # y-axis label size

                # Axis tick labels font size
                axis.text.x = element_text(size = 17),     # x-axis ticks (labels) size
                axis.text.y = element_text(size = 17),     # y-axis ticks (labels) size

                # Legend title and text font size
                legend.title = element_text(size = 17),    # Legend title size
                legend.text = element_text(size = 16)      # Legend labels size
            )
    ggsave(plot_path, p, width=15, height=10, dpi=720)
}

normalization_edgeR <- function(counts, samples, design, features) {
    # Normalize data from of the counts and sample and put them into a DGE List before going into the DGE Analysis
    normalized_data <- DGEList(counts = counts, samples = samples)

    # This reduces the amount of rows to 17.000, however this in order with the EdgeR vigenette
    keep <- filterByExpr(y = normalized_data, design)
    normalized_data <- normalized_data[keep, ,keep.lib.sizes = FALSE]

    normalized_data <- calcNormFactors(normalized_data, method = "TMM")
    return(normalized_data)
}

make_save_table <- function(result, features, filepath_name) {
    # Saves the result of all three libraries to an excel file and a csv file
    data <- merge(as.data.frame(result), features, by="row.names")
    names(data)[names(data) %in% c("padj", "FDR", "adj.P.Val")] <- "P.Adjust" # DESeq2, EdgeR, LimmaVoom
    names(data)[names(data) %in% c("log2FoldChange", "logFC")] <- "Log2FoldChange" # DESeq2, EdgeR & LimmaVoom
    data <- data[order(data$P.Adjust), ]
    data <- data %>%
        mutate( col = case_when(
            P.Adjust <= CUT_OFF & Log2FoldChange >= LOG2_OFF ~ "up",
            P.Adjust <= CUT_OFF & Log2FoldChange < -LOG2_OFF ~ "down",
            TRUE ~ "not"
        ) )

    write.csv(data, file = paste0(filepath_name, ".csv"))
    writexl::write_xlsx(data, paste0(filepath_name, ".xlsx"))
    return(data)
}