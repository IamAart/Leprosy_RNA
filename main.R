library(DESeq2)

# GLOBAL VARIABLES
COUNT_DATA_FILENAME <- "sasc326.rda"

load_one_rdata_object <- function(file_path) {
    # Gets one object from rdata (in our case "gte") and returns it
    # Benefit: variable name can be defined and used easily
    res <- local({
        load(file_path)
        return(get(ls()))
    })
    return(res)
}

count_data <- load_one_rdata_object(COUNT_DATA_FILENAME)
