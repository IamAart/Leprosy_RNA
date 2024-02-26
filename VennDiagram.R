# EVENTUALLY NOT USED DUE TO AMOUNT OF DGE GENES

library(VennDiagram)
library(stringr)

source("common.R")

# Get Gene names or Other column info (Ensemble names) from the saved data
limma <- read.csv("./data/LimmaVoom/all_nc_data-First-HHC.csv")[[VENN_NAMES]]
edger <- read.csv("./data/EdgeR/all_nc_data-First-HHC.csv")[[VENN_NAMES]]
deseq <- read.csv("./data/DESeq2/new_all_nc_data-First-HHC.csv")[[VENN_NAMES]]

print(length(unique(deseq[1:183])))

# Create Venn Diagram
v <- venn.diagram(
  x = list(limma, edger, deseq),
  category.names = c("Limma Voom" , "Edge R" , "DESeq2"),
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow"),
  disable.logging = TRUE,
  filename = NULL
)

# Create Names on Diagram
overlaps <- VennDiagram::calculate.overlap(list(limma, edger, deseq))
indx <- as.numeric(substr(names(overlaps),2,2))
for (i in seq_along(overlaps)){
  v[[6 + indx[i] ]]$label <- paste(overlaps[[i]], collapse = "\n")
}

# Draw Venn Diagram
grid.newpage()
grid.draw(v)