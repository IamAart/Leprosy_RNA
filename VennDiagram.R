library(VennDiagram)
library (purrr)
library(stringr)
library(readxl)

source("common.R")

# Get Gene names or Other column info (Ensemble names) from the saved data
limma <- readxl::read_excel("./data/LimmaVoom/data-First-HHC.xlsx")[[VENN_NAMES]]
edger <- readxl::read_excel("./data/EdgeR/data-First-HHC.xlsx")[[VENN_NAMES]]
deseq <- readxl::read_excel("./data/DESeq2/data-First-HHC.xlsx")[[VENN_NAMES]]

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
overlaps <- calculate.overlap(list(limma, edger, deseq))
indx <- as.numeric(substr(names(overlaps),2,2))
for (i in seq_along(overlaps)){
  v[[6 + indx[i] ]]$label <- paste(overlaps[[i]], collapse = "\n")
}

# Draw Venn Diagram
grid.newpage()
grid.draw(v)