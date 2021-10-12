## compare results DE
# libraries and setup
library(here) # paths reproducible
library(edgeR) # DE import from edgeR
library(dplyr) #manipulatedf
project <- "fire-mice"
# import
de_results_edgeR <- readRDS(here("processed", project, "DE_oligo_edgeR_de_results.RDS"))
de_results_edgeR_dr <- readRDS(here("processed", project, "DE_oligo_edgeR_detRate_de_results.RDS"))
top_mast <- readRDS(here("processed", project, "DE_oligo_mast_de_results.RDS"))
top_dr <- data.frame(topTags(de_results_edgeR_dr, n = 120))
top_edge <- data.frame(topTags(de_results_edgeR, n = 120))

# similar
# edge_dr
sum(row.names(top_dr) %in% row.names(top_edge))
sum(row.names(top_dr) %in% row.names(top_mast))
#mast
row.names(top_mast)[row.names(top_mast) %in% row.names(top_edge)]
row.names(top_mast)[row.names(top_mast) %in% row.names(top_dr)]
# edge
(row.names(top_edge) %in% row.names(top_dr))
(row.names(top_edge) %in% row.names(top_mast))