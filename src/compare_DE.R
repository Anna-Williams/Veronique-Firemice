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
top_mast_batch <- readRDS(here("processed", project, "DE_oligo_mast_batchcorrect_de_results.RDS"))
top_edge_dr <- data.frame(topTags(de_results_edgeR_dr, n = 120))
top_edge <- data.frame(topTags(de_results_edgeR, n = 120))

# extract just the gene names, and add helpful names
# b = batch correct
# dr = detection rate correct

 mast_dr <- row.names(top_mast)
 mast_dr_b <- row.names(top_mast_batch)
 edgeR_b <- row.names(top_edge)
 edgeR_dr_b <- row.names(top_edge_dr)
 
# edge_dr
# compared with the other edgeR
sum(edgeR_dr_b %in% edgeR_b)
edgeR_dr_b[edgeR_dr_b %in% edgeR_b]
# compared with the masts
sum(edgeR_dr_b %in% mast_dr_b)
edgeR_dr_b[edgeR_dr_b %in% mast_dr_b]
sum(edgeR_dr_b %in% mast_dr)


#mast
mast_dr[row.names(top_mast) %in% row.names(top_edge)]
row.names(mast_dr)[row.names(top_mast) %in% row.names(top_dr)]
# edge
(row.names(top_edge) %in% row.names(top_dr))
(row.names(top_edge) %in% row.names(mast_dr))