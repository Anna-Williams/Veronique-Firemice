# Author:NadineBestard
# Date 26/02/21

#Here I import the count matrices generated by CellRanger into R
#using Bioconductor packages. I will use the filtered matrices where an initial
#cellcalling has been done by CellRanger. In Cellranger version 5 the cellcalling
#is done with a very similar algorithm to dropletutilities one (ref).


# Set up ------
library(here)
library(DropletUtils)

# Create file path for each sample matrix and metadata ------

# Select only the young  samples, WT and KO
metadata <- read.csv(here("data/metadata_scRNAseq.csv"))
metadata <- metadata[metadata$age == "young" & metadata$genotype %in% c("WT", "KO"), ]
matrices <- here("data/filtered/",  metadata$Sample, "filtered_feature_bc_matrix")


# Create object -----
sce <- read10xCounts(matrices, sample.names = metadata$Sample, version = "auto", col.names = TRUE)
# add the metadata
coldata <- colData(sce)
# keep the rownames (same as colnames of sce)
names <- row.names(coldata)
coldata <- merge(x= coldata, y= metadata, by = "Sample", all.x=TRUE)
# reset the rownames to previous value, as there're lost with merge
row.names(coldata) <- names
# assign new coldata back to sce 
colData(sce) <- coldata
# The chip is detected as an integer instead of a character
sce$chip <- as.character(sce$chip)
# Save the object
saveRDS(sce, here("processed", "sce_young.RDS"))


