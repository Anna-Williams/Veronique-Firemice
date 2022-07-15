library(SingleCellExperiment)
library(ShinyCell)
library(here)
library(stringr) # modify metadata

project <- "fire-mice" 

sce <- readRDS(here("processed", project, "sce_oligo_clusters_01.RDS"))

# remove unwanted dimensional reductions
 reducedDim(sce, "PCA_all") <- NULL
 reducedDim(sce, "PCA_coldata") <- NULL


conf <- createConfig(sce)
#Delete some unnecessary metadata
conf <- delMeta(conf, meta.to.del = c( "subsets_mt_sum", "subsets_mt_detected",
                       "outlier", "ratio_detected_sum", "outlier_ratio", "chip", "age",
                       "discard", "sizeFactor", "total", "discard", "original_sample_name",
                       "cluster_k20", "cluster_k40", "cluster_k60", "cluster_k80", "cluster_oligo_k80" )
                )
clusters_oligo <- grep("cluster_oligo", conf$ID, value = TRUE)
clusters <- str_replace(clusters_oligo,"cluster_oligo", "clusters")

# Change name of some metadata
conf <- modMetaName(conf,
                    meta.to.mod = c("sum", "detected","Sample", "subsets_mt_percent", clusters_oligo),
                    new.name = c("umi counts", "detected genes", "sample", "% mt genes", clusters))
# reorder metadata
conf <- reorderMeta(conf, conf$ID[c(5,13,3,2,4,1,6:12)])
conf <- modDefault(conf, "genotype", "cluster_oligo")

# change colours
source(here("src/colours.R"))
conf <- modColours(conf, meta.to.mod = "genotype", 
                   new.colours = col_magenta_green)



metadata_categorical <- c("Sample", "tissue", "mouse", "batch", clusters_oligo)
for(metadata in metadata_categorical){
  print(metadata)
  n_levels <- length(levels(as.factor(sce[[metadata]])))
  conf <- modColours(conf, meta.to.mod = metadata,
                     new.colours = cols[1:n_levels])
}


# citation
citation = list(
  author  = "author",
  title   = "title",
  journal = "journal",
  volume  = "vol",
  page    = "page",
  year    = "year", 
  doi     = "doi",
  link    = "link")

makeShinyApp(sce, conf, gene.mapping = TRUE, 
             shiny.title = "FIRE mice Oligodendrocytes",
             default.dimred = c("TSNE1", "TSNE2"),
             default.gene1 = "Serpina3n",
             default.gene2 = "C4b",
             default.multigene = c( "Tagln2", "C4b", "Cldn11", "Flnc", "Snx10", "Bcl3", "Emp3", "Plvap", "Anxa2", "Trf"),
             shiny.footnotes = citation )



