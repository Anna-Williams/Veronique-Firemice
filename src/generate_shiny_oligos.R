library(SingleCellExperiment)
library(ShinyCell)
library(here)

project <- "fire-mice" 

sce <- readRDS(here("processed", project, "sce_oligo_clusters_01.RDS"))

# remove unwanted dimensional reductions
 reducedDim(sce, "PCA_all") <- NULL
 reducedDim(sce, "PCA_coldata") <- NULL


conf <- createConfig(sce)
#Delete some unnecessary metadata
conf <- delMeta(conf, meta.to.del = c( "subsets_mt_sum", "subsets_mt_detected",
                       "outlier", "ratio_detected_sum", "outlier_ratio",
                       "discard", "sizeFactor", "total", "discard", "original_sample_name",
                       "cluster_k20", "cluster_k40", "cluster_k60", "cluster_k80" )
                )
# Change name of some metadata
conf <- modMetaName(conf,
                    meta.to.mod = c("sum", "detected", "subsets_mt_percent"),
                    new.name = c("umi counts", "detected genes", "% mt genes"))
# 

makeShinyApp(sce, conf, gene.mapping = TRUE, 
             shiny.title = "ShinyCell oligos app",
             default.dimred = c("TSNE1", "TSNE2"))



