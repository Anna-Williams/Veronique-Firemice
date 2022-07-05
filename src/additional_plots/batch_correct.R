library(here)
library(scater)
library(patchwork)

project <- "fire-mice"
sce <- readRDS(here("processed",project,"sce_oligo_clusters_01.RDS"))


# batch correction( not really umi_detected_mt, should move to another script)
batch <- plotReducedDim(sce, "TSNE_uncorrected", colour_by="chip") + ggtitle("Dimensional reduction before batch correction")
batch_correct <- plotReducedDim(sce, "TSNE", colour_by="chip") + ggtitle("Dimensional reduction after batch correction")

pdf(here("outs",project, "plots","batch_correction.pdf"), height = 7, width = 14)
print(batch + batch_correct)
dev.off()