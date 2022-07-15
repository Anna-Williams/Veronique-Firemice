library(here)
library(scater)
library(patchwork)

project <- "fire-mice"


# from before QC
sce <- readRDS(here("processed", project, "sce_clusters_01.RDS"))
sce <- sce[, sce$cluster_k60 %in% c(5, 12, 15)]


umi <- plotColData(sce, x = "Sample", y = "sum") +
  ylab("UMI count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

detected <- plotColData(sce, x = "Sample", y = "detected") +
  ylab("Detected genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

mt <- plotColData(sce, x = "Sample", y = "subsets_mt_percent") +
  ylab("Mitochondrial percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(here("outs",project, "plots","umi_detected_mt_oligo_beforeQC.pdf"), height = 7, width = 14)
print(umi + detected + mt)
dev.off()

# TSNE before
detected <- plotTSNE(sce, colour_by = "detected", point_size = 0.2) +
  ggtitle("Detected Genes")
umi <- plotTSNE(sce, colour_by = "sum", point_size = 0.2) +
  ggtitle("UMI count")
mt <- plotTSNE(sce, colour_by = "subsets_mt_percent", point_size = 0.2) +
  ggtitle("Mitochondrial percentage")

# mark the cells filtered out
sce$filter <-  !( sce$sum >= 5000 &
  sce$detected >= 2000 &
  sce$subsets_mt_percent <= 10 &
  !(sce$cluster_k60 %in% c(15)))

filter <- plotTSNE(sce, colour_by = "filter") +
  ggtitle("Cells filtered out")

# save combination in pdf
pdf(here("outs",project, "plots","umi_detected_mt_oligo_TSNE.pdf"), height = 7, width = 8)
print((umi + mt)/( detected + filter))
dev.off()

# tiff
tiff(here("outs",project, "plots","umi_detected_mt_oligo_TSNE.tiff"), height = 7, width = 8, res=300, units = "in")
print((umi + mt)/( detected + filter))
dev.off()


# from after QC (same axis as before QC)
sce <- readRDS(here("processed",project,"sce_oligo_clusters_01.RDS"))

umi <- plotColData(sce, x = "Sample", y = "sum") +
  ylab("UMI count") + ylim(c(0,NA)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

detected <- plotColData(sce, x = "Sample", y = "detected") +
   ylab("Detected genes") + ylim(c(0,NA)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

mt <- plotColData(sce, x = "Sample", y = "subsets_mt_percent") +
   ylab("Mitochondrial percentage") + ylim(c(0,15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(here("outs",project, "plots","umi_detected_mt_oligo_afterQC_matchaxis.pdf"), height = 7, width = 14)
print(umi + detected + mt)
dev.off()

# with free axis, to be used if we don't show the before QC
umi <- plotColData(sce, x = "Sample", y = "sum") +
  ylab("UMI count") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

detected <- plotColData(sce, x = "Sample", y = "detected") +
  ylab("Detected genes") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

mt <- plotColData(sce, x = "Sample", y = "subsets_mt_percent") +
  ylab("Mithocondrial percentage") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(here("outs",project, "plots","umi_detected_mt_oligo_afterQC_freeaxis.pdf"), height = 7, width = 14)
print(umi + detected + mt)
dev.off()

