library(here)
library(scater)
library(clustree)
library(patchwork)

project <- "fire-mice"
source(here("src/colours.R"))

sce <- readRDS(here("processed",project,"sce_oligo_clusters_01.RDS"))

sce$cluster_oligo_k80 <- NULL
tree <- clustree(sce, prefix = "cluster_oligo_k", edge_arrow = FALSE)

pdf(here("outs",project,"plots", "cluster_tree_oligo.pdf"), height = 10, width = 7)
print(tree)
dev.off()

# test custom reverse tree, by adding number for ordering before and a 1 to keep 0s
sce$cluster_oligo_sorted_3.101 <- sce$cluster_oligo_k10
sce$cluster_oligo_sorted_0.1501 <- sce$cluster_oligo_k150
sce$cluster_oligo_sorted_2.1001 <- sce$cluster_oligo_k100
sce$cluster_oligo_sorted_1.1201 <- sce$cluster_oligo_k120


tree_reverse <- clustree(sce, prefix = "cluster_oligo_sorted_", edge_arrow = FALSE) +
  scale_color_discrete(name="clustering k", labels=c("150","120","100","10"))

# resolutions 
theme <-  theme( legend.position = "none") 
colour <- scale_colour_manual(values = cols) 

k100 <- plotTSNE(sce, colour_by = "cluster_oligo_k100", text_by="cluster_oligo_k100" ) + 
  theme + colour + ggtitle("Clustering k100")
k150 <- plotTSNE(sce, colour_by = "cluster_oligo_k150", text_by="cluster_oligo_k150" )+
  theme + colour + ggtitle("Clustering k150")
k120 <- plotTSNE(sce, colour_by = "cluster_oligo_k120", text_by="cluster_oligo_k120") +
  theme + colour + ggtitle("Clustering k120")
k10 <- plotTSNE(sce, colour_by = "cluster_oligo_k10", text_by="cluster_oligo_k10") +
  theme + colour + ggtitle("Clustering k10")


pdf(here("outs", project, "plots", "cluster_resolutions_and_tree.pdf"), height = 7, width = 15)
((k10+k100)/(k120+k150))  | tree_reverse
dev.off()

tiff(here("outs", project, "plots", "cluster_resolutions_and_tree.tiff"), height = 7, width = 15, res = 300, units = "in")
((k10+k100)/(k120+k150))  | tree_reverse
dev.off()
