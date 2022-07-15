#lipid/cholesterol: "Fdps", "Hmgcs1", "Abca5", "Sqle", "Scd2", "Idi1", "Hsd17b7", "Mvk", "Fdft1", "Cyp51", "Ncor2", "Me1",  "Clock", "Fads1", "Scd1", "Plpp3", "Abca2", "Asah2", "Gltp", "Gm2a", "Arsb", "Hexa", "Stard10", "Pcyt2", "Slc44a1", "Gpd1", "Chpt1", "Agpat4", "Mtmr2", "Plekha1", "Enpp6", "Mgll", "Mid1ip1", "Plin3", "Ptgds", "Txnrd1", "Hsd17b11", "Ppp1cc", "Hexa", "Fabp5", "Gpx4", "Chpt1", "Msmo1", "Elovl5", "Cyp2j6", "Srd5a1"
# TGFB1: "adi1", "aif1l", "aldh2", "anxa2", "apoe", "arc", "bcl3", "bin1", "bmp1", "cadm1", "camk2n1", "ccdc85b", "cdkn1a", "csf1", "cst3", "ctsb", "daam1", "dnaja1", "dock9", "dynll1", "eef1a1", "fabp5", "fgfr2", "flnb", "fnbp1", "fos", "foxj1", "fscn1", "fth1", "gadd45a", "gatm", "gprc5b", "gsn", "hexa", "hsd17b10", "hspa5", "htra1", "iars1", "idi1", "ifrd1", "il33", "jun", "junb", "klf9", "lamp2", "ldha", "ly6a", "malat1", "map1b", "masp1", "mdm2", "msmo1", "msn", "ndst1", "nme2", "pak3", "parp3", "pcolce2", "pdgfa", "phgdh", "pink1", "ppp2r2a", "psat1", "ptgds", "rasgrp3", "rasl11b", "s100a6", "sbno2", "scd", "scg5", "serpina3", "serpinb1", "slc1a2", "slc39a1", "slc3a2", "slc4a2", "socs3", "sparc", "sparcl1", "stat3", "timp2", "tnfaip6", "trim2", "tsc22d1", "tspan7", "tuba1a", "tubb3", "txnrd1", "znf365"

library(scater)
library(here)
project <- "fire-mice"
sce <- readRDS(here("processed", project, "sce_oligo_clusters_01.RDS"))
# I think I removed all duplicates already, I keep the unique to be sure
lipid <- unique(c("Fdps", "Hmgcs1", "Abca5", "Sqle", "Scd2", "Idi1", "Hsd17b7", "Mvk", "Fdft1", "Cyp51", "Ncor2", "Me1",  "Clock", "Fads1", "Scd1", "Plpp3", "Abca2", "Asah2", "Gltp", "Gm2a", "Arsb",  "Stard10", "Pcyt2", "Slc44a1", "Gpd1", "Agpat4", "Mtmr2", "Plekha1", "Enpp6", "Mgll", "Mid1ip1", "Plin3", "Ptgds", "Txnrd1", "Hsd17b11", "Ppp1cc", "Hexa", "Fabp5", "Gpx4", "Chpt1", "Msmo1", "Elovl5", "Cyp2j6", "Srd5a1"))
# the notation should be with first letter capital
tgfb1 <- stringr::str_to_sentence(c("adi1", "aif1l", "aldh2", "anxa2", "apoe", "arc", "bcl3", "bin1", "bmp1", "cadm1", "camk2n1", "ccdc85b", "cdkn1a", "csf1", "cst3", "ctsb", "daam1", "dnaja1", "dock9", "dynll1", "eef1a1", "fabp5", "fgfr2", "flnb", "fnbp1", "fos", "foxj1", "fscn1", "fth1", "gadd45a", "gatm", "gprc5b", "gsn", "hexa", "hsd17b10", "hspa5", "htra1", "iars1", "idi1", "ifrd1", "il33", "jun", "junb", "klf9", "lamp2", "ldha", "ly6a", "malat1", "map1b", "masp1", "mdm2", "msmo1", "msn", "ndst1", "nme2", "pak3", "parp3", "pcolce2", "pdgfa", "phgdh", "pink1", "ppp2r2a", "psat1", "ptgds", "rasgrp3", "rasl11b", "s100a6", "sbno2", "scd", "scg5", "serpina3", "serpinb1", "slc1a2", "slc39a1", "slc3a2", "slc4a2", "socs3", "sparc", "sparcl1", "stat3", "timp2", "tnfaip6", "trim2", "tsc22d1", "tspan7", "tuba1a", "tubb3", "txnrd1", "znf365"))
# some genes were not found in our dataset, subset to only keep if they are there
tgfb1_sub <- tgfb1[tgfb1 %in% row.names(sce)]


# heatmaps
pdf(here("outs",project, "plots","lipid_genes.pdf"), height = 14, width = 7)
plotGroupedHeatmap(sce, features=lipid, group="cluster_oligo", 
                   center=TRUE, scale=TRUE, cluster_cols=FALSE)#, zlim=c(-3, 3))
dev.off()

pdf(here("outs",project, "plots","tgfb1_genes.pdf"), height = 14, width = 7)
plotGroupedHeatmap(sce, features=tgfb1_sub, group="cluster_oligo", 
                   center=TRUE, scale=TRUE, cluster_cols=FALSE)#, zlim=c(-3, 3))
dev.off()

# tiff
tiff(here("outs",project, "plots","lipid_genes.tiff"), height = 14, width = 7, res=300, units = "in")
plotGroupedHeatmap(sce, features=lipid, group="cluster_oligo", 
                   center=TRUE, scale=TRUE, cluster_cols=FALSE)#, zlim=c(-3, 3))
dev.off()

tiff(here("outs",project, "plots","tgfb1_genes.tiff"), height = 14, width = 7, res=300, units = "in")
plotGroupedHeatmap(sce, features=tgfb1_sub, group="cluster_oligo", 
                   center=TRUE, scale=TRUE, cluster_cols=FALSE)#, zlim=c(-3, 3))
dev.off()

# I tested with violins, but it is not helpful
pdf(here("outs",project, "plots","lipid_genes_violin.pdf"), height = 14, width = 7)
plotExpression(sce, features=lipid, 
               x="genotype", colour_by="genotype")
dev.off()

