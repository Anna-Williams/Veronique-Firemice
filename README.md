# Veronique-Firemice

This repository includes the source code used for the analysis of scRNAseq from  Reference: McNamara N. Microglia regulate central nervous system myelin growth and integrity Nature , (2022) doi: 10.1038/s41586-022-05534-y [Link](https://www.nature.com/articles/s41586-022-05534-y)

The analysed dataset is available as a ShinyApp: https://annawilliams.shinyapps.io/shinyApp_oligos_VM/

The raw data is available in GEO: GSE215440

The order the scripts in this repository were run is: 

- import.R
- QC_01.Rmd
- normalise_01.Rmd
- feature_selection_dimred_01.Rmd
- batch_correct.Rmd
- clustering_01.Rmd
- annottion_01.Rmd
- oligos.Rmd
- clustering_oligos_01.Rmd
- DA_WT_KO_oligos.Rmd
- DE_WT_KO_oligos*.Rmd
- compare_DE.R*
- generate_shiny_oligos.R
- additional_plots
