###############
# Steps taken to convert .rds (Seurat object) to Anndata h5 (.h5ad)
# Noah Burget
###############
library(Seurat)
library(SeuratDisk)

# Read in Seurat object saved as rds
PLN_RNAAssay <- readRDS("/mnt/alvand/VahediLab/Projects/Multiome_HPAP/Shared/Spleen_RNAAssay_SCT.rds")

#Write as Seurat h5
SaveH5Seurat(PLN_RNAAssay, "/mnt/alvand/noah/Spleen_RNAAssay_SCT.h5Seurat")

# Convert Seurat h5 to Anndata h5
Convert("/mnt/alvand/noah/Spleen_RNAAssay_SCT.h5Seurat", dest="h5ad")
