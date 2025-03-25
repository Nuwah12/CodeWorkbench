##########
# Script for processing 10X scRNA-seq data with Seurat
# Noah Burget
##########
library(Seurat); library(SeuratDisk)

dir <- ""
setwd(dir)

mcl.patients <- Read10X_h5(paste0(dir,""))
mcl.seurat <- CreateSeuratObject(counts=mcl.patients,
                                 project="",
                                 min.cells=3,
                                 min.features=200)
rm(mcl.patients); gc()

### QC 
mcl.seurat <- PercentageFeatureSet(mcl.seurat, "^MT-", col.name="percent.mito")
mcl.seurat <- PercentageFeatureSet(mcl.seurat, "^RP[SL]", col.name = "percent_ribo")
