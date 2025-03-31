##########
# Script for processing 10X scRNA-seq data with Seurat
# Noah Burget
##########
library(Seurat); library(SeuratDisk); library(DoubletFinder)

### Gather data
work.dir <- ""
setwd(work.dir)
data.dir <- paste0(work.dir,"/")
aligned.dirs <- list.files(data.dir, pattern="SRR")

### Read in h5 and create Seurat Object
mcl.patients.seurat <- setNames(lapply(aligned.dirs, function(x){
  print(paste0(data.dir,"/",x,"/outs/count/filtered_feature_bc_matrix.h5"))
  patient <- Read10X_h5(paste0(data.dir,"/",x,"/outs/filtered_feature_bc_matrix.h5"))
  seuratObj <- CreateSeuratObject(counts=patient, project=x)
  print(paste0("Cells in sample ",x," = ",ncol(seuratObj@assays$RNA@layers$counts)))
  return(seuratObj)
}), aligned.dirs)

### Calculate QC metadata and plot
seurat.qc <- list()
for (name in names(mcl.patients.seurat)) {
  obj <- mcl.patients.seurat[[name]]
  print(name)
  
  obj <- PercentageFeatureSet(obj, "^MT-", col.name = "percent.mito")
  obj <- PercentageFeatureSet(obj, "^RP[SL]", col.name = "percent.ribo")
  
  pdf(paste0(work.dir, "/qc_figs/", name, "_qc_violin.pdf"))
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4)) # Wrapping in print() solved the corrupt PDF issue
  dev.off()
  
  seurat.qc[[name]] <- obj
}

### Quality control (QC)
apply.qc.filters <- function(sobj, totalRNA.thresh.lo, totalRNA.thresh.hi, featureRNA.thresh.lo, featureRNA.thresh.hi, mitoPerc.thresh, riboPerc.thresh){
  return(subset(sobj, nCount_RNA > totalRNA.thresh.lo & nCount_RNA < totalRNA.thresh.hi 
                & nFeature_RNA > featureRNA.thresh.lo & nFeature_RNA < featureRNA.thresh.hi
                & percent.mito < mitoPerc.thresh & percent.ribo < riboPerc.thresh))
}

# Apply QC to all datasets
#all.qc.thresh <- data.frame("SampleName"=character(), "totalRNA.low"=integer(), "totalRNA.hi"=integer(), "featRNA.low"=integer(), "featRNA.hi"=integer(), "percentMito.hi"=integer(), "riboPercent.hi"=integer())

# Normalize and Scale the data for each sample separately
SCTransform.norm <- function(obj){
  obj <- SCTransform(obj, vars.to.regress=c("percent.mito","percent.ribo"))
  return(obj)
}

### Apply same QC to all samples (simple)
seurat.post.qc <- lapply(seurat.qc, apply.qc.filters, 500, 25000, 50, 7500, 25, 20)

### Apply SCTransform to each sample separately
normalized.list <- lapply(seurat.post.qc, SCTransform.norm) # SCtransform normalize

### For comparison, do an unintegrated run 
s1 <- normalized.list$SRR15871909
no.integrate <- merge(x=s1, y=as.vector(normalized.list[-1])) # Just merge the objects, no integration
no.integrate <- SelectIntegrationFeatures(object.list = normalized.list, nfeatures = 3000) # PCA needs variable features, and we can't use FindVariableFeatures because the merged object contains multiple SCTransform models, one for each sample
no.integrate <- RunPCA(no.integrate) # Run PCA
no.integrate <- RunUMAP(no.integrate, dims = 1:30) # Run UMAP

DimPlot(no.integrate) # Plot unintegrated 

### Integrate samples
features <- SelectIntegrationFeatures(object.list = normalized.list, nfeatures = 3000) # Select integration features
normalized.list <- PrepSCTIntegration(object.list = normalized.list, anchor.features = features) # Prepare object for integration
anchors <- FindIntegrationAnchors(object.list = normalized.list, normalization.method = "SCT", anchor.features = features) # Find integration anchors

# Integrate data
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)

DimPlot(integrated)


