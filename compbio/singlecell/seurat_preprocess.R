##########
# Script for processing 10X scRNA-seq data with Seurat
# Noah Burget
##########

preprocess.Seurat <- function(work.dir,
                              sample.dir.prefix="SRR",
                              totalRNA.thresh.lo=500,
                              totalRNA.thresh.hi=25000,
                              featureRNA.thresh.lo=50,
                              featureRNA.thresh.hi=7500,
                              mitoPerc.thresh=25,
                              riboPerc.thresh=20){
  ### Gather data
  setwd(work.dir)
  data.dir <- paste0(work.dir,"/")
  aligned.dirs <- list.files(data.dir, pattern=sample.dir.prefix)
  
  ### Read in h5 and create Seurat Object
  sample.list <- setNames(lapply(aligned.dirs, function(x){
    print(paste0(data.dir,"/",x,"/outs/count/filtered_feature_bc_matrix.h5"))
    patient <- Read10X_h5(paste0(data.dir,"/",x,"/outs/filtered_feature_bc_matrix.h5"))
    seuratObj <- CreateSeuratObject(counts=patient, project=x)
    print(paste0("Cells in sample ",x," = ",ncol(seuratObj@assays$RNA@layers$counts)))
    return(seuratObj)
  }), aligned.dirs)
  
  ### Calculate QC metadata and plot
  seurat.qc <- list()
  for (name in names(sample.list)) {
    obj <- sample.list[[name]]
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
  
  ### Integrate samples
  features <- SelectIntegrationFeatures(object.list = normalized.list, nfeatures = 3000) # Select integration features
  normalized.list <- PrepSCTIntegration(object.list = normalized.list, anchor.features = features) # Prepare object for integration
  anchors <- FindIntegrationAnchors(object.list = normalized.list, normalization.method = "SCT", anchor.features = features) # Find integration anchors
  
  # Integrate data
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  integrated <- RunPCA(integrated)
  integrated <- RunUMAP(integrated, dims = 1:30)
  
  return(intgegrated)
}

