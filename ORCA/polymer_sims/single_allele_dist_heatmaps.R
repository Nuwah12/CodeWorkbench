####################
# Script for generating single-allele distance heatmaps like in Boettiger et. al
# Noah Burget
# 4/18/24
####################
library(data.table)
library(pheatmap)
setwd('/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims')
N <- 900
PROMOTER <- 735
# Read in files
dir <- '3D_simulation/scenario19_v8/confs_txt/'
files <- list.files(dir)
xyz <- lapply(files, function(x) {
  read.table(paste0(dir,x), sep = " ", skip = 1)
})
# Get distance matrices for traces
#xyz <- sample(xyz,1000)
##########
#xyz <- Granta519_untr_filled_poly
#xyz <- sample(xyz,7000)
dist <- lapply(xyz, FUN=dist, upper=F, diag=T) # compute distance matrices for x,y,x coordinates
dist <- lapply(dist, FUN=as.matrix) # coerce to matrix
rm(xyz) # remove large object to save space
#s1.dist.median <- lapply(dist, median) # calculate median distance for all 
#s1.dist.median <- median(unlist(s1.dist.median)) # get median 
# Binarize distances
STRINGENT_THRESH <- 12 # set a stringent diatance threshold (units are arbitrary here)
s1.bin.dists <- do.call(rbind,lapply(dist, function(x){ifelse(x[,PROMOTER] <= STRINGENT_THRESH, 1, 0)})) # Binarize; 1 if probe in contact w/ promoter. else 0
# sort
ord <- cbind(s1.bin.dists, apply(s1.bin.dists, 1, function(x){which(x==1)[1]})) # find first occurence of '1' in table
s1.bin.dists.wt <- ord[order(ord[,N+1],decreasing=TRUE),] # sort according to ord
# make heatmap
pdf('060624_simulatedTraces_populationHeatmap_toPromoter_Perturbed.pdf')
pheatmap(s1.bin.dists.wt[,-(N+1)], cluster_rows=F, cluster_cols=F, show_colnames=FALSE, 
         scale="none", legend = FALSE, color = c('white','red'), border_color = NA)
dev.off()
##### Differential heatmap
pdf('060624_simulatedTraces_populationHeatmap_toPromoter_Perturbed-Minus-WT.pdf')
pheatmap(s1.bin.dists.diff[,-(N+1)], cluster_rows=F, cluster_cols=F, show_colnames=FALSE, 
         scale="none", legend = FALSE, color = c('blue','white','red'), border_color = NA)
dev.off()
