###############
# Script for clustering the simulation results
# 5/31/24
# Noah Burget
##############
library(abind)
library(pheatmap)

setwd('/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims')

dir.wt <- '3D_simulation/scenario19_v8/confs_txt/'
dir.per <- '3D_simulation/scenario21_v8/confs_txt/'
wt.full <- lapply(list.files(dir.wt), function(x) {
  read.table(paste0(dir.wt,x), sep = " ", skip = 1)
})
wt <- lapply(wt.full,dist,upper=T,diag=T)
wt <- lapply(wt,as.matrix)
wt.distToProm <- lapply(wt, function(x){return(x[735,])})
#per <- lapply(list.files(dir.per), function(x) {
  read.table(paste0(dir.per,x), sep = " ", skip = 1)
})

########## Clustering WT
wt.cluster <- do.call(rbind, wt.distToProm)
wt.full <- lapply(wt.full,dist)
wt.kmeans <- kmeans(wt.cluster, centers = 10, iter.max = 20)
wt.clusters <- list()
for (i in 1:10){
  wt.clusters[[i]] <- wt.full[which(wt.kmeans$cluster==i)]
}
rm(wt, wt.full)
wt.clusters <- lapply(wt.clusters, function(x){lapply(x,as.matrix)})
map <- abind(wt.clusters[[4]], along=3)
map <- apply(map.c1, c(1,2), mean)
pheatmap(map, cluster_cols = F, cluster_rows = F, border_color = NA, labels_row = NA, labels_col = NA)

