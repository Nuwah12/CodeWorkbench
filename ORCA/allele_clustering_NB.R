library(factoextra)
library(abind)
library(ggplot2)
library(NMF)
library(umap)
library(data.table)
library(pheatmap)
library(clusterSim)
library(fpc)
library(cluster)
library(SMFilter)
library(sClust)
library(plotly)
library(stackoverflow)
library(MCL)
library(h2o)
library(dbscan)
library(ggpubr)
library(gridExtra)
###################################################################################################################
# Allele Clustering: Indentifying Populations of Alleles via Distance Matrices from ORCA experiments.
# 09-15-23
# Noah Burget
###################################################################################################################
setwd("/mnt/data0/noah/analysis/ORCA_analysis/allele_clustering")
source("processChrTracer.R")
source("plotContactFrequency.R")
###################################################################################################################
str.untreated <- '231020_Granta519cl27_untreated_4phBl_30minHyb_30step_allFits'
str.treated <- '230902_Granta519cl27_24hdTAG_MYC5p_30mHyb_4phBl_30step_allfits'
raw.untreated <- paste0("ChTracer3_output/biological_conditions/",str.untreated,".csv")
raw.treated <- paste0("ChTracer3_output/biological_conditions/",str.treated,".csv")
n <- 30
Forward_all_distances_adj.untreated <- process.chrtracer(raw.untreated, n, T)
Forward_all_distances_adj.treated <- process.chrtracer(raw.treated, n, T)
full_dist.untreated <- make.distTable(Forward_all_distances_adj.untreated, useDiag = FALSE, fullSize = (n*(n-1))/2, probe.to.fix=19)
full_dist.treated <- make.distTable(Forward_all_distances_adj.treated, useDiag = FALSE, fullSize = (n*(n-1))/2, probe.to.fix=19)
dist.upper.lens.untreated <- unlist(lapply(Forward_all_distances_adj.untreated, function(x) length(which(!is.na(x)))))
dist.upper.lens.treated <- unlist(lapply(Forward_all_distances_adj.treated, function(x) length(which(!is.na(x)))))
###################################################################################################################
# Form of input data:
# TRACES/ROWS X DISTANCES/COLS

##### 0: KMeans on distance table #####
fviz_nbclust(full_dist, kmeans, "wss")
kmean.res <- kmeans(full_dist.noDTag, 3)
full_dist.noDTag$kmeans.cluster <- kmean.res$cluster


##### 1.1: Nonnegative Matrix Factorization (NMF) #####
set.seed(123)
# Estimate rank for input
nmf.res <- nmf(x = as.matrix(full_dist.treated), 
               rank = 2, 
               method = 'lee',
               seed = 'random')
# Extract Basis matrix 
nmf.basis <- as.data.frame(basis(nmf.res))
# Plot basis matrix
ggplot(data=nmf.basis, aes(x = V1, y = V2))+
  geom_point()
#pdf(paste0(str,'snmf-l_basisHeatmap.pdf'))
#basismap(nmf.res[["snmf/l"]])
#dev.off()

##### 1.1.a: PCA #####
set.seed(123)
pca.obj <- prcomp(full_dist.all, center=F, scale=T)
pca.res <- as.data.frame(pca.obj$x[,1:2])
# Plot first 2 PCs
ggplot(data = pca.res, aes(x = PC1, y = PC2))+
  geom_point()

##### 1.2: KMeans clustering of NMF basis matrix #####
set.seed(123)
# Estimate optimal k
fviz_nbclust(nmf.basis, kmeans, method='silhouette')
# KMeans clustering
kmeans.res <- kmeans(nmf.basis,3)
nmf.basis$kmeans.cluster <- kmeans.res$cluster
# Plot basis matrix by cluster
ggplot(data = nmf.basis, aes(x = V1, y = V2))+
  geom_point(aes(color = factor(kmeans.cluster)))
# Use factoextra package for better cluster viz
pdf(paste0(str.untreated,'_NMF_clusterPlot_30steps.pdf'))
fviz_cluster(kmeans.res,
             data = nmf.basis[,-3],
             nmf.basispalette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw())
dev.off()
##### 1.2.a: Assess clustering quality
# Davies-Bouldin's cluster separation measure / index
index.DB(full_dist.DTag, kmeans.res$cluster)
# Calinski-Harabasz index
calinhara(full_dist.DTag, kmeans.res$cluster, cn=3)

##### 1.3: KMeans clustering of first 2 PCs #####
set.seed(123)
# Estimate optimal k
#fviz_nbclust(pca.res, kmeans, method='wss')
kmeans.res <- kmeans(as.matrix(pca.res), 3)
pdf(paste0(str,'_PCA_clusterPlot.pdf'))
fviz_cluster(kmeans.res,
             data = pca.res[,-3],
             nmf.basispalette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw())
dev.off()
pca.res$kmeans.cluster <- kmeans.res$cluster

##### 2: UMAP #####
set.seed(123)
umap.res <- as.data.frame(umap(full_dist)$layout)
# Plot UMAP embedding
ggplot(data = umap.res, aes(x = V1, y = V2))+
  geom_point()
##### 2.1: KMeans clustering of UMAP embedding #####
# Estimate optimal k
#fviz_nbclust(umap.res, kmeans, method='silhouette')
# Run
kmeans.res <- kmeans(umap.res, 3)
pdf(paste0(str,'_UMAP_clusterPlot_22step.pdf'))
fviz_cluster(kmeans.res,
             data = umap.res[,-3],
             nmf.basispalette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw())
dev.off()
umap.res$kmeans.cluster <- kmeans.res$cluster

#### Convert distances in each cluster to contact frequency #####
group1.bind <- plotContactFrequency(Forward_all_distances_adj.treated[which(dist.upper.lens.treated==n^2)[which(nmf.basis.treated$kmeans.cluster == 1)]],n)
group2.bind <- plotContactFrequency(Forward_all_distances_adj.treated[which(dist.upper.lens.treated==n^2)[which(nmf.basis.treated$kmeans.cluster == 2)]],n)
group3.bind <- plotContactFrequency(Forward_all_distances_adj.treated[which(dist.upper.lens.treated==n^2)[which(nmf.basis.treated$kmeans.cluster == 3)]],n)
breaks <- c(seq(0.35, 0.85, length.out=ceiling(1000)))

##### Fix probe 19
group1.bind[,19] <- rowMeans(cbind(group1.bind[,18], group1.bind[,20]))
group1.bind[19,] <- colMeans(rbind(group1.bind[18,], group1.bind[20,]))
group1.bind[19,19] <- 1
group2.bind[,19] <- rowMeans(cbind(group2.bind[,18], group2.bind[,20]))
group2.bind[19,] <- colMeans(rbind(group2.bind[18,], group2.bind[20,]))
group2.bind[19,19] <- 1
group3.bind[,19] <- rowMeans(cbind(group3.bind[,18], group3.bind[,20]))
group3.bind[19,] <- colMeans(rbind(group3.bind[18,], group3.bind[20,]))
group3.bind[19,19] <- 1

pdf(paste0(str.treated,'_NMF_cluster1_',n,'step.pdf'))
pheatmap(group1.bind, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#FFFFFF","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000),breaks=breaks,border_color=NA)
dev.off()
pdf(paste0(str.treated,'_NMF_cluster2_',n,'step.pdf'))
pheatmap(group2.bind, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#FFFFFF","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000),breaks=breaks,border_color=NA)
dev.off()
#breaks <- c(seq(0.4, 0.7, length.out=ceiling(500)))
pdf(paste0(str.treated,'_NMF_cluster3_',n,'step.pdf'))
pheatmap(group3.bind, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#FFFFFF","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000),breaks=breaks,border_color=NA)
dev.off()

########## OLIVE input file ##########
all_xyz.untreated <- process.chrtracer(raw.untreated, 30, T, coords=TRUE)
all_xyz.treated <- process.chrtracer(raw.treated, 30, T, coords=TRUE)
#all_xyz <- lapply(seq_along(all_xyz), function(i){
#  new <- rep(i, nrow(all_xyz[[i]]))
#  cbind(all_xyz[[i]],s=new)
#})
#all_xyz <- lapply(seq_along(all_xyz), function(i){
#  new <- 1:22
#  cbind(all_xyz[[i]],readout=new)
#})
#all_xyz <- lapply(seq_along(all_xyz), function(i){
#  cbind(all_xyz[[i]],fov=1)
#})
#all_xyz <- do.call(rbind, all_xyz)

########## HClust with RMSD matrix ##########
rmsd.dist <- as.dist(rmsd.mat)
pam.res <- pam(rmsd.dist, 3)
plot(silhouette(pam.res))
full_dist <- as.data.frame(full_dist)
full_dist$kmeans.cluster <- pam.res$clustering

########## Autoencoder ##########
h2o.init(max_mem_size = "50g")
features <- as.h2o(as.data.frame(full_dist)[,1:ncol(full_dist)-1])
hyper_grid <- list(hidden = list(
  c(50),
  c(100), 
  c(300, 100, 300),
  c(100, 50, 100),
  c(250, 100, 50, 100, 250)
))
ae_grid <- h2o.grid(
  algorithm = 'deeplearning',
  x = seq_along(features),
  training_frame = features,
  grid_id = 'autoencoder_grid',
  autoencoder = TRUE,
  activation = 'Tanh',
  hyper_params = hyper_grid,
  sparse = TRUE,
  ignore_const_cols = FALSE,
  seed = 123
)
h2o.getGrid('autoencoder_grid', sort_by = 'mse', decreasing = FALSE)
ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = 5,
  activation = 'Tanh',
  sparse = TRUE
)
encodings <- as.data.frame(h2o.deepfeatures(ae1, features, layer = 1))
ggplot(encodings, aes(x=DF.L1.C1, y=DF.L1.C2))+
  geom_point()
# KMeans
fviz_nbclust(encodings, kmeans, method='silhouette')

trace <- as.data.frame(all_xyz[[1000]])
trace$cluster <- kmeans(trace, centers = 2)[["cluster"]]
trace[23,] <- c(mean(trace$x),mean(trace$y),mean(trace$z),-1)
p <- plotly::plot_ly(x=trace$x,y=trace$y,z=trace$z, mode="markers", type='scatter3d', color = trace$cluster) 
plotly::add_text(p,x=trace$x,y=trace$y,z=trace$z,text=1:nrow(trace))

########## Calculate distances to center ##########
# Get center of each matrix
# Calculate distance between each probe and the trace's center
all_xyz <- lapply(all_xyz.untreated,as.data.frame)
dist_to_center_per_probe <- function(matrix){
  center <- colMeans(matrix)
  d <- apply(matrix,1,function(x){
    distance <- sqrt((center[1] - x[1])^2+(center[2] - x[2])^2+(center[3] - x[3])^2)
  })
  matrix['dist_to_center'] <- d
  return(matrix)
}
# Make table of traces x distances to center (2065x22)
avg_dist_per_probe <- lapply(all_xyz, dist_to_center_per_probe) # Apply func to xyz tables 
avg_dist_per_probe <- lapply(avg_dist_per_probe, function(x){return(x[,'dist_to_center'])})
###### !!! CORRECT PROBE 19 BY AVERAGING ADJACENT PROBES
TO_FIX <- 19
avg_dist_per_probe <- lapply(avg_dist_per_probe, function(x){
  x[TO_FIX] <- mean(c(x[TO_FIX-1],x[TO_FIX+1]))
  return(x)
})
# CLUSTER/CONFIG 1 = LOWER-HALF INTER
# CLUSTER/CONFIG 2 = UPPER-HALF INTER
# CLUSTER/CONFIG 3 = MIDDLE-CONF
#### Untreated = 3,1,2
#### Treated 2,1,3
avg_dist_per_probe <- do.call(rbind, avg_dist_per_probe)
### Plotting without cluster
#avg_dist_per_probe.1 <- colMedians(avg_dist_per_probe[which(nmf.basis$kmeans.cluster==1),])
#avg_dist_per_probe.2 <- colMedians(avg_dist_per_probe[which(nmf.basis$kmeans.cluster==2),])
#avg_dist_per_probe.3 <- colMedians(avg_dist_per_probe[which(nmf.basis$kmeans.cluster==3),])
#### Order clusters here!
#avg_dists_from_center <- melt(data.frame(clus1=avg_dist_per_probe.2,clus2=avg_dist_per_probe.1,clus3=avg_dist_per_probe.3))
avg_dists_from_center.noCluster <- colMedians(avg_dist_per_probe)
#avg_dists_from_center$x <- 1:30
# Line plot for median distances
ggplot(data=data.frame(x=avg_dists_from_center.noCluster), aes(x=1:30, y=x))+
  geom_line(size=1.5)+
  ylab("Distance to center (nm)")+
  xlab("Probe")+
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"), labels=c(1,2,3))+
  labs(color="Trace Configuration")+
  theme(legend.position=c(0.9,0.9))+
  scale_x_continuous(breaks=1:30, labels=1:30)+
  theme(panel.grid.minor.x = element_blank())+
  ylim(c(200, 400))
ggsave(paste0(str.untreated,"_medianDistToCenter_withCorrectedProbe19_noCluster.pdf"))

########## "Tightness' of TAD probes ##########
# tightness = avg. distance from probe to geometric center of TAD
tad1 <- c(1:5)
tad2 <- c(8:18)
tad3 <- c(19:22)
all_xyz_tad1 <- unlist(lapply(lapply(lapply(lapply(all_xyz,function(x){return(x[tad1,])}), dist_to_center_per_probe),function(x){return(x[,'dist_to_center'])}),mean))
all_xyz_tad2 <- unlist(lapply(lapply(lapply(lapply(all_xyz,function(x){return(x[tad2,])}), dist_to_center_per_probe),function(x){return(x[,'dist_to_center'])}),mean))
all_xyz_tad3 <- unlist(lapply(lapply(lapply(lapply(all_xyz,function(x){return(x[tad3,])}), dist_to_center_per_probe),function(x){return(x[,'dist_to_center'])}),mean))

##### TAD 1 tightness by cluster #####
all_xyz_tad1.c1 <- all_xyz_tad1[which(nmf.basis$kmeans.cluster==1)]
all_xyz_tad1.c2 <- all_xyz_tad1[which(nmf.basis$kmeans.cluster==2)]
all_xyz_tad1.c3 <- all_xyz_tad1[which(nmf.basis$kmeans.cluster==3)]

c1 <- ggplot(data=data.frame(x=all_xyz_tad1.c1), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
c2 <- ggplot(data=data.frame(x=all_xyz_tad1.c2), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
c3 <- ggplot(data=data.frame(x=all_xyz_tad1.c3), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
grid.arrange(c1,c2,c3,nrow=2)
pdf(paste0(str,'_TAD1_distsToGeomCenter.pdf'))
grid.arrange(c1,c2,c3,nrow=2)
dev.off()
# plot CDF
all_xyz_tad1.c1 <- data.frame(x=all_xyz_tad1.c1,y="c1")
all_xyz_tad1.c2 <- data.frame(x=all_xyz_tad1.c2,y="c2")
all_xyz_tad1.c3 <- data.frame(x=all_xyz_tad1.c3,y="c3")
tad1.allclus <- rbind(all_xyz_tad1.c1,all_xyz_tad1.c2,all_xyz_tad1.c3)
ggplot(as.data.frame(tad1.allclus), aes(x=x,col=y))+
  stat_ecdf()+
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"), labels=c("1","2","3"),)+
  labs(color="Trace Configuration")+
  theme(legend.position=c(0.85,0.75))+
  xlab("Distance from Geometric Center (nm)")+
  ylab("Proportion")
ggsave(paste0(str,"_TAD1_allClusters_distFromTADcenter.pdf"))

##### TAD 2 tightness by cluster #####
all_xyz_tad2.c1 <- all_xyz_tad2[which(nmf.basis$kmeans.cluster==1)]
all_xyz_tad2.c2 <- all_xyz_tad2[which(nmf.basis$kmeans.cluster==2)]
all_xyz_tad2.c3 <- all_xyz_tad2[which(nmf.basis$kmeans.cluster==3)]
c1 <- ggplot(data=data.frame(x=all_xyz_tad2.c1), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
c2 <- ggplot(data=data.frame(x=all_xyz_tad2.c2), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
c3 <- ggplot(data=data.frame(x=all_xyz_tad2.c3), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
grid.arrange(c1,c2,c3,nrow=2)
pdf(paste0(str,'_TAD2_distsToGeomCenter.pdf'))
grid.arrange(c1,c2,c3,nrow=2)
dev.off()
# plot CDF
all_xyz_tad2.c1 <- data.frame(x=all_xyz_tad2.c1,y="c1")
all_xyz_tad2.c2 <- data.frame(x=all_xyz_tad2.c2,y="c2")
all_xyz_tad2.c3 <- data.frame(x=all_xyz_tad2.c3,y="c3")
tad2.allclus <- rbind(all_xyz_tad2.c1,all_xyz_tad2.c2,all_xyz_tad2.c3)
ggplot(as.data.frame(tad2.allclus), aes(x=x,col=y))+
  stat_ecdf()+
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"), labels=c("1","2","3"),)+
  labs(color="Trace Configuration")+
  theme(legend.position=c(0.85,0.75))+
  xlab("Distance from Geometric Center (nm)")+
  ylab("Proportion")
ggsave(paste0(str,"_TAD2_allClusters_distFromTADcenter.pdf"))

##### TAD 3 tightness by cluster #####
all_xyz_tad3.c1 <- all_xyz_tad3[which(nmf.basis$kmeans.cluster==1)]
all_xyz_tad3.c2 <- all_xyz_tad3[which(nmf.basis$kmeans.cluster==2)]
all_xyz_tad3.c3 <- all_xyz_tad3[which(nmf.basis$kmeans.cluster==3)]
c1 <- ggplot(data=data.frame(x=all_xyz_tad3.c1), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
c2 <- ggplot(data=data.frame(x=all_xyz_tad3.c2), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
c3 <- ggplot(data=data.frame(x=all_xyz_tad3.c3), aes(x=x))+
  geom_histogram()+
  xlim(c(0,800))+
  xlab("Distance (nm)")
pdf(paste0(str,'_TAD3_distsToGeomCenter.pdf'))
grid.arrange(c1,c2,c3,nrow=2)
dev.off()
# plot CDF
all_xyz_tad3.c1 <- data.frame(x=all_xyz_tad3.c1,y="c1")
all_xyz_tad3.c2 <- data.frame(x=all_xyz_tad3.c2,y="c2")
all_xyz_tad3.c3 <- data.frame(x=all_xyz_tad3.c3,y="c3")
tad3.allclus <- rbind(all_xyz_tad3.c1,all_xyz_tad3.c2,all_xyz_tad3.c3)
ggplot(as.data.frame(tad3.allclus), aes(x=x,col=y))+
  stat_ecdf()+
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"), labels=c("1","2","3"),)+
  labs(color="Trace Configuration")+
  theme(legend.position=c(0.85,0.75))+
  xlab("Distance from Geometric Center (nm)")+
  ylab("Proportion")
ggsave(paste0(str,"_TAD3_allClusters_distFromTADcenter.pdf"))

######### Modelling loop extrusion ##########
promoter <- 25
ebf1.bind <- c(8,11)
# Extract all vectors of distances from the promoter
promoter.dists <- lapply(Forward_all_distances_adj[which(dist.upper.lens==n^2)], function(x){return(x[8:17,promoter])})
promoter.closest.probe <- unlist(lapply(promoter.dists, which.min))+ebf1.bind[1]-1
ggplot(data=data.frame(x=promoter.closest.probe), aes(x=factor(x)))+geom_bar()+
  xlab("Probe")

######### Measuring mean distances from E1, E2, and Promoter to all other probes
all_xyz <- lapply(all_xyz, as.data.frame)
# E1
e1 <- 8 #E1 probe number
all_xyz <- lapply(all_xyz, function(matrix){
  e1 <- (as.matrix(dist(matrix))[,e1])
  return(cbind(matrix,e1))
})
# E2
e2 <- 11 #E2 probe number
all_xyz.untreated <- lapply(all_xyz.untreated, function(matrix){
  e2 <- (as.matrix(dist(matrix))[,e2])
  return(cbind(matrix,e2))
})
all_xyz.treated <- lapply(all_xyz.treated, function(matrix){
  e2 <- (as.matrix(dist(matrix))[,e2])
  return(cbind(matrix,e2))
})
# Promoter
prom <- 25 #Promoter probe number
all_xyz.untreated <- lapply(all_xyz.untreated, function(matrix){
  prom <- (as.matrix(dist(matrix))[,prom])
  return(cbind(matrix,prom))
})
all_xyz.treated <- lapply(all_xyz.treated, function(matrix){
  prom <- (as.matrix(dist(matrix))[,prom])
  return(cbind(matrix,prom))
})
##### Weird promoter distance values??????????
e1.to.prom.untreated <- unlist(lapply(all_xyz.untreated, function(x){return(x[8,4])}))
e1.to.prom.treated <- unlist(lapply(all_xyz.treated, function(x){return(x[8,4])}))
e2.to.prom.untreated <- unlist(lapply(all_xyz.untreated, function(x){return(x[11,4])}))
e2.to.prom.treated <- unlist(lapply(all_xyz.treated, function(x){return(x[11,4])}))
##### Subset by cluster
all_xyz.c1 <- lapply(all_xyz[which(nmf.basis$kmeans.cluster==1)],function(x){return(x[,c('e1','e2','prom')])})
all_xyz.c2 <- lapply(all_xyz[which(nmf.basis$kmeans.cluster==2)],function(x){return(x[,c('e1','e2','prom')])})
all_xyz.c3 <- lapply(all_xyz[which(nmf.basis$kmeans.cluster==3)],function(x){return(x[,c('e1','e2','prom')])})
##### Get mean values for each column in each cluster
# NOTE: HERE WE SWAP THE CLUSTER ASSIGNMENTS AROUND FOR CONTINUITY'S SAKE
# CLUSTER/CONFIG 1 = LOWER-HALF INTER
# CLUSTER/CONFIG 2 = UPPER-HALF INTER
# CLUSTER/CONFIG 3 = MIDDLE-CONF
#### Untreated = 3,1,2
#### Treated 2,1,3
element <- 'prom'
c1 <- colMeans(do.call(rbind,lapply(all_xyz.c2, function(x){return(x[,element])})))
c2 <- colMeans(do.call(rbind,lapply(all_xyz.c1, function(x){return(x[,element])})))
c3 <- colMeans(do.call(rbind,lapply(all_xyz.c3, function(x){return(x[,element])})))
dists <- melt(data.frame(cluster1=c1, cluster2=c2, cluster3=c3))
dists$probe <- 1:30
##### Plot
ggplot(data=dists, aes(x=probe, y=value, group=variable, color=variable))+
  geom_line(size=1.5)+
  ylab("Distance to MYC Promoter (Probe 25) (nm)")+
  xlab("Probe")+
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"), labels=c(1,2,3))+
  labs(color="Trace Configuration")+
  theme(legend.position=c(0.2,0.2))+
  ylim(0,1000)+
  scale_y_continuous(breaks=seq(0,1000, by=100))
ggsave(paste0(str,"_medianDistToMYCPromoter-Probe23.pdf"))


########## By-trace interaction map #########
prom <- 25 #Promoter probe number
all_xyz.untreated <- lapply(all_xyz.untreated, function(matrix){
  prom <- (as.matrix(dist(matrix))[,prom])
  return(cbind(matrix,prom))
})
all_xyz.treated <- lapply(all_xyz.treated, function(matrix){
  prom <- (as.matrix(dist(matrix))[,prom])
  return(cbind(matrix,prom))
})
###### !!! CORRECT PROBE 19 BY AVERAGING ADJACENT PROBES
TO_FIX <- 19
all_xyz.untreated <- lapply(all_xyz.untreated, function(x){
  x[TO_FIX,4] <- mean(c(x[TO_FIX-1,4],x[TO_FIX+1,4]))
  return(x)
})
all_xyz.treated <- lapply(all_xyz.treated, function(x){
  x[TO_FIX] <- mean(c(x[TO_FIX-1],x[TO_FIX+1]))
  return(x)
})
# Determine distance cutoffs
source("getContactFreqDistanceCutoff.R")
untreated.cutoff <- getContactFreqDistanceCutoff(Forward_all_distances_adj.untreated[which(dist.upper.lens.untreated==n^2)])
treated.cutoff <- getContactFreqDistanceCutoff(Forward_all_distances_adj.treated[which(dist.upper.lens.treated==n^2)])
# Binarize distances to promoter
all_xyz.untreated.promContact <- lapply(lapply(all_xyz.untreated, function(x){ifelse(x[,'prom'] < untreated.cutoff, 1, 0)}),as.vector)
all_xyz.treated.promContact <- lapply(lapply(all_xyz.treated, function(x){ifelse(x[,'prom'] < treated.cutoff, 1, 0)}),as.vector)
all_xyz.untreated.promContact <- do.call(rbind, all_xyz.untreated.promContact)
all_xyz.treated.promContact <- do.call(rbind, all_xyz.treated.promContact)
# Sort matrices by rowSums
all_xyz.untreated.promContact <- all_xyz.untreated.promContact[order(rowSums(all_xyz.untreated.promContact)),]
all_xyz.treated.promContact <- all_xyz.treated.promContact[order(rowSums(all_xyz.treated.promContact)),]

### Make heatmaps
pheatmap(all_xyz.untreated.promContact, cluster_rows = F, cluster_cols = F, color = c("white", "red"), labels_col = 1:30)
