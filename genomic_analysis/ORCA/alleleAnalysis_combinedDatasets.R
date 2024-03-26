library(NMF)
library(zoo)
setwd('/mnt/data0/noah/analysis/ORCA_analysis/allele_clustering')
source('processChrTracer.R')
# Read in x,y,z coords for untreated and treated 
untreated.1 <- process.chrtracer('ChTracer3_output/comparisonRuns/230718_Grant519cl27_0hr_reverse_boxSize.csv', 
                                 30, T, coords=TRUE)
untreated.2 <- process.chrtracer('ChTracer3_output/comparisonRuns/230718_Granta519cl27_0hr_forward_boxsize.csv', 
                                 30, T, coords=TRUE)
untreated.3 <- process.chrtracer('ChTracer3_output/comparisonRuns/231020_Granta519cl27_untreated_4phBl_30minHyb_30step_allFits.csv', 
                                 30, T, coords=TRUE)
treated.1 <- process.chrtracer('ChTracer3_output/comparisonRuns/230718_Granta519cl27_24hdTAG_box_size.csv', 
                                 30, T, coords=TRUE)
treated.2 <- process.chrtracer('ChTracer3_output/comparisonRuns/230902_Granta519cl27_24hdTAG_MYC5p_30mHyb_4phBl_30step_allfits.csv', 
                               30, T, coords=TRUE)
untreated <- c(untreated.1, untreated.2, untreated.3)
treated <- c(treated.1, treated.2)
rm(untreated.1,untreated.2,untreated.3,treated.1,treated.2)
########## Keep traces with at least 60% non-NA probes #####
untreated <- lapply(untreated, function(x){
  if (sum(rowSums(is.na(x)) == ncol(x)) <= 12){
    return(x)
  }
  else{
    return(NULL)
  }
})
treated <- lapply(treated, function(x){
  if (sum(rowSums(is.na(x)) == ncol(x)) <= 12){
    return(x)
  }
  else{
    return(NULL)
  }
})
# Remove 
untreated <- Filter(function(x) !is.null(x), untreated)
treated <- Filter(function(x) !is.null(x), treated)

########## Fill in NA values ########
# 1. Fit missing points with adjacent points, if they exist
untreated <- lapply(untreated, function(na){apply(na, 2, na.approx, method = "linear")})
treated <- lapply(treated, function(na){apply(na, 2, na.approx, method = "linear")})
# 2. Use predictions based on linear regression
linear_predict <- function(v){
  vv <- data.table(c(1,2), v) ## two-by-two matrix
  colnames(vv) <- c("x", "y")
  model <- lm(y ~ x, data = vv)
  n <- data.table(1:Nhyb) ## I keep the first two values for checking.
  colnames(n) = "x"
  predict(model, newdata = n)
}
Nhyb <- 30
untreated <- lapply(untreated, function(pp){
     pp_length <- nrow(pp)
     if(pp_length < 30){
       y <- pp[(pp_length - 1):pp_length,]  ### the last two filled values
       yy <- apply(y, 2, linear_predict) ### predict xyz with cats and dogs
       pp_fill <- 3:(2 + (Nhyb - pp_length)) ### from the 3rd to whatever is needed to reach 30
       yyy <- data.table(rbind(pp, yy[pp_fill,]))
       return(yyy)
     }else{
       return(pp)
       }
  })
treated <- lapply(treated, function(pp){
  pp_length <- nrow(pp)
  if(pp_length < 30){
    y <- pp[(pp_length - 1):pp_length,]  ### the last two filled values
    yy <- apply(y, 2, linear_predict) ### predict xyz with cats and dogs
    pp_fill <- 3:(2 + (Nhyb - pp_length)) ### from the 3rd to whatever is needed to reach 30
    yyy <- data.table(rbind(pp, yy[pp_fill,]))
    return(yyy)
  }else{
    return(pp)
  }
})
untreated <- lapply(untreated, as.matrix)
treated <- lapply(treated, as.matrix)
########## Make dist matrices ##########
untreated.dists <- lapply(lapply(untreated, dist, diag=T, upper=T),as.matrix)
treated.dists <- lapply(lapply(treated, dist, diag=T, upper=T),as.matrix)

########## Name list
names(untreated) <- paste0('untreated.',seq(length(untreated)))
names(treated) <- paste0('treated.',seq(length(treated)))
########## Extract only column 25 (vector of distances to promoter)
untreated.dists.df <- do.call(rbind, lapply(untreated, function(x){return(x[,25])}))
treated.dists.df <- do.call(rbind, lapply(treated, function(x){return(x[,25])}))
untreated.dists.df <- as.data.frame(untreated.dists.df)
treated.dists.df <- as.data.frame(treated.dists.df)
########## Combine untreated and treated
all.traces <- rbind(untreated.dists.df, treated.dists.df)
rm(untreated.dists.df, treated.dists.df)
#all.traces.list <- c(untreated.dists, treated.dists)
########## NMF  ##########
set.seed(123)
# Estimate rank for input
nmf.res <- nmf(x = as.matrix(all.traces[,-25]), 
               rank = 2, 
               method = 'lee',
               seed = 'random')
# Extract Basis matrix 
nmf.basis <- as.data.frame(basis(nmf.res))
# Plot basis matrix
ggplot(data=nmf.basis, aes(x = V1, y = V2))+
  geom_point()
########## K-Means clustering
set.seed(123)
# Estimate optimal k
#fviz_nbclust(nmf.basis, kmeans, method='silhouette')
# KMeans clustering, do a large number of clusters (15)
kmeans.res <- kmeans(nmf.basis,15)
nmf.basis$kmeans.cluster <- kmeans.res$cluster
# Plot basis matrix by cluster
ggplot(data = nmf.basis, aes(x = V1, y = V2))+
  geom_point(aes(color = factor(kmeans.cluster)))
# Use factoextra package for better cluster viz
#pdf(paste0(str.untreated,'_NMF_clusterPlot_30steps.pdf'))
fviz_cluster(kmeans.res,
             data = nmf.basis[,-3],
             nmf.basispalette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw())
#dev.off()
########## Separate clusters 
clustered.alleles <- list()
for(i in 1:15){
  clustered.alleles[[i]] <- nmf.basis[nmf.basis$kmeans.cluster == i,]
}
########## Calculate proportions of untreated/treated in each cluster
clustered.alleles <- lapply(clustered.alleles, as.data.frame)
for(i in 1:length(clustered.alleles)){
  print(paste('Cluster',i))
  cluster.size <- nrow(clustered.alleles[[i]])
  untreated.num <- length(which(row.names(clustered.alleles[[i]]) %like% 'untreated'))
  treated.num <- cluster.size - untreated.num
  prop <- c((untreated.num/cluster.size),(treated.num/cluster.size))
  print(prop)
}
# Untreated clusters <- 1,3,4,6,8,9,10,12,13,14,15
# Treated clusters <- 2,5,7,11
########## Plot contact frequency for each cluster
source('plotContactFrequency.R')
cluster.dists <- list()
for(i in 1:15){
  cluster.dists[[i]] <- all.traces.list[rownames(clustered.alleles[[i]], match(names(all.traces.list)))]
}

########## Calculate interaction distance cutoffs ##########
getContactFreqDistanceCutoff <- function(x){
  all_distances_adj <- lapply(x, data.table)
  raw_dist_flat <- do.call(cbind, lapply(all_distances_adj, unlist))
  neighborDist <- median(raw_dist_flat[, !is.na(colSums(raw_dist_flat))],na.rm=T)
  return(neighborDist)
}

untreated.cutoff <- getContactFreqDistanceCutoff(untreated.dists)
treated.cutoff <- getContactFreqDistanceCutoff(treated.dists)
########## Binarize promoter distances (column 25) ##########
untreated <- do.call(rbind, lapply(untreated.dists, function(x){ifelse(x[,25] < untreated.cutoff, 1, 0)}))
treated <- do.call(rbind, lapply(treated.dists, function(x){ifelse(x[,25] < treated.cutoff, 1, 0)}))
# Sort by index of first 1, descending
# Get index of first 1 in each row
untreated <- cbind(untreated, apply(untreated, 1, function(x){which(x==1)[1]}))
treated <- cbind(treated, apply(treated, 1, function(x){which(x==1)[1]}))
# Sort
untreated <- untreated[order(untreated[,31], decreasing=T),]
treated <- treated[order(treated[,31], decreasing=T),]

########## Make heatmap ##########
set.seed(123)
size.diff <- nrow(untreated)-nrow(treated)
to.drop <- sample.int(nrow(untreated), size.diff)
untreated <- untreated[-to.drop,]
pdf("Granta519_Untreated_MYCpromoterProbe25Contact.pdf", width=3)
pheatmap(untreated[,-c(1:7,26:31)], scale="none", main="Untreated", color=c("white","red"),
        cluster_rows=F, cluster_cols=F, legend=F)
dev.off()
pdf("Granta_Treated_MYCpromoterProbe25Contact.pdf", width=3)
pheatmap(treated[,-c(1:7,26:31)], scale="none", main="Treated", color=c("white","red"),
         cluster_rows=F, cluster_cols=F, legend=F)
dev.off()
########## Subtract treated from untreated and plot heatmap
pdf('Granta519_Untreated-Minus-Treated_heatmap_8To25.pdf', width=2.5)
untreated.diff.treated <- untreated[,-31] - treated[,-31]
pheatmap(untreated.diff.treated[,-c(1:7,26:31)], scale="none", main="Untreated - Treated", color=c("blue","white","red"),
         cluster_rows=F, cluster_cols=F, legend=T)
dev.off()
pdf('Granta519_Untreated-PLUS-Treated_heatmap_8To25.pdf', width=2.5)
untreated.sum.treated <- untreated[,-31] + treated[,-31]
pheatmap(untreated.sum.treated[,-c(1:7,26:31)], scale="none", main="Untreated + Treated", color=c("white","red","black"),
         cluster_rows=F, cluster_cols=F, legend=T)
dev.off()

########## Plot mean distance to promoter ##########

source("plotContactFrequency.R")

untreated.contactFreq <- plotContactFrequency(untreated.dists, 30)[,11]
treated.contactFreq <- plotContactFrequency(treated.dists, 30)[,11]

# Fix probe 19 by replacing with mean of adjacent probes
untreated.contactFreq[19] <- mean(c(untreated.contactFreq[18],untreated.contactFreq[20]))
treated.contactFreq[19] <- mean(c(treated.contactFreq[18],treated.contactFreq[20]))

ggplot(data=data.frame(x=untreated.contactFreq), aes(x=1:30, y=x))+
  geom_line(color="blue", linewidth=1.5)+
  geom_line(color="red", data=data.frame(x=treated.contactFreq), aes(y=x))+
  scale_x_continuous(breaks=1:30, labels=1:30)

########## Idea: calculate ratio of 8 --> 25 distance / avg pairwise distance ##########
untreated.ratios <- unlist(lapply(untreated.dists, function(x){return(x[11,25]/mean(x[,25]))}))
treated.ratios <- unlist(lapply(treated.dists, function(x){return(x[11,25]/mean(x[,25]))}))










