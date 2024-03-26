##########
# Noah Burget
# 9/6/2023
# This script will determine the most DEP (diff. expressed protein) from a matrix w/ group/cluster labels
# Originally used for AnnoSpat paper analysis
##########
setwd('/mnt/data0/noah/analysis/analysis_for_other_people/aanchal')

proteinExpr.clusterLabeled <- read.table('ctrl_metaCluster.tsv',sep='\t',header=T)[,-36] # File of cells x protein expression

clusters <- unique(proteinExpr.clusterLabeled$meta.clus)
clusters <- sort(clusters)
test.results.foldChange <- data.frame()
test.results.wilcox <- data.frame()
for(i in clusters){ # separate cells in target cluster from all other cells
  target.cluster <- subset(proteinExpr.clusterLabeled, meta.clus == i)[,-36] # Get cells in current target cluster
  other.clusters <- subset(proteinExpr.clusterLabeled, meta.clus != i)[,-36] # Get cells in every other cluster
  for(j in colnames(target.cluster)){ # For every protein measured, compare 
    print(paste0("Testing ",j," in cluster ",i," against all other clusters"))
    target <- target.cluster[,j] # target cluster expression
    other <- other.clusters[,j] # all other cluster expression
    log2FC <- log2(mean(target) / mean(other)) # Calculate log2 fold change from mean expression values
    test.results.foldChange[i+1,j] <- log2FC
    res <- wilcox.test(target, other) # Use wilcox rank-sum test to determine p-value
    test.results.wilcox[i+1,j] <- res$p.value
  }
}

test.results.wilcox[test.results.wilcox == 0] <- .Machine$double.xmin # p-values = 0 are set to lowest double val
test.results.wilcox <- -log10(test.results.wilcox)
test.results.wilcoxFoldChange <- test.results.wilcox * test.results.foldChange
