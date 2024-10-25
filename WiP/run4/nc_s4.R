####################
# Script for analyzing polymers from polymer simulation
# Noah Burget
####################
library(robustbase)
library(pheatmap)
setwd('/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims/MYC/3D_simulation/new_comparison/run4')
#source('/mnt/data0/noah/analysis/misc-analysis/genomic_analysis/ORCA/centrality.R')
#source('/mnt/data0/noah/analysis/misc-analysis/genomic_analysis/ORCA/radius_of_gyration.R')

# Plot matrix heatmap
mat <- read.table('matrix.txt')
labs <- rep('',900); labs[165] <- 'Bound.'; labs[180] <- 'E1.1'; labs[210] <- 'E1.2'; labs[315] <- 'E2'; labs[405] <- 'B1'; labs[555] <- 'B3'; labs[495] <- 'B2'; labs[735] <- 'MYC' 
color <- colorRampPalette(c("#FFFFFF","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000)
pdf('Granta519_MYC_simulation_randomLoading_NO-EBF1-Block.pdf')
pheatmap(mat, cluster_cols = F, cluster_rows = F, border_color = NA, color = color, breaks =  c(seq(0, 8000, length.out=1000)), labels_row = labs, labels_col = labs)
dev.off()

# Read in configurations
files <- list.files('confs_txt')
confs <- lapply(files, function(x) {
  read.table(paste0('confs_txt/',x), sep = " ", skip = 1)
})
# centrality
cent <- trace.centrality(confs)
cent.median <- colMedians(cent)
plot(cent.median, type="l",xlab="Monomer Position", ylab="Dist. from Center", main="Scenario #10.4 Centrality")
abline(v = 195, col="red")
abline(v = 315, col="red")
abline(v = 735, col="red")
abline(v = 405, col="blue")
abline(v = 495, col="blue")
