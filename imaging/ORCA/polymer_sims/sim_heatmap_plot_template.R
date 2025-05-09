####################
# Script for analyzing polymers from polymer simulation
# Noah Burget
# This is a template file for making annotated contact maps from simulation output matrix.txt (10,000 conformation averaged heatmap)
# Replace all '...' strings
# A distance of '10' has been used for distance cutoff.
####################
library(robustbase)
library(pheatmap)
setwd('...') # Must contain 'matrix.txt', define working directory

# Read contact map
mat <- read.table('matrix.txt')
# ZEB2 simulation blocking regions:
labs.zeb2 <- rep('',900); labs[221] <- 'B1'; labs[408] <- 'ZEB2'; labs[593] <- 'B2' 
# MYC simulation blocking regions:
labs.myc <- rep('',900); labs[165] <- 'Bound.'; labs[180] <- 'E1.1'; labs[210] <- 'E1.2'; labs[315] <- 'E2'; labs[405] <- 'B1'; labs[555] <- 'B3'; labs[495] <- 'B2'; labs[735] <- 'MYC'; labs[568] <- 'B3 blocking start'; labs[570] <- 'B3 blocking end'
# Define color scale for heatmap
color <- colorRampPalette(c("#FFFFFF","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000)

# Plot
pdf('...') # define filename for contact map PDF
pheatmap(mat, 
         cluster_cols = F, 
         cluster_rows = F, 
         border_color = NA, 
         color = color, 
         breaks =  c(seq(0, 8000, length.out=1000)), 
         labels_row = labs.myc, 
         labels_col = labs.myc)
dev.off()