library("matrixStats")
library("dplyr")
library('ggplot2')
library('data.table')
library('R.matlab')
library('str2str')
library('parallel')
library('robustbase')
library('pheatmap')
library(parallel)

# Path to simulated conformations
mainDir <- "/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims/MYC/3D_simulation/112524_onlyEBF1Blocking_randomLoading/perturb4_120424"
setwd(mainDir)
#### functions for Rg, centrality, kmeans
source("/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims/orca_analysis_functions.R")

Nhyb <- 30
Npol <- 900

#### Read in configurations
files <- list.files('confs_txt')
confs <- lapply(files, function(x) {
  read.table(paste0('confs_txt/',x), sep = " ", skip = 1)
})

#### For each polymer, calculate pair-wise Euclidean distance. 900 x 900 matrices
all_distances <- lapply(confs, dist, diag = T, upper = T)

#### make each dist a matrix 
all_distances_adj <- lapply(all_distances, as.matrix)

aggregate_matrix <- function(mat, new_dim) {
  # Calculate the block size
  block_size <- nrow(mat) / new_dim
  if (block_size != floor(block_size)) {
    stop("Matrix dimensions must be divisible by new dimensions!")
  }
  
  # Use row and column indices to group and summarize
  aggregated <- matrix(
    tapply(mat, 
           list((row(mat) - 1) %/% block_size, (col(mat) - 1) %/% block_size), 
           mean), 
    nrow = new_dim, 
    ncol = new_dim
  )
  return(aggregated)
}

cores <- 16

agg.confs.dists <- mclapply(all_distances_adj, aggregate_matrix, new_dim=15, mc.cores=cores)

# Calculate centrality vectors for every simulated polymer
cent <- trace.centrality(confs)

neighborDist <- 10

### E1
E1_MYC_contact <- lapply(all_distances_adj, function(matt){
  E1_MYC <- matt[210:219,737:741]
  prob <- sum(E1_MYC < neighborDist)
  E1_pr <- ifelse(prob == 0, 0, 1)
  return(E1_pr)
})

#### E2
E2_MYC_contact <- lapply(all_distances_adj, function(matt){
  E2_MYC <- matt[304:308,737:741]
  prob <- sum(E2_MYC < neighborDist)
  E2_pr <- ifelse(prob == 0, 0, 1)
  return(E2_pr)
})

E1_MYC_prob <- sum(unlist(E1_MYC_contact))/10000
E2_MYC_prob <- sum(unlist(E2_MYC_contact))/10000

#### get the polymers
E1_E2_MYC <- data.table(cbind(unlist(E1_MYC_contact), unlist(E2_MYC_contact)))
E1_E2_MYC[,`:=`(freq = ifelse(V1 == 1 & V2 == 1, "Both", ifelse(V1 == 1 & V2 == 0, "E1_only", ifelse(V1 == 0 & V2 == 1, "E2_only", "Neither"))))] 

Both_poly <- which(E1_E2_MYC$freq == "Both")
Both_poly_cent <- cent[Both_poly,]
Both_cent_med <- apply(Both_poly_cent, 2, median)

Neither_poly <- which(E1_E2_MYC$freq == "Neither")
Neither_poly_cent <- cent[Neither_poly,]
Neither_cent_med <- apply(Neither_poly_cent, 2, median)

E1only_poly <- which(E1_E2_MYC$freq == "E1_only")
E1only_poly_cent <- cent[E1only_poly,]
E1only_cent_med <- apply(E1only_poly_cent, 2, median)

E2only_poly <- which(E1_E2_MYC$freq == "E2_only")
E2only_poly_cent <- cent[E2only_poly,]
E2only_cent_med <- apply(E2only_poly_cent, 2, median)

cent_plot <- data.table(Both_cent_med, Neither_cent_med, E1only_cent_med, E2only_cent_med, poly = c(1:900))

ggplot(cent_plot) +
  geom_line(aes(x = poly, y = Both_cent_med), color = "darkred") +
  geom_line(aes(x = poly, y = Neither_cent_med), color = "black") +
  geom_line(aes(x = poly, y = E1only_cent_med), color = "orange") +
  geom_line(aes(x = poly, y = E2only_cent_med), color = "darkgreen") +
  scale_x_continuous(limits = c(0,900), breaks = seq(0,900,30))


### Median centrality
cent.median <- colMedians(cent)
group_indices <- rep(1:(length(cent.median) / 30), each = 30)
aggregated <- tapply(cent.median, group_indices, median)
plot(aggregated, type='b')

cent_plot$probe <- group_indices

cent_plot.agg <- aggregate(cent_plot, FUN=median, by=list(group_indices))

ggplot(cent_plot.agg) +
  geom_line(aes(x = poly, y = Both_cent_med, color = "Both")) +
  geom_line(aes(x = poly, y = Neither_cent_med, color = "Neither")) +
  geom_line(aes(x = poly, y = E1only_cent_med, color = "E1 only")) +
  geom_line(aes(x = poly, y = E2only_cent_med, color = "E2 only")) +
  scale_x_continuous(limits = c(0, 900), breaks = seq(0, 900, 30)) +
  scale_color_manual(
    name = "Strat.",
    values = c("Both" = "darkred", 
               "Neither" = "black", 
               "E1 only" = "orange", 
               "E2 only" = "darkgreen")
  ) +
  theme(legend.position = 'right')

