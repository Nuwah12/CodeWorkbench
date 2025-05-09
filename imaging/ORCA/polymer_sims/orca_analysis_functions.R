
library("matrixStats")
library("dplyr")
library('ggplot2')
library('data.table')
library('R.matlab')
library('str2str')
library('parallel')

#### functions for Rg, centrality, kmeans

#### ggplot cosmatics
p0 <- theme_bw() +
  theme( plot.title = element_text(hjust = 0.5),
         panel.background = element_rect(fill = "white", colour = NA),
         panel.border = element_rect( fill = NA, colour = "black"),
         legend.position="none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#### Using 7/25 and 8/25 as enhancer 1; when either is in contact, count as 1
add_NA <- function(v){
  x <- v[1]
  y <- v[2]
  if (is.na (x) & is.na (y)){
    return (NA)
  } else if (is.na (x)) {
    return (y)
  } else if (is.na(y)) {
    return (x)
  } else {
    return (x + y)
  }
}

#### function for E1 - E2 - MYC three way interactions for each allele
#### doesn't calculate frequency
E1E2_contact <- function(distMat, v1, v2, v3, distCutOff) {
  #### binary contacts of 7/8 25 
  v1_mat <- distMat[v1,]
  v1_contact <- v1_mat < distCutOff
  
  v2_mat <- distMat[v2,]
  v2_contact <- v2_mat < distCutOff
  
  #### calculate either or of v1 and v2; only when both points are missing, value is NA
  v1_v2 <- rbind(as.vector(v1_contact), as.vector(v2_contact))
  v1_v2_contact <- data.table(apply(v1_v2, 2, add_NA) > 0)
  
  #### binary contacts of 11 25
  v3_mat <- distMat[v3,]
  v3_contact <- data.table(v3_mat < distCutOff)
  
  #### combine and only keep alleles with both 
  E1_E2 <- cbind(v1_v2_contact, v3_contact)
  colnames(E1_E2) <- c("E1", "E2")
  E1_E2[,`:=`(E1E2MYC = ifelse(E1 == T & E2 == T, "Both", ifelse(E1 == T & E2 == F, "E1_only", ifelse(E1 == F & E2 == T, "E2_only", ifelse(E1 == F & E2 == F, "Neither", "NA")))))]
  
  return(E1_E2)
}

###### Func for calculating the squared Radius of gyration for a given chromatin trace's X,Y,Z components
# Noah Burget
# get vector of distances from each (X,Y,Z) to the center
# For MYC in Granta519, step 19 is always imputed by averaging 18 and 20
# Do not use this for other locus
dist_to_center_per_probe_impute_19 <- function(matrix){
  center <- colMeans(matrix)
  d <- apply(matrix,1,function(x){
    dist_raw <- sqrt((center[1] - x[1])^2 + (center[2] - x[2])^2 + (center[3] - x[3])^2)
    return(dist_raw)
  })
  d[19] <- (d[18] + d[20])/2
  return(d)
}

trace.centrality <- function(dat){
  dat <- lapply(dat, as.matrix)
  centrality.vectors <- lapply(dat, dist_to_center_per_probe_impute_19)
  # Make centrality table (Traces X distances, y X 30)
  return(do.call(rbind, centrality.vectors))
}

trace.centrality <- function(dat){
  dat <- lapply(dat, as.matrix)
  centrality.vectors <- lapply(dat, dist_to_center_per_probe_impute_19)
  # Make centrality table (Traces X distances, y X 30)
  return(do.call(rbind, centrality.vectors))
}

radius_of_gyration <- function(dat){
  dat <- lapply(dat, as.matrix)
  r_sqrt <- lapply(dat, dist_to_center_per_probe_impute_19)
  radii <- lapply(r_sqrt,mean)
  return(radii)
}

#### confidence interval
sample_CI <- function(v){
  mean_value <- mean(v)
  n <- length(v)
  standard_deviation <- sd(v)
  
  # Find the standard error
  standard_error <- standard_deviation / sqrt(n)
  alpha = 0.05
  degrees_of_freedom = n - 1
  t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
  margin_error <- t_score * standard_error
  
  # Calculating lower bound and upper bound
  lower_bound <- mean_value - margin_error
  upper_bound <- mean_value + margin_error
  
  # Print the confidence interval
  return(c(lower_bound,upper_bound))
}

#### calculate centrality of subsampling
sample_centrality <- function(i, matt, poly, size, cutt) {
  set.seed(i)
  samp <- sample(1:ncol(matt), size) ### which columns of dist matrix
  samp_contact <- E1E2_contact(distMat = matt[,samp], v1 = a1, v2 = a2, v3 = b, distCutOff = cutt) ### whether in contact
  samp_poly <- poly[samp] ### which polymers
  samp_cent <- data.table(trace.centrality(samp_poly)) ### centrality
  
  #### categorize
  samp_Both_cent <- apply(samp_cent[which(samp_contact$E1E2MYC == "Both"),], 2, median)
  samp_E1_cent <- apply(samp_cent[which(samp_contact$E1E2MYC == "E1_only"),], 2, median)
  samp_E2_cent <- apply(samp_cent[which(samp_contact$E1E2MYC == "E2_only"),], 2, median)
  samp_Neither_cent <- apply(samp_cent[which(samp_contact$E1E2MYC == "Neither"),], 2, median)
  
  return(data.table(samp_Both_cent, samp_E1_cent, samp_E2_cent, samp_Neither_cent))
}

#### calculate centrality without sampling
hub_centrality <- function(matt, poly, cutt) {
  samp_contact <- E1E2_contact(distMat = matt, v1 = a1, v2 = a2, v3 = b, distCutOff = cutt) ### whether in contact
  samp_cent <- data.table(trace.centrality(poly)) ### centrality
  
  #### categorize
  samp_Both_cent <- apply(samp_cent[which(samp_contact$E1E2MYC == "Both"),], 2, median)
  samp_E1_cent <- apply(samp_cent[which(samp_contact$E1E2MYC == "E1_only"),], 2, median)
  samp_E2_cent <- apply(samp_cent[which(samp_contact$E1E2MYC == "E2_only"),], 2, median)
  samp_Neither_cent <- apply(samp_cent[which(samp_contact$E1E2MYC == "Neither"),], 2, median)
  
  return(data.table(samp_Both_cent, samp_E1_cent, samp_E2_cent, samp_Neither_cent))
}

