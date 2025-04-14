####################
# Functions for common normalization techniques of genomic count data (RPM, RPKM)
# Noah Burget
# 7/16/24
####################

rpkm <- function(x){
  if(ncol(x) < 3){
    stop('Table must contain three columns in this order: Gene ID, Gene length in bp, and raw counts.')
  }
  totalReads <- sum(x[,3])
  factor <- (x[,2]/1E3)*(totalReads/1E6)
  x$rpkm <- x[,3] / factor
  return(x)
}

rpm <- function(x){
  if(ncol(x) < 2){
    stop('Table must contain two columns in this order: Gene ID, and raw counts.')
  }
  totalReads <- sum(x[,2])
  factor <- totalReads/1E6
  x$rpm <- x[,2] / factor
  return(x)  
}
