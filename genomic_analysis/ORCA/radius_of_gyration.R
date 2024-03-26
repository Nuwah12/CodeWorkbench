##########
# Func for calculating the squared Radius of gyration for a given chromatin trace's X,Y,Z components
# Noah Burget
# 1/5/24
##########

radius.of.gyration <- function(dat){
  dat <- lapply(dat, as.matrix)
  # get vector of distances from each (X,Y,Z) to the center
  dist_to_center_per_probe <- function(matrix){
    center <- colMeans(matrix)
    d <- apply(matrix,1,function(x){
      return((center[1] - x[1])^2 + (center[2] - x[2])^2 + (center[3] - x[3])^2)
    })
    return(d)
  }
  radii.of.gyration <- lapply(dat, dist_to_center_per_probe)
  radii.of.gyration <- lapply(radii.of.gyration,mean)
  return(lapply(radii.of.gyration,sqrt))
}
