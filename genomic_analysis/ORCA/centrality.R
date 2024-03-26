####################
# Script for calculating centrality - the distances of each ORCA probe to the center of the trace
# Noah Burget
# 1/3/24
####################

# dat should be a list of XYZ coordinates with no NA values
trace.centrality <- function(dat){
	dat <- untreated
	dat <- lapply(dat, as.matrix)
	# get vector of distances from each (X,Y,Z) to the center
	dist_to_center_per_probe <- function(matrix){
	  center <- colMeans(matrix)
	  d <- apply(matrix,1,function(x){
	    return(sqrt((center[1] - x[1])^2 + (center[2] - x[2])^2 + (center[3] - x[3])^2))
	  })
	  return(d)
	}
	centrality.vectors <- lapply(dat, dist_to_center_per_probe)
	# Make centrality table (Traces X distances, y X 30)
	return(do.call(rbind, centrality.vectors))
}

# Plot mean centrality of all traces
mean.centrality <- colMeans(centrality)
plot(mean.centrality, type="l")
