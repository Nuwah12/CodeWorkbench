##########
#Function(s) to process ChrTracer3 output file.
##########
############### NOTE ####################
# This function (process.chrtracer) assumes that the desired bridges are the ODD numbered readouts, that is 1,3,5,7,9, etc. 
#########################################
process.chrtracer <- function(file, Nhyb, toehold, coords=FALSE) {
  all_locations <- fread(file, header = T, stringsAsFactors = F)
  all_locations[,`:=`(trace_id = paste("fov", fov, "s", s, sep = "_"))] # Create an ID for the trace
  
  #### count total BRIDGE steps, bridges are the ODD numbered readouts
  to.call <- ifelse(toehold, Nhyb*2, Nhyb)
  bridge_only <- all_locations[readout%in%seq(1,to.call,2),] #bridges are odd numbered readouts
  bridge_count <- bridge_only[, .N, by = trace_id] #count each bridge
  
  Forward_CT_out <- bridge_only
  Forward_CT_out$readout <- (bridge_only$readout + 1)/2
  
  #### count number of hybes per spot per FOV for checking
  Forward_Nspot <- unique(Forward_CT_out$trace_id) # Hybridizations (unique spots)
  
  #### 1-column table of all hybes
  full_xyz <- data.table(c(1:Nhyb)) # 1-col table of all hybridizations
  
  fill_missing_with_NA <- function(i, mat){
    S1 <- mat[trace_id == i, c("x", "y", "z", "s", "readout")] # Extract x,y,z,spot num, and readout num for current trace_id
    S2 <- merge(S1, full_xyz, by.x = "readout", by.y = "V1", all.x = T, all.y = T) # Merge w/ full_xyz dataframe to get xyz at each readout
    return (S2)
  }
  Forward_filled_xyz_NA <- lapply(Forward_Nspot, fill_missing_with_NA, Forward_CT_out)
  if(coords){
    # Mark traces with NA in them as NA
    #Forward_filled_xyz_NA <- Forward_filled_xyz_NA
    # Remove NA elements
    Forward_filled_xyz <- Forward_filled_xyz_NA
    # Remove extra columns in the data.tables
    Forward_filled_xyz <- lapply(Forward_filled_xyz, function(x){
      x <- x[,readout:=NULL]
      x <- x[,s:=NULL]
      return(x)
    })
    return(Forward_filled_xyz)
  }
  else{
    #### change NA to Inf because dist() will ignore it.
    Forward_filled_xyz_Inf <- lapply(1:length(Forward_Nspot), function(f){
      ff <- Forward_filled_xyz_NA[[f]]
      ff[is.na(ff)] <- Inf
      return(ff)
    })
    Forward_all_distances <- lapply(Forward_filled_xyz_Inf, dist, diag = T, upper = T)
    
    #### Step2: correct distances calculated on Inf values  -------------------------------//
    set_Inf_to_NA <- function(j, mat1, mat2){
      ### whether hyb is missing in original table
      missing_na <- ifelse(is.na(mat1[[j]]$x), NA, 1) 
      ### individual distance matrix; is.infinite only works in matrix
      dist_s <- as.matrix(mat2[[j]])
      ### correct Inf
      dist_s[is.infinite(dist_s)] <- NA
      dist_na <- apply(dist_s, 1, function(k){k * missing_na})
      return(dist_na)
    }
    
    Forward_all_distances_adj <- lapply(1:length(Forward_Nspot), set_Inf_to_NA, Forward_filled_xyz_NA, Forward_all_distances)
    return(Forward_all_distances_adj) 
  }
}

make.distTable <- function(AllDistMatrix, fullSize=45, useDiag=F, probe.to.fix=NA) {
  if (!useDiag){
    all_dist_upper <- lapply(AllDistMatrix, function(m){
      mm <- m[upper.tri(m)]
      mm <- mm[is.na(mm) == F]
      return(mm)
    })
    
    #### full trace upper triangle distances -> one matrix
    full_dist <- do.call(rbind, lapply(all_dist_upper, function(i){if (length(i) == fullSize){return(i)}}))    
  }
  else {
    all_dist_upper <- lapply(AllDistMatrix, function(x){
      sn <- c()
      d <- (row(x) - col(x)) * -1
      r <- split(x,d)
      for(i in seq_along(r)[1:9]){
        if(i %% 2 == 1){
          sn <- c(sn, r[[i]])
        }
        else{
          sn <- c(sn, rev(r[[i]]))
        }
      }
      return(sn)
    })
    full_dist <- na.omit(do.call(rbind, lapply(all_dist_upper, function(i){if (length(i) == fullSize){return(i)}})))
  }

  return(full_dist)
}
