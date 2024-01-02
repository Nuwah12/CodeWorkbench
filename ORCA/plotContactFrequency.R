plotContactFrequency <- function(x, Nhyb) {
  all_distances_adj <- lapply(x, data.table)
  raw_dist_flat <- do.call(cbind, lapply(all_distances_adj, unlist))
  
  # Median 'neighboring' distances
  neighborDist <- median(raw_dist_flat[, !is.na(colSums(raw_dist_flat))],na.rm=T)
  # Binarize matrix with neighboring distance as cutiff
  raw_contact_count <- ifelse(raw_dist_flat < neighborDist, 1, 0)
  # Count frequency of valid contact for 
  raw_contact_sum <- rowSums(raw_contact_count, na.rm = T)
  
  ### count total number of non-NA bins
  raw_non_NA <- ncol(raw_dist_flat) - rowSums(is.na(raw_contact_count))
  
  ### contact frequency
  raw_contact_Freq <- raw_contact_sum/raw_non_NA
  raw_contact_Freq_map <- matrix(raw_contact_Freq, nrow = Nhyb, ncol = Nhyb, byrow = T)

  return(raw_contact_Freq_map)
}
