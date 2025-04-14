###################
# Functions for manipulating dataframes/matrices
# Noah Burget
###################

#' Average together every n columns in dataframe df
#' @param df - the DataFrame to be operated on
#' @param n - the number of consecutive columns to be averaged
#' 
#' @returns averaged dataframe
average_every_n <- function(df, n) 
{
  num_cols <- ncol(df)
  message(paste0("Averaging every ", n," columns on a dataframe with ", num_cols, " columns."))
  
  if ((num_cols %% n) != 0){stop("n is not a multiple of ncol(n).")}
  
  # Determine how many groups of n we have
  group_starts <- seq(1, num_cols, by = n)
  
  averaged <- lapply(group_starts, function(start) 
  {
    end <- min(start + n - 1, num_cols)
    rowMeans(df[, start:end, drop=FALSE], na.rm = TRUE)
  })
  
  averaged <- as.data.frame(averaged)
  colnames(averaged) <- paste0("avg_", group_starts, "_", group_starts+1)
  
  return(averaged)
}
