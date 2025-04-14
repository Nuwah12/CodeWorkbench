library(stringr)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(plotgardener)
library(dplyr)
library(reshape2)

process.rna.counts <- function(count.dir, norm="rpkm")
{
  if (!grepl("^rpm$|^rpkm$", norm)){stop(paste0("Norm must be one of ('rpm', 'rpkm'), not ", norm))}
  
  source('/mnt/data0/noah/analysis/misc-analysis/compbio/NGS-Analysis/DiffExp/normalize_counts.R')
  setwd(count.dir)
  
  rna.outs <- list.files(count.dir, pattern=".counts", full.names=TRUE)
  rna.names <- list.files(count.dir, pattern=".counts", full.names=FALSE)
  
  rna.counts <- setNames(lapply(rna.outs, function(x){
    read.table(x, sep='\t', header=TRUE, col.names=c("gene",".",".",".",".","length","counts"))[,c(1,6,7)]
  }), rna.names)
  
  if (norm == "rpkm"){rna.norm <- lapply(rna.counts, rpkm)}
  else if (norm == "rpm"){rna.norm <- lapply(rna.counts, rpm)}
  
  # Refactor list to only include norm values
  norm_only <- setNames(
    lapply(names(rna.norm), function(n) rna.norm[[n]][[norm]]),
    names(rna.norm)
  )
  
  # Combine gene column and count list from previous statement
  combined <- data.frame(
    gene = rna.norm[[1]]$gene,
    do.call(cbind, norm_only)
  )  
  
  return(combined)
}


# Average every 2 replicates to get avg RPKM
average_every_2 <- function(df) {
  num_cols <- ncol(df)
  print(num_cols)
  col_pairs <- seq(2, num_cols, by = 2)
  averaged <- sapply(col_pairs, function(i) {
    if (i < num_cols) {
      rowMeans(df[, c(i, i + 1)], na.rm = TRUE)
    } else {
      df[, i]
    }
  })
  averaged <- as.data.frame(averaged)
  colnames(averaged) <- paste0("avg_", col_pairs)
  cbind(df[, 1, drop=FALSE], averaged)  # keep gene column if needed
}