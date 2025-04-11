library(stringr)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(plotgardener)
library(dplyr)
library(reshape2)
setwd("/mnt/data0/noah/analysis/misc-analysis/compbio/NGS-Analysis")
source('../DiffExp/normalize_counts.R')
source("../../scripts/plotting.R")

COUNT.DIR <- "/mnt/data0/noah/processing/dockerize-workflows/workflows/rna_seq/outputs"
EBF1.ENSG <- "ENSG00000164330"
SOX11.ENSG <- "ENSG00000176887"
CD28.ENSG <- "ENSG00000178562"
B2M.ENSG <- "ENSG00000166710"
CD5.ENSG <- "ENSG00000110448"

rna.outs <- list.files(COUNT.DIR, pattern=".counts", full.names=TRUE)
rna.names <- list.files(COUNT.DIR, pattern=".counts", full.names=FALSE)

rna.counts <- setNames(lapply(rna.outs, function(x){
  read.table(x, sep='\t', header=TRUE, col.names=c("gene",".",".",".",".","length","counts"))[,c(1,6,7)]
}), rna.names)

rna.rpkm <- lapply(rna.counts, rpkm)

# Refactor list to only include RPKM values
counts_only <- setNames(
  lapply(names(rna.rpkm), function(n) rna.rpkm[[n]][["rpkm"]]),
  names(rna.rpkm)
)

# Combine gene column and count list from previous statement
combined <- data.frame(
  gene = rna.rpkm[[1]]$gene,
  do.call(cbind, counts_only)
)

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