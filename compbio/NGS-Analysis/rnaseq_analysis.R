library(stringr)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(plotgardener)
setwd("/mnt/data0/noah/analysis/misc-analysis/compbio/NGS-Analysis")
source('../DiffExp/normalize_counts.R')

COUNT.DIR <- "/mnt/data0/noah/processing/dockerize-workflows/workflows/rna_seq/outputs"

rna.outs <- list.files(COUNT.DIR, pattern=".counts", full.names=TRUE)
rna.names <- list.files(COUNT.DIR, pattern=".counts", full.names=FALSE)

rna.counts <- setNames(lapply(rna.outs, function(x){
  read.table(x, sep='\t', header=TRUE, col.names=c("gene",".",".",".",".","length","counts"))[,c(1,6,7)]
}), rna.names)

rna.rpkm <- lapply(rna.counts, rpkm)
