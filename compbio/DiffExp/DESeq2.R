###############
# Func for performing Differential Expression testing using DESeq2 on arbitrary samples
# Noah Burget
# 1/8/24
###############
do.deseq2 <- function(counts, colData, comparison){
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast=comparison)
  return(res)
}

make.coldata <- function(samples, conditions, conditionColName="condition") {
  if (length(samples) != length(conditions)) {stop("Samples and conditions must be same length")}
  d <- data.frame(row.names = samples, conditionColName = conditions)
  return(d)
}
