library(stringr)
library(dplyr)

args = commandArgs(trailingOnly = T)

consensus.loops <- "tmp_consensus.txt"
dotscores <- args[2]
outf <- args[3]
filter.na <- args[4]

if(filter.na == "TRUE"){
  filter.na=T
}else{filter.na=F}

### Build final output file
consensus.loops <- scan(file = consensus.loops, character())
final <- as.data.frame(do.call(rbind, str_split(consensus.loops, "_")))
dotscores <- read.table(file = dotscores, sep = ",")
final$raw <- NA
final$balanced_value <- NA
final$raw <- dotscores[,1]
final$balanced_value <- dotscores[,2]
colnames(final) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "count", "balanced")

#final$raw <- mutate_all(final$raw, function(x) as.numeric(as.chatacter(x)))
#final <- subset(final, count > -1)

if(filter.na){
  print("Filtering out loops with balanced = NA")
  final <- final[complete.cases(final),]
}

#Write to file
write.table(file = outf, x = final, quote=F, row.names=F, sep = "\t")
