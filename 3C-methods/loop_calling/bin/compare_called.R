args = commandArgs(trailingOnly = T)

chromosight <- args[1]
mustache <- args[2]
cool <- args[3]
overlap_type <- args[4] #overlap
assays_to_combine <- args[5] #datasets

library(UpSetR)
library(stringr)

########## Make upset plot for loop intersections across methods ##########
chromosight <- read.table(chromosight, sep = "\t", header = T)
mustache <- read.table(mustache, sep = "\t", header = T)
cool = read.table(cool, sep = "\t", header = T)
  
cool <- cool[,-1]
  
chromosight.loops <- paste(chromosight[,1],chromosight[,2],chromosight[,3],
                            chromosight[,4],chromosight[,5],chromosight[,6],
                            sep = "_")
mustache.loops <- paste(mustache[,1],mustache[,2],mustache[,3],
                        mustache[,4],mustache[,5],mustache[,6],
                        sep = "_")
cool.dots <- paste(cool[,1],cool[,2],cool[,3],
                    cool[,4],cool[,5],cool[,6],
                    sep = "_")

all.dots <- list(Cooltools_DotFinder = cool.dots,
		Mustache = mustache.loops,
		Chromosight = chromosight.loops)

#Decide which datasets we want to combine, and how we will combine them
datasets <- str_split(assays_to_combine, ",")

#Collect dots (probably a better way to do this)
chosen.dots <- c()
for (i in datasets){
  i <- as.integer(i)
  chosen.dots <- all.dots[i]
}

if (overlap_type == "intersect"){
  consensus.loops <- Reduce(intersect, chosen.dots)
} else if (overlap_type == "union"){
  consensus.loops <- Reduce(union, chosen.dots)
}

pdf(paste0("./sharedDots_preClustering_upset.pdf"), paper = "a4r")
upset(fromList(all.dots), order.by = "freq", sets.x.label = "Loops identified", 
      mainbar.y.label = "Loop intersections")
dev.off()
############################################################################
  
########## Generate consensus (union) set of loops ##########
#consensus.loops <- Reduce(union, list(chromosight.loops, mustache.loops, cool.dots))
write.table(file="tmp_consensus.txt", consensus.loops, col.names = F, row.names = F, quote = F)

