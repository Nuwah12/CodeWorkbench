library(ggplot2)
source("/mnt/data0/noah/analysis/CodeWorkbench/scripts/plotting.R")

setwd("/mnt/data0/noah/analysis/CodeWorkbench/imaging/FISH/piscis/piscesPipeline")
spots <- do.call(rbind,lapply(list.files("./test_spots"), function(x){
  table <- read.table(paste0("./test_spots/",x),sep="\t")
  count <- nrow(table)
  return(data.frame(x, count))
}))
spots$cond <- ifelse(grepl("well2", spots$x), "Control", 
                     ifelse(grepl("well1", spots$x), "Promoter KD",
                            ifelse(grepl("well3", spots$x), "Enhancer 1 KD", "Enhancer 2 KD")))

ggplot(data=spots, aes(x=cond, y=count))+
  geom_boxplot()+
  min.theme()+
  xlab("ICOS Condition")+
  ylab("Total Spots / FOV")
