###############
# Analyzing DepMap data for EBF1 Project
# 8/20/24
# Noah Burget
# Data source: https://depmap.org/portal/data_page/?tab=allData
###############
rm(list=ls())
gc()
setwd('/mnt/data0/noah/analysis/EBF1_Project')
library(ggplot2)

############### EBF1 essentiality scores ###############
# Read DepMap essentiality scores and Model index
ess.raw <- read.table('DepMap_data/CRISPRGeneEffect.csv',sep=',',header=TRUE)
model.idx <- read.csv('DepMap_data/Model.csv',sep=',')[,c('ModelID','OncotreeSubtype')]

# Subset models by Lymphomas
model.lymph <- model.idx[grepl("Lymphoma",model.idx$OncotreeSubtype),] # Lymphoma models
ess.raw <- ess.raw[ess.raw$X %in% model.lymph$ModelID,] # Subset essentiality scores for only Lymphomas
models.keep <- ess.raw$X # save models to use 
model.lymph <- model.lymph[model.lymph$ModelID %in% models.keep,] # model metadata
row.names(ess.raw) <- ess.raw$X

# Subset genes we are interested in (EBF1)
ess.scores <- ess.raw[,grepl("EBF1",colnames(ess.raw))]
ess.scores$subtype <- model.lymph$OncotreeSubtype
ess.scores <- ess.scores[,-2]
ess.scores$subtype <- as.factor(ess.scores$subtype)
# Remove Primary Effusion Lymphoma, as it is not in expression dataset
ess.scores <- subset(ess.scores, subtype!="Primary Effusion Lymphoma")

pdf('EBF1_essentialityScores_LymphomaSubtypes_MCL_DepMap.pdf')
ggplot(data=ess.scores, aes(x=subtype, y=EBF1..1879.)) +
  geom_point(size = -1, aes(color=subtype)) + 
  geom_boxplot(aes(fill=subtype), show.legend=FALSE) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("EBF1 Essentiality Score") +
  xlab("") +
  labs(fill="Lymphoma", color="Lymphoma") + 
  guides(fill=guide_legend(override.aes=list(alpha = 1, shape = 22, size = 10))) 
dev.off()

############### EBF1 expression ###############
exp <- read.csv('DepMap_data/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv')
row.names(exp) <- exp$X
exp.lymph <- exp[exp$X %in% models.keep,]
exp.lymph.ebf1 <- exp.lymph[,grepl("EBF1",colnames(exp.lymph))]
exp.lymph.ebf1$model <- row.names(exp.lymph.ebf1)
exp.lymph.ebf1$subtype <- apply(exp.lymph.ebf1, 1, function(x){
  model.lymph[model.lymph$ModelID==x[3],2]
})
exp.lymph.ebf1 <- exp.lymph.ebf1[,-c(2,3)]

pdf('EBF1_expression_LymphomaSubtypes_MCL_DepMap.pdf')
ggplot(data=exp.lymph.ebf1, aes(x=subtype, y=EBF1..1879.)) +
  geom_point(size = -1, aes(color=subtype)) + 
  geom_boxplot(aes(fill=subtype), show.legend=F) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("EBF1 Essentiality Score") +
  xlab("") +
  labs(color="Lymphoma", fill="Lymphoma") +  # Use both color and fill labels
  guides(color=guide_legend(override.aes=list(fill=NA, shape = 22, size = 10)),  # Color legend with hollow points
        )
dev.off()



