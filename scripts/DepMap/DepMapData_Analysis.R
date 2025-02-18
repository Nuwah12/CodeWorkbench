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
model.idx <- read.csv('DepMap_data/Model.csv',sep=',')

### 2/17/25: Changing figure such that there are two boxes: one for B-cell lymphomas ('Mature B-Cell Neoplasm') and all others
get.essScore.forGene <- function(gene, essScores, model.index, negate=FALSE)
{
  model.bLymph <- model.index[model.index$OncotreePrimaryDisease=="Mature B-Cell Neoplasms" | model.index$OncotreeSubtype=="Mature B-Cell Neoplasms",]
  model.allOther <- model.index[model.index$OncotreePrimaryDisease!="Mature B-Cell Neoplasms" & model.index$OncotreeSubtype!="Mature B-Cell Neoplasms" & model.index$OncotreeCode!="BLL",]
  
  gene.col <- grep(gene, colnames(essScores))
  print(gene.col)
  ess.bLymph <- essScores[essScores$X%in%model.bLymph$ModelID,c(1, gene.col)]
  colnames(ess.bLymph) <- c("ModelID", paste0(gene,".ess"))
  ess.allOther <- essScores[essScores$X%in%model.allOther$ModelID,c(1, gene.col)]
  colnames(ess.allOther) <- c("ModelID", paste0(gene,".ess"))
  
  ess.bLymph.complete <- merge(ess.bLymph, model.bLymph, on=c("ModelID"))
  ess.allOther.complete <- merge(ess.allOther, model.allOther, on=c("ModelID"))
  
  # Combine
  ess.bLymph.df <- data.frame('ess.scores'=ess.bLymph, 'bLymph'=TRUE)
  ess.allOther.df <- data.frame('ess.scores'=ess.allOther, 'bLymph'=FALSE)
  ess.all <- rbind(ess.bLymph.df, ess.allOther.df)
  
  bLymph.ebf1.ess <- ess.all[ess.all$bLymph,"ess.scores"]
  allOther.ebf1.ess <- ess.all[!ess.all$bLymph,"ess.scores"] 
  
  if(gene=="EBF1"){
    ess.all <- ess.all[,c(1,2,4)]
  }
  
  colnames(ess.all) <- c("ModelID", paste0(gene,".ess"), "bLymph")
  
  if(negate){
    ess.all[,paste0(gene,".ess")] <- ess.all[,paste0(gene,".ess")] * -1
  }
  
  return(ess.all)
}
# EBF1
ebf1.ess <- get.essScore.forGene('EBF1', ess.raw, model.idx, negate=T)
# FOXO1
foxo1.ess <- get.essScore.forGene('FOXO1', ess.raw, model.idx, negate=T)
# IRF4
irf4.ess <- get.essScore.forGene('IRF4', ess.raw, model.idx, negate=T)
irf4.full <- merge(irf4.ess, model.idx, on="ModelID")
# PAX5
pax5.ess <- get.essScore.forGene('PAX5', ess.raw, model.idx, negate=T)

pdf("250217_EBF1_essentialityScores_BCellLymphoma-vs-allOther_DepMap_noBLL_negatedScore.pdf")
ggplot(data=ebf1.ess, aes(x=bLymph, y=EBF1.ess)) +
  geom_boxplot(aes(fill=bLymph)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("EBF1 Chronos Score") +
  guides(fill=guide_legend(title="B-Cell Lymphoma")) +
  scale_fill_manual(labels = c("No", "Yes"), values = c("brown2", "royalblue1")) +
  xlab("")
dev.off()
print(paste("Yes:",nrow(ebf1.ess[ebf1.ess$bLymph,]),"No:",paste("Yes:",nrow(ebf1.ess[!ebf1.ess$bLymph,]))))
ks.test(ebf1.ess[ebf1.ess$bLymph,"EBF1.ess"], ebf1.ess[!ebf1.ess$bLymph,"EBF1.ess"], alternative="less")

pdf("250217_FOXO1_essentialityScores_BCellLymphoma-vs-allOther_DepMap_noBLL_negatedScore.pdf")
ggplot(data=foxo1.ess, aes(x=bLymph, y=FOXO1.ess)) +
  geom_boxplot(aes(fill=bLymph)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("FOXO1 Chronos Score") +
  guides(fill=guide_legend(title="B-Cell Lymphoma")) +
  scale_fill_manual(labels = c("No", "Yes"), values = c("brown2", "royalblue1")) +
  xlab("")
dev.off()
print(paste("Yes:",nrow(foxo1.ess[foxo1.ess$bLymph,]),"No:",paste("Yes:",nrow(foxo1.ess[!foxo1.ess$bLymph,]))))
ks.test(foxo1.ess[foxo1.ess$bLymph,"FOXO1.ess"], foxo1.ess[!foxo1.ess$bLymph,"FOXO1.ess"], alternative="less")

pdf("250217_IRF4_essentialityScores_BCellLymphoma-vs-allOther_DepMap_noBLL_negatedScore.pdf")
ggplot(data=irf4.ess, aes(x=bLymph, y=IRF4.ess)) +
  geom_boxplot(aes(fill=bLymph)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("IRF4 Chronos Score") +
  guides(fill=guide_legend(title="B-Cell Lymphoma")) +
  scale_fill_manual(labels = c("No", "Yes"), values = c("brown2", "royalblue1")) +
  xlab("")
dev.off()
print(paste("Yes:",nrow(irf4.ess[irf4.ess$bLymph,]),"No:",paste("Yes:",nrow(irf4.ess[!irf4.ess$bLymph,]))))
ks.test(irf4.ess[irf4.ess$bLymph,"IRF4.ess"], irf4.ess[!irf4.ess$bLymph,"IRF4.ess"], alternative="less")

pdf("250217_PAX5_essentialityScores_BCellLymphoma-vs-allOther_DepMap_noBLL_negatedScore.pdf")
ggplot(data=pax5.ess, aes(x=bLymph, y=PAX5.ess)) +
  geom_boxplot(aes(fill=bLymph)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("PAX5 Chronos Score") +
  guides(fill=guide_legend(title="B-Cell Lymphoma")) +
  scale_fill_manual(labels = c("No", "Yes"), values = c("brown2", "royalblue1")) +
  xlab("")
dev.off()
print(paste("Yes:",nrow(pax5.ess[pax5.ess$bLymph,]),"No:",paste("Yes:",nrow(pax5.ess[!pax5.ess$bLymph,]))))
ks.test(pax5.ess[pax5.ess$bLymph,"PAX5.ess"], pax5.ess[!pax5.ess$bLymph,"PAX5.ess"], alternative="less")
#########################################

############### EBF1 expression ###############
exp <- read.csv('DepMap_data/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv')

exp.ebf1 <- get.essScore.forGene("EBF1", exp, model.idx)
exp.ebf1$EBF1.ess <- ifelse(exp.ebf1$EBF1.ess < 0, 0, exp.ebf1$EBF1.ess)

exp.foxo1 <- get.essScore.forGene("FOXO1", exp, model.idx)
exp.foxo1$FOXO1.ess <- ifelse(exp.foxo1$FOXO1.ess < 0, 0, exp.foxo1$FOXO1.ess)

exp.irf4 <- get.essScore.forGene("IRF4", exp, model.idx)
exp.irf4$IRF4.ess <- ifelse(exp.irf4$IRF4.ess < 0, 0, exp.irf4$IRF4.ess)

exp.pax5 <- get.essScore.forGene("PAX5", exp, model.idx)
exp.pax5$PAX5.ess <- ifelse(exp.pax5$PAX5.ess < 0, 0, exp.pax5$PAX5.ess)

pdf("250218_EBF1_ExpressionTPM_BCellLymphoma-vs-allOther_DepMap_noBLL.pdf")
ggplot(data=exp.ebf1, aes(x=bLymph, y=EBF1.ess)) +
  geom_boxplot(aes(fill=bLymph), outlier.shape=NA) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("EBF1 Expression (TPM)") +
  guides(fill=guide_legend(title="B-Cell Lymphoma")) +
  scale_fill_manual(labels = c("No", "Yes"), values = c("brown2", "royalblue1")) +
  xlab("")
dev.off()

pdf("250218_FOXO1_ExpressionTPM_BCellLymphoma-vs-allOther_DepMap_noBLL.pdf")
ggplot(data=exp.foxo1, aes(x=bLymph, y=FOXO1.ess)) +
  geom_boxplot(aes(fill=bLymph), outlier.shape=NA) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("FOXO1 Expression (TPM)") +
  guides(fill=guide_legend(title="B-Cell Lymphoma")) +
  scale_fill_manual(labels = c("No", "Yes"), values = c("brown2", "royalblue1")) +
  xlab("")
dev.off()

pdf("250218_IRF4_ExpressionTPM_BCellLymphoma-vs-allOther_DepMap_noBLL.pdf")
ggplot(data=exp.irf4, aes(x=bLymph, y=IRF4.ess)) +
  geom_boxplot(aes(fill=bLymph), outlier.shape=NA) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("IRF4 Expression (TPM)") +
  guides(fill=guide_legend(title="B-Cell Lymphoma")) +
  scale_fill_manual(labels = c("No", "Yes"), values = c("brown2", "royalblue1")) +
  xlab("")
dev.off()

pdf("250218_PAX5_ExpressionTPM_BCellLymphoma-vs-allOther_DepMap_noBLL.pdf")
ggplot(data=exp.pax5, aes(x=bLymph, y=PAX5.ess)) +
  geom_boxplot(aes(fill=bLymph), outlier.shape=NA) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  ylab("PAX5 Expression (TPM)") +
  guides(fill=guide_legend(title="B-Cell Lymphoma")) +
  scale_fill_manual(labels = c("No", "Yes"), values = c("brown2", "royalblue1")) +
  xlab("")
dev.off()



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

##### 1/31/25 - Plotting EBF1 essentiality scores of MCL cell lines (JVM-2, MAVER-1, and Z-138)
ess.mcl <- ess.scores[row.names(ess.scores)%in%c("ACH-000106", "ACH-002485", "ACH-002500"),]
colnames(ess.mcl) <- c("score", "subtype")
ess.mcl$cell.line <- c("JVM-2", "MAVER-1", "Z-138")

pdf("DepMap_EBF1_essentialityScore_MCL.pdf")
ggplot(data=ess.mcl, aes(x=cell.line, y=score))+
  geom_bar(stat="identity", fill="#575454")+
  theme(axis.ticks.y = element_line(color="black"),
        axis.text.x = element_text(color="black", size=11),
        axis.text.y = element_text(color="black"),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  xlab("Cell Line")+
  ylab("EBF1 essentiality score")
dev.off()
