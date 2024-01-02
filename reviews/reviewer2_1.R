##########
# This script addresses Reviewer #2's comment #1 for SOX9 paper
# Noah Burget 
# 12/15/23
##########
gc()
library(stringr)
library(pheatmap)
library(latticeExtra)
setwd("/mnt/data0/noah/analysis/analysis_for_other_people/jr/1223_review")

## Read in expression data and Model reference
exp <- read.table('data/reviewer2/CRISPRGeneEffect.csv', sep=',', header=TRUE)
models <- read.table('data/reviewer2/Model.csv', sep=',', fill=TRUE, header=TRUE)
# Genes included in CRISPR screen
screen.genes <- c(read.table('data/reviewer2/MB157.gene_summary.txt', sep='\t', header=T)[,1], 'X') 
exp.genes <- sapply(strsplit(colnames(exp), "\\."), `[`, 1)
# !!! Not all genes from the CRISPR screen (screen.genes) are in the dataset (exp.genes) !!!
exp.rel_genes <- exp[,exp.genes %in% screen.genes]
#rm(exp, screen.genes, exp.genes)

# Map model to cancer subtype
models <- models[,c('ModelID','LegacyMolecularSubtype')]
colnames(models) <- c('X','LecacyMolecularSubtype')
exp.models <- merge(models, exp.rel_genes, by='X')

# Subset the dataframe with the cancer subtypes we want to keep
keep <- c('basal_A','basal_B','luminal','HER2_amp')
exp.models.subtypes <- subset(exp.models, exp.models$LecacyMolecularSubtype %in% keep)
rm(exp.rel_genes, models)
# order matrix based on subtype factor order
exp.models.subtypes <- exp.models.subtypes[order(as.factor(exp.models.subtypes$LecacyMolecularSubtype)),]
###### Remove columns (genes) with any NA expression values 
exp.models.subtypes <- exp.models.subtypes[ , colSums(is.na(exp.models.subtypes))==0]
subtypes <- exp.models.subtypes$LecacyMolecularSubtype
exp.models.subtypes <- as.data.frame(t(exp.models.subtypes[,-c(1,2)]))
row.names(exp.models.subtypes) <- make.unique(sapply(strsplit(row.names(exp.models.subtypes), "\\."), `[`, 1))

# Make heatmap
my.breaks <- c(seq(-1,1,by=0.01))
redblue <- round(0.475 * length(my.breaks))
ramp1 <- colorRampPalette(c("blue", "white"))(redblue)
ramp2 <- colorRampPalette(c("white", "red"))(redblue)
my.colors <- c(ramp1, rep('#FFFFFF', (length(my.breaks)-redblue*2)-1),ramp2[-1])
## Order matrix by most negative basal scores
exp.models.subtypes$basal.sum <- rowSums(exp.models.subtypes[c(grep('basal',subtypes))])
exp.models.subtypes <- exp.models.subtypes[order(exp.models.subtypes$basal.sum),]

### SUBSET: keep only 'essential' genes from CRISPR screen
essential.genes <- read.table('data/reviewer1/essential_and_expressed.tsv',sep='\t')
exp.models.subtypes <- exp.models.subtypes[rownames(exp.models.subtypes) %in% essential.genes$id,]

### SUBSET: keep only genes in hyperconnected hubs AND that are essential
#hyperconnected.genes <- read.table('data/reviewer1/essential_and_expressed_hyper.tsv',sep='\t')
#exp.models.subtypes <- exp.models.subtypes[rownames(exp.models.subtypes) %in% hyperconnected.genes$id,]

## Subset matrix for better view of essential genes (basal score > 2 | basal score < -2)
#exp.models.subtypes <- subset(exp.models.subtypes, basal.sum > 3 | basal.sum < -3)

colann <- data.frame(Subtype=exp.models$LecacyMolecularSubtype)
rowann <- data.frame(rep("", nrow(exp.models.subtypes)))

rowann[grep('SOX9', rownames(exp.models.subtypes)),1] <- 'SOX9' #####
rowann[rownames(exp.models.subtypes)=='MYC',1] <- 'MYC' #####
pdf('figures/heatmap_breastCancerEssentialGenes_CRISPRscreen_010224_essential_and_expressed.pdf')
pheatmap(exp.models.subtypes[,-23], 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         color=my.colors, 
         breaks=my.breaks, 
         show_colnames=FALSE, 
         show_rownames=TRUE,
         annotation_col = colann, 
         labels_row = rowann[,1],
         annotation_names_col = FALSE
         #treeheight_col = 0,
         #fontsize_col = 10
         )
dev.off()

##### Plot CDFs of SOX9 score distributions for basal/non-basal(other)
# Negate scores to make positive = more essential & get CDFunctions
basal <- as.matrix(sort(exp.models.subtypes['SOX9',1:12]*-1))
other <- t(as.matrix(sort( c(unlist(exp.models.subtypes['SOX9',13:22]*-1),0.20394466,0.20394466))))
props <- t(rbind(basal,other))
colnames(props) <- c('basal','other')

basal.cdf <- ecdf(basal)
other.cdf <- ecdf(other)
x.vals <- seq(min(basal, other), max(basal, other), length.out=100)
props <- data.frame(x=x.vals,basal=basal.cdf(x.vals),other=other.cdf(x.vals))
props <- (props)

pdf('SOX9EssentialityScore_BreastCancer_BasalvsOthers_123123.pdf', height=7.5, width=10)
ggplot(data=props, aes(x=x))+
  geom_line(aes(y = basal, color = "Basal"), linewidth = 0.5)+
  geom_line(aes(y = other, color = "Other"), linewidth = 0.5)+
  ylab('Proportion')+
  xlab('Essentiality Score * -1')+
  theme_minimal()+
  scale_color_manual(values = c("Basal" = "blue", "Other" = "red")) +
  labs(color = "Breast Cancer Subtype")+
  ggtitle('(Negative) Essentialiy Score of SOX9 in Basal/Other Breast Cancer Subtypes')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
  

ecdfplot(~ basal + other, data=as.data.frame(props), xlab='-1 * (Essentiality Score)', ylab='Proportion')



