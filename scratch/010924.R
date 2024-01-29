setwd("/mnt/data0/noah/analysis/misc-analysis")
source('DiffExp/DESeq2.R')

counts <- list.files('/mnt/data0/noah/testing/dockerize-workflows/workflows/steps/quantify/out', full.names=T)
ctrl.1 <- read.table(counts[1],sep='\t')[,c(1:4)]
ctrl.2 <- read.table(counts[2],sep='\t')[,c(1:4)]
ko.1 <- read.table(counts[3],sep='\t')[,c(1:4)]
ko.2 <-read.table(counts[4],sep='\t')[,c(1:4)]
all.chip <- cbind(ctrl.1, ctrl.2[,4], ko.1[,4], ko.2[,4])
row.names(all.chip) <- paste0(all.chip[,1],'_',all.chip[,2],'_',all.chip[,3])
all.chip <- all.chip[,-c(1:3)]
colnames(all.chip) <- c('ctrl_rep1','ctrl_rep2','ko_rep1','ko_rep2')
rm(counts, ctrl.1, ctrl.2, ko.1, ko.2)

colData <- data.frame(row.names = colnames(all.chip), 'condition' = c('control','control','KO','KO'))

FCthresh <- 0.5
pvalThresh <- 0.01

r <- as.data.frame(do.deseq2(counts=all.chip, colData=colData, comparison=c('condition', 'KO','control')))
r$sig <- ifelse((r$log2FoldChange <= -FCthresh | r$log2FoldChange >= FCthresh) & (r$padj <= pvalThresh), TRUE, FALSE)
r$diff <- ifelse(r$sig, ifelse(r$log2FoldChange <= -FCthresh, 'Down in KO', 'Up in KO'), 'No Change')

down.ko <- nrow(subset(r, diff=='Down in KO'))
up.ko <- nrow(subset(r, diff=='Up in KO'))

pdf('010924_MB1576Cas9_H3K27ac-CnR_KO-over-Control_volcano.pdf')
ggplot(data=r, aes(x=log2FoldChange, y=-log10(padj), col=diff))+
  geom_point(alpha=0.75)+
  annotate("text", label="1209", x=-1.5, y=65, size=5)+
  annotate("text", label="613", x=1.5, y=65, size=5)+
  scale_color_manual(values=c('red','gray','blue'))+
  geom_hline(yintercept=-log10(pvalThresh), linetype="dashed")+
  geom_vline(xintercept=c(-FCthresh, FCthresh), linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.50),
        axis.line = element_line(size = 0.75, color = "black", linetype=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "transparent"))+
  ggtitle('SOX9 KO/Control H3K27ac CnR')+
  scale_x_continuous(name='Log2 Fold Change', breaks = c(-3:3))+
  scale_y_continuous(name='-log(Adj. p-value)', breaks = seq(0, 60, by=10))
dev.off()
