##########
# This script addresses Reviewer #1's comments for SOX9 paper
# Noah Burget 
# 12/21/23
##########
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(shadowtext)
setwd("/mnt/data0/noah/analysis/analysis_for_other_people/jr/1223_review")

### Load hub data
hubs <- read.table('data/reviewer1/clique_info_091522.tsv',sep='\t',header=TRUE)[,c(1:6,13)]

### Load SOX9 ChIP
sox9 <- read.table('data/reviewer1/s01_170929_MB157_untreted_Sox9_ChIPseq_36160142_S1-MACS2-PVALUE_1e-3-SHIFT_108_peaks.bed',sep='\t')[,c(1:3)]
colnames(sox9) <- c('chr','start','end')
sox9 <- makeGRangesFromDataFrame(sox9)

### Load H3K27ac data
h3k <- read.table('data/reviewer1/180123_MB157_H3K27ac_WO_summit.txt',sep='\t',header=TRUE)[,c(1:4)]
### Load promoter data
prom <- read.table('data/reviewer1/190816_hg19_ensembl_ensg_TSS_5k.bed',sep='\t')
# Combine
colnames(prom) <- colnames(h3k)
elements <- rbind(h3k,prom)
rm(h3k, prom)

### Mark elements as bound by SOX9 (T) or not (F)
elements.granges <- makeGRangesFromDataFrame(elements)
overlapping.elements <- queryHits(findOverlaps(elements.granges, sox9))
elements$x <- seq(1:nrow(elements))
elements$sox9.bound <- ifelse(elements$x %in% overlapping.elements, TRUE, FALSE)
sox9.bound.elements <- subset(elements, sox9.bound)$elm.id.summit

### Insert binding status to hub table
hubs$from.sox9.bound <- ifelse(hubs$from_name %in% sox9.bound.elements, TRUE, FALSE)
hubs$to.sox9.bound <- ifelse(hubs$to_name %in% sox9.bound.elements, TRUE, FALSE)

##### Only keep hubs that have at least 1 SOX9 binding site
hubs.sox9 <- hubs %>% group_by(rank) %>% filter(any(from.sox9.bound | to.sox9.bound)) %>% ungroup()
hubs <- hubs.sox9
### Separate into 2 tables: superhubs and normal hubs
super <- subset(hubs, hyperconnected)
normal <- subset(hubs, !hyperconnected)

### Make stacked bar chart for proportions, and add T/F and F/T rows
super <- as.data.frame(table(super[,c(8,9)]) / nrow(super))
super[2,3] <- super[2,3] + super[3,3]
super <- super[-3,]
super$status <- factor(c('Neither','One', 'Both'),levels=c('Neither','One','Both'))
super <- super[,-c(1,2)]
super$type <- 'Super Hub'
normal <- as.data.frame(table(normal[,c(8,9)]) / nrow(normal))
normal[2,3] <- normal[2,3] + normal[3,3]
normal <- normal[-3,]
normal$status <- factor(c('Neither','One', 'Both'),levels=c('Neither','One','Both'))
normal <- normal[,-c(1,2)]
normal$type <- 'Normal Hub'

all.props <- rbind(super, normal)
#all.props <- all.props[order(all.props$type),]

pdf('figures/SOX9-BindingProportion_superHubsVsNormalHub_byAnchor_012024.pdf')
ggplot(all.props, aes(x=type, y=Freq, fill=status))+
  geom_bar(stat="identity",position = position_fill(reverse = TRUE))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title="Proportion of Anchor Pairs Bound by SOX9", x="Hub Type", y="Proportion")+
  scale_fill_discrete(name = "SOX9 Anchor Binding")+
  geom_shadowtext(aes(label=scales::percent(round(Freq, 3))), y=c(0.24,0.665,0.92,0.23,0.685,0.94),color="black",size = 14/.pt, bg.colour = "white", bg.r = .1)+
  theme(plot.title = element_text(hjust = 0.20),
        axis.line = element_line(size = 1.2, color = "black", linetype=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "transparent"))
dev.off()
  


