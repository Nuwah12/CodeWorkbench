##########
# This script addresses Reviewer #1's comments (1.2.2) for SOX9 paper
# Noah Burget 
# 1/4/23
##########
library(GenomicRanges)
library(ggplot2)
setwd("/mnt/data0/noah/analysis/analysis_for_other_people/jr/1223_review")

### Load hub data
hubs <- read.table('data/reviewer1/clique_info_091522.tsv',sep='\t',header=TRUE)

########## Only considering loops that were LOST (up_WO) ##########
#hubs <- subset(hubs, loop.diff.status=='up_WO')[,c(1:4,13)]
hubs <- hubs[,c(1:4,13)]

### Load Notch1 ChIP
notch.on <- read.table('data/reviewer1/S07_170114_MB157_DMSO_Notch1_22375365_R1_001-MACS2-PVALUE_1e-3-SHIFT_123_peaks.bed',sep='\t')[,c(1:3)]
notch.off <- read.table('data/reviewer1/S08_170114_MB157_GSI_Notch1_22375365_R1_001-MACS2-PVALUE_1e-3-SHIFT_121_peaks.bed',sep='\t')[,c(1:3)]
colnames(notch.on) <- c('chr','start','end')
colnames(notch.off) <- c('chr','start','end')
notch.on <- makeGRangesFromDataFrame(notch.on)
notch.off <- makeGRangesFromDataFrame(notch.off)

### Load H3K27ac data
h3k <- read.table('data/reviewer1/180123_MB157_H3K27ac_WO_summit.txt',sep='\t',header=TRUE)[,c(1:4)]
### Load promoter data
prom <- read.table('data/reviewer1/190816_hg19_ensembl_ensg_TSS_5k.bed',sep='\t')
# Combine
colnames(prom) <- colnames(h3k)
elements <- rbind(h3k,prom)
rm(h3k, prom)

### Mark elements as bound by Notch1 (ON condition) (T) or not (F)
elements.granges <- makeGRangesFromDataFrame(elements)
overlapping.elements <- queryHits(findOverlaps(elements.granges, notch.on))
elements$x <- seq(1:nrow(elements))
elements$notch1.bound <- ifelse(elements$x %in% overlapping.elements, TRUE, FALSE)
notch1.bound.elements <- subset(elements, notch1.bound)$elm.id.summit

### Insert binding status to hub table
hubs$from.notch1.bound <- ifelse(hubs$from_name %in% notch1.bound.elements, TRUE, FALSE)
hubs$to.notch1.bound <- ifelse(hubs$to_name %in% notch1.bound.elements, TRUE, FALSE)

##### Only keep hubs that have at least 1 Notch1 binding site
hubs.notch <- hubs %>% group_by(rank) %>% filter(any(from.notch1.bound | to.notch1.bound)) %>% ungroup()
hubs <- hubs.notch

### Superhubs vs normal hubs
super <- subset(hubs, hyperconnected)
normal <- subset(hubs, !hyperconnected)

### Prepare proportions for plotting
super <- as.data.frame(table(super[,c(6,7)]) / nrow(super))
super[2,3] <- super[2,3] + super[3,3]
super <- super[-3,]
super$status <- factor(c('Neither','One', 'Both'),levels=c('Neither','One','Both'))
super <- super[,-c(1,2)]
super$type <- 'Super Hub'
normal <- as.data.frame(table(normal[,c(6,7)]) / nrow(normal))
normal[2,3] <- normal[2,3] + normal[3,3]
normal <- normal[-3,]
normal$status <- factor(c('Neither','One', 'Both'),levels=c('Neither','One','Both'))
normal <- normal[,-c(1,2)]
normal$type <- 'Normal Hub'

all.props <- rbind(super, normal)

pdf('figures/Notch1-BindingProportion_superHubsVsNormalHub_byAnchor_LoopsDownInGSI_012024.pdf')
ggplot(all.props, aes(x=type, y=Freq, fill=status))+
  geom_bar(stat="identity",position = position_fill(reverse = TRUE))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title="Proportion of Anchor Pairs Bound by Notch1", x="Hub Type", y="Proportion")+
  scale_fill_discrete(name = "Notch1 Anchor Binding")+
  geom_shadowtext(aes(label=scales::percent(round(Freq, 3))), y=c(0.325,0.785,0.97,0.335,0.79,0.98),color="black",size = 11/.pt, bg.colour = "white", bg.r = .1)+
  theme(plot.title = element_text(hjust = 0.20),
        axis.line = element_line(size = 1.2, color = "black", linetype=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "transparent"))
dev.off()


