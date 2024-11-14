####################
# Script for making a table of monomers x linear genome features
# Noah Burget
# 11/13/24
####################
rm(list=ls())
gc()

MONOMERSIZE <- 1000

setwd('/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims/MYC/monomerFeatureTable')
library(GenomicRanges)
library(AnnotationHub)
library(plyranges)
library(stringr)

### Read in peaks
atac.peaks <- read.table('peaks/221212_Granta519EBF1KI_cl27_0hr_ATACseq_merged.bed',sep='\t',col.names=c('chrom','start','end'))
ctcf.peaks <- read.table('peaks/s01_240610_Granta519EBF1KI_cl27_0hr_CTCF_merged_peaks.bed',sep='\t',col.names=c('chrom','start','end'))
h3k27ac.peaks <- read.table('peaks/S01_Granta519EBF1KI_cl27_0hr_H3K27ac_merged_peaks.bed',sep='\t',col.names=c('chrom','start','end','id'))[,c(1:3)]
smc1.peaks <- read.table('peaks/S24_230111_Granta519EBF1KI_clone27_0hr_SMC1_merged.bed',sep='\t',col.names=c('chrom','start','end','id'))[,c(1:3)]
ebf1.peaks <- read.table('peaks/221212_Granta519_cl27_EBF1_0hr_union.bed',sep='\t',col.names=c('chrom','start','end','id'))[,c(1:3)]

### CTCF motif data
ctcf.motif.dir <- read.table('/mnt/data0/noah/analysis/misc-analysis-local/ORCA/polymer_sims/CTCF_directionality/multiple_motifs/fimo_out/best_site.narrowPeak',sep='\t')[,c(1:3,6)]
colnames(ctcf.motif.dir) <- c('chrom','start','end','strand')
ctcf.motifs.granges <- makeGRangesFromDataFrame(ctcf.motif.dir)

### Define chromosomal coordinates to bin
start <- 128010000
end <- 128910000 # 900 kb region
monomers <- seq(from=start, to=end, by=MONOMERSIZE)

### Scan through monomer regions and count overlapping peaks
monomers <- data.frame('chrom'='chr8', start=monomers[-901], end=monomers[-901]+MONOMERSIZE)
monomer.granges <- makeGRangesFromDataFrame(monomers)

ctcf.motif.overlap <- findOverlaps(monomer.granges, ctcf.motifs.granges)
ctcf.motif.hits <- ctcf.motif.dir[unique(ctcf.motif.overlap@to),'strand']
strand <- rep(0, 900)
strand[unique(ctcf.motif.overlap@from)] <- ctcf.motif.hits

atac.granges <- makeGRangesFromDataFrame(atac.peaks)
overlap.atac.count <- as.data.frame(countOverlaps(monomer.granges, atac.granges))

ctcf.granges <- makeGRangesFromDataFrame(ctcf.peaks)
overlap.ctcf.count <- as.data.frame(countOverlaps(monomer.granges, ctcf.granges))

h3k27ac.granges <- makeGRangesFromDataFrame(h3k27ac.peaks)
overlap.h3k27ac.count <- as.data.frame(countOverlaps(monomer.granges, h3k27ac.granges))

smc1.granges <- makeGRangesFromDataFrame(smc1.peaks)
overlap.smc1.count <- as.data.frame(countOverlaps(monomer.granges, smc1.granges))

ebf1.granges <- makeGRangesFromDataFrame(ebf1.peaks)
overlap.ebf1.count <- as.data.frame(countOverlaps(monomer.granges, ebf1.granges))

overlap.table <- cbind(overlap.atac.count, overlap.ctcf.count, overlap.h3k27ac.count, overlap.smc1.count, overlap.ebf1.count, strand)

colnames(overlap.table) <- c('ATAC','CTCF','H3K27ac','SMC1', 'EBF1', 'CTCF_direction')

### Binarize the table
overlap.table[,'ATAC'] <- ifelse(overlap.table[,'ATAC']  >= 1, 1, 0)
overlap.table[,'CTCF'] <- ifelse(overlap.table[,'CTCF']  >= 1, 1, 0)
overlap.table[,'H3K27ac'] <- ifelse(overlap.table[,'H3K27ac']  >= 1, 1, 0)
overlap.table[,'SMC1'] <- ifelse(overlap.table[,'SMC1']  >= 1, 1, 0)
overlap.table[,'EBF1'] <- ifelse(overlap.table[,'EBF1']  >= 1, 1, 0)

overlap.table$CTCF_direction <- ifelse(overlap.table$CTCF_direction=='-', -1, ifelse(overlap.table$CTCF_direction=='+', 1, 0))

### Make it a bed file
overlap.table$chrom <- 'chr8'
overlap.table$start <- monomers$start
overlap.table$end <- overlap.table$start+MONOMERSIZE
overlap.table <- overlap.table[,c(7,8,9,1,2,3,4,5,6)]

### Write all the features separately so we can upload them to UCSC as genome tracks
write.table(subset(overlap.table[,c(1:3,4)], ATAC==1),'Granta519_ATAC_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(subset(overlap.table[,c(1:3,5,9)], CTCF==1),'Granta519_CTCF_withMotifOrient_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(subset(overlap.table[,c(1:3,6)], H3K27ac==1),'Granta519_H3K27ac_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(subset(overlap.table[,c(1:3,7)], SMC1==1),'Granta519_SMC1_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(subset(overlap.table[,c(1:3,8)], EBF1==1),'Granta519_EBF1_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)

