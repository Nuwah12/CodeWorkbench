####################
# Script for making a table of monomers x linear genome features
# Noah Burget
# 11/13/24
####################
rm(list=ls())
gc()

MONOMERSIZE <- 1000
STEPSIZE <- 30000
MONOMERS.PER.PROBE <- STEPSIZE / MONOMERSIZE

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
probes <- seq(from=start, to=end, by=STEPSIZE)[-31]

### Scan through monomer regions and count overlapping peaks
monomers <- data.frame('chrom'='chr8', start=monomers[-901], end=monomers[-901]+MONOMERSIZE)
monomer.granges <- makeGRangesFromDataFrame(monomers)

ctcf.motif.overlap <- findOverlaps(monomer.granges, ctcf.motifs.granges)
ctcf.motif.hits <- ctcf.motif.dir[subjectHits(ctcf.motif.overlap),'strand']
ctcf.motif.overlap.list <- split(ctcf.motif.hits, ctcf.motif.overlap@from)

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

overlap.table <- cbind(overlap.atac.count, overlap.ctcf.count, overlap.h3k27ac.count, overlap.smc1.count, overlap.ebf1.count)

colnames(overlap.table) <- c('ATAC','CTCF','H3K27ac','SMC1', 'EBF1')

### Binarize the table
overlap.table[,'ATAC'] <- ifelse(overlap.table[,'ATAC']  >= 1, 1, 0)
overlap.table[,'CTCF'] <- ifelse(overlap.table[,'CTCF']  >= 1, 1, 0)
overlap.table[,'H3K27ac'] <- ifelse(overlap.table[,'H3K27ac']  >= 1, 1, 0)
overlap.table[,'SMC1'] <- ifelse(overlap.table[,'SMC1']  >= 1, 1, 0)
overlap.table[,'EBF1'] <- ifelse(overlap.table[,'EBF1']  >= 1, 1, 0)

### Make it a bed file
overlap.table$chrom <- 'chr8'
overlap.table$start <- monomers$start
overlap.table$end <- overlap.table$start+MONOMERSIZE
overlap.table <- overlap.table[,c(7,8,9,1,2,3,4,5,6)-1]

### Make CTCF w/ strands UCSC compatible 
ctcf <- overlap.table[,c('chrom','start','end','CTCF')]
ctcf$id <- paste(ctcf$chrom,ctcf$start,ctcf$end,sep='_')
ctcf <- ctcf[,c(1,2,3,6,4,5)]

### Write all the features separately so we can upload them to UCSC as genome tracks
write.table(subset(overlap.table[,c(1:3,4)], ATAC==1),'Granta519_ATAC_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(subset(ctcf, CTCF==1),'Granta519_CTCF_withMotifOrient_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(subset(overlap.table[,c(1:3,6)], H3K27ac==1),'Granta519_H3K27ac_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(subset(overlap.table[,c(1:3,7)], SMC1==1),'Granta519_SMC1_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(subset(overlap.table[,c(1:3,8)], EBF1==1),'Granta519_EBF1_simulationLocus_byMonomer_binarizedPeaks.tsv',sep='\t',row.names=F,col.names=F,quote=F)

### Find proportions of monomers w/ peaks in each probe by summing the overlap table every 30 rows
probe <- rep(1:(nrow(overlap.table) / MONOMERS.PER.PROBE), each = MONOMERS.PER.PROBE)
overlap.table.perprobe <- aggregate(. ~ probe, data = cbind(probe, overlap.table[,c(4:8)]), FUN = sum)[,-c(1)]
overlap.table.perprobe$totalPeaks <- rowSums(overlap.table.perprobe)
overlap.table.perprobe[,c('ATAC','CTCF','H3K27ac','SMC1','EBF1')] <- overlap.table.perprobe[,c('ATAC','CTCF','H3K27ac','SMC1','EBF1')] / overlap.table.perprobe$totalPeaks
overlap.table.perprobe <- lapply(overlap.table.perprobe, function(x) {
  if (is.numeric(x)) {
    x[is.nan(x)] <- 0
  } else if (is.data.frame(x)) {
    x[is.nan(x)] <- 0
  }
  return(x)
})
### write peak proportion table
write.table(overlap.table.perprobe[-6], 'Granta519_perProbe_peaksInMonomerProportion.tsv',sep='\t',row.names=F,quote=F)


####################
# Calculating coverage for each monomer region (1kb windows)
# 11/18/24 
####################
atac.cov <- read.table('Granta519_simulationMonomerBins_ATACcoverage.tsv',sep='\t')[,c(1:4)]
smc1.cov <- read.table('Granta519_simulationMonomerBins_SMC1coverage.tsv',sep='\t')[,c(1:4)]
ctcf.cov <- read.table('Granta519_simulationMonomerBins_CTCFcoverage.tsv',sep='\t')[,c(1:4)]
ebf1.cov <- read.table('Granta519_simulationMonomerBins_EBF1coverage.tsv',sep='\t')[,c(1:4)]
atac.reads <- 76213064
smc1.reads <- 109433687
ctcf.reads <- 210110958
ebf1.reads <- 245337975
# Normalize to RPM
atac.cov$V4 <- (atac.cov$V4/atac.reads)*1000000
smc1.cov$V4 <- (smc1.cov$V4/smc1.reads)*1000000
ctcf.cov$V4 <- (ctcf.cov$V4/ctcf.reads)*1000000
ebf1.cov$V4 <- (ebf1.cov$V4/ebf1.reads)*1000000

monomer.coverage.table <- data.frame('chrom'='chr8',
                                     'start'=atac.cov$V2,
                                     'end'=atac.cov$V3,
                                     'atac.cov'=atac.cov$V4,
                                     'smc1.cov'=smc1.cov$V4,
                                     'ctcf.cov'=ctcf.cov$V4,
                                     'ebf1.cov'=ebf1.cov$V4)

par(mfrow=c(4,1))

## E1
probe.to <- 6
probe.from <- 8 

view <- c((30*(probe.to-1)):(30*probe.from))

plot(monomer.coverage.table[view,'atac.cov'], type='l',ylab='ATAC')
abline(v=35, col='green')
abline(v=59, col='green')
plot(monomer.coverage.table[view,'smc1.cov'], type='l',ylab='SMC1')
plot(monomer.coverage.table[view,'ctcf.cov'], type='l',ylab='CTCF')
abline(v=31, col='red')
abline(v=34, col='red')
abline(v=60, col='red')
abline(v=70, col='red')
plot(monomer.coverage.table[view,'ebf1.cov'], type='l',ylab='EBF1')

## E2
probe.to <- 10
probe.from <- 12 

view <- c((30*(probe.to-1)):(30*probe.from))

plot(monomer.coverage.table[view,'atac.cov'], type='l',ylab='ATAC')
plot(monomer.coverage.table[view,'smc1.cov'], type='l',ylab='SMC1')
plot(monomer.coverage.table[view,'ctcf.cov'], type='l',ylab='CTCF')
abline(v=31, col='red')
abline(v=34, col='red')
plot(monomer.coverage.table[view,'ebf1.cov'], type='l',ylab='EBF1')











