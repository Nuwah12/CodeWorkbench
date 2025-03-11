library(GenomicRanges)
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
########## 
#Plotting aethetics
p0 <- theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                          panel.background = element_rect(fill = "white", colour = NA),
                          panel.border = element_rect( fill = NA, colour = "gray50", size = 1),
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          text= element_text(family="sans", size=10) )
###############################################################
# Assigns promoters (TSS) and enhancers (H3K27ac peaks) to EBF1 HiCHIP anchors
# Noah Burget, 6/8/23
###############################################################

# Setup enhancer and promoter regions + gene expression
setwd("/mnt/data0/noah/analysis/033023_loopAnalysis/EBF1_HiChIP_cliques")
source('graph_libs.R')
source('modularity.R')
h3k27ac.slopped <- read.table('S01_170114_MB157_DMSO_H3K27ac_22375365_R1_001-MACS2-FDR_1e-8-SHIFT_93_summits_slopped.bed', sep = '\t', 
                              col.names = c('chrom','start','end'))
h3k27ac.slopped$id <- paste(h3k27ac.slopped$chrom,h3k27ac.slopped$start,h3k27ac.slopped$end,sep='_')
tss.slopped <- read.table('190816_hg19_ensembl_ensg_TSS_5k.bed', sep = '\t', 
                              col.names = c('chrom','start','end','gene'))
rnaseq <- read.table('SOX9KO1_FDR0.01_log2FC0.50_all_sorted.tsv', sep = '\t',
                     header = T)
connections <- read.table('171006_MB157_SMC1_MACS2_0.05_TAD_0.05_PET_4.mango',sep='\t',header=T)[,c(1:7,10)]
#connections <- DMSO.anchor.keep[,c("p_chr","p_start","p_end","q_chr","q_start","q_end","pair.id","score")]
##########
# Function to annotate anchors of a 3D genomic assay 
# Parameters:
# tss.slopped - data.frame of tss positions (chrom,start,end,ensembl_id)
# h3k27ac.slopped - data.frame of h3k27ac peaks for enhancer marks (chrom,start,end,id)
# rnaseq - data.frame of RNAseq results (ensembl_id, RPKM)
# connections - data.frame of interacting anchor positions ('chrom1','start1','end1','chrom2','start2','end2','id','PET_count')
##########
annotate.anchors <- function(tss.slopped, h3k27ac.slopped, rnaseq, connections, PET_thresh=4, RPKM.thresh=1){
  # Enforce column names
  colnames(connections) <- c('chrom1','start1','end1','chrom2','start2','end2','id','PET')
  colnames(rnaseq) <- c('id','initial.expression')
  
  # Filter connections by PET count
  connections <- subset(connections, PET >= 4)
  
  # Filter genes out of TSS list that are not in RNAseq
  rnaseq <- subset(rnaseq, initial.expression >= 1)
  tss.slopped <- subset(tss.slopped, gene %in% rnaseq$id)

  # Identify which H3K27ac peaks are in TSS regions of expressed genes
  h3k.granges <- makeGRangesFromDataFrame(h3k27ac.slopped)
  tss.expressed.granges <- makeGRangesFromDataFrame(tss.slopped)
  h3k.slopped.enhancers <- as.data.frame(subsetByOverlaps(h3k.granges, tss.expressed.granges, invert = TRUE))
  h3k.slopped.enhancers$id <- paste(h3k.slopped.enhancers$seqnames, h3k.slopped.enhancers$start,
                                    h3k.slopped.enhancers$end, sep = "_")
  # Now, we have a list of enhancers that are NOT in an expressed TSS, and a list of expressed TSS
  
  # Make GRanges for p anchors and q anchors
  granges.p <- makeGRangesFromDataFrame(connections[,c(1,2,3,7)], seqnames.field = 'chrom1', start.field = 'start1', end.field = 'end1')
  granges.q <- makeGRangesFromDataFrame(connections[,c(4,5,6,7)], seqnames.field = 'chrom2', start.field = 'start2', end.field = 'end2')
  
  # Overlap TSSs and Enhancers with p and q anchors
  tss.granges <- makeGRangesFromDataFrame(tss.slopped)
  h3k27ac.granges <- makeGRangesFromDataFrame(h3k.slopped.enhancers)
  
  p.enh <- as.data.frame(findOverlaps(granges.p, h3k27ac.granges))
  p.tss <- as.data.frame(findOverlaps(granges.p, tss.granges))
  q.tss <- as.data.frame(findOverlaps(granges.q, tss.granges))
  q.enh <- as.data.frame(findOverlaps(granges.q, h3k27ac.granges))
  
  p.tss$feature <- tss.slopped[p.tss$subjectHits, "gene"]
  p.enh$feature <- h3k.slopped.enhancers[p.enh$subjectHits, "id"]
  q.tss$feature <- tss.slopped[q.tss$subjectHits, "gene"]
  q.enh$feature <- h3k.slopped.enhancers[q.enh$subjectHits, "id"]
  
  # all p and q anchors
  p <- rbind(p.tss, p.enh)
  q <- rbind(q.tss, q.enh)
  # Remove regions that are assigned to 2 different features 
  p <- subset(p, !duplicated(queryHits))
  q <- subset(q, !duplicated(queryHits))
  
  # match p and q anchors by index
  anchors.annotated <- merge(p, q, by = 'queryHits')
  anchors.annotated <- anchors.annotated[,-1]
  # dump duplicate interactions
  anchors.annotated <- distinct(anchors.annotated)
  # dump self-edges
  anchors.annotated <- subset(anchors.annotated, feature.x!=feature.y)
  
  # Clean up table
  anchors.annotated <- anchors.annotated[,-c(1,3)]
  colnames(anchors.annotated) <- c('from','to')
  
  return(anchors.annotated)
}
graph_communities_sizes <- function(annotated.anchors, cluster_resolution = 1, sh_out_name='super_hubs', seed=100){
  set.seed(seed=seed)
  network <- graph_from_data_frame(d = annotated.anchors, directed = FALSE)
  # Get communities
  communities <- getClustersCommunity(cluster_louvain(network, resolution=cluster_resolution))
  # Get community sizes
  community_sizes = lapply(communities, function(x) getClusterSize(network, x))
  
  # Get modularity by community
  mat <- getSimMat(network)
  comm_mod <- lapply(communities, function(x) getModularityContribution(mat, x))
  mod.by.comm = data.frame(comm.number = 1:length(community_sizes), 
                           modularity = unlist(comm_mod))
  mod.by.comm$rank <- rank(mod.by.comm$modularity, ties.method = 'random')
  
  sizes = data.frame(comm.number = 1:length(community_sizes), 
                     size = unlist(community_sizes))
  # Add ranking
  sizes$rank <- rank(sizes$size, ties.method = 'random')
  
  # df of community members with their community index #
  community.members <- data.frame(name=cbind(unlist(lapply(communities,names))),
                                  comm.num=rep(seq_along(communities),lapply(communities,length)))
  super.hubs <- find_SHs_from_size(sizes, sh_out_name)
  super.hubs <- read.table(paste0(sh_out_name,'SH_plot.tsv'), sep = '\t', header = TRUE)
  return(list("sizes"=sizes,"community.members"=community.members,"super_hubs"=super.hubs))
}
make_labeled_output <- function(annotated, membership, superhubs=NA, outname='annotated_withHubNum'){
  colnames(annotated) <- c('from','to')
  annotated$from.hub <- membership[match(annotated$from,membership$name),'comm.num']
  annotated$to.hub <- membership[match(annotated$to,membership$name),'comm.num']
  print(annotated)
  annotated <- annotated[annotated[,'from.hub']==annotated[,'to.hub'],]
  
  annotated <- subset(annotated,annotated$from.hub %in% subset(superhubs,hubClass=='SH')$comm.num)
  
  write.table(annotated[,-4], quote=F, col.names=F, row.names=F, file=paste0(outname,'.tsv'), sep='\t')
}

##########
# Compares edges between hubs
# Input format: hubNum1 hubNum2 , 
# overlapping hubs by genomic coordinates
##########
compare_hub_edges <- function(hubs, hub1.bedpe, hub2.bedpe){
  print(paste0(hub1.bedpe,'/community_hub',hubs[1],'.bedpe'))
  hub1 <- read.table(paste0(hub1.bedpe,'/community_hub',hubs[1],'.bedpe'),sep='\t',header=T)[,c('old_from','old_to')]
  hub2 <- read.table(paste0(hub2.bedpe,'/community_hub',hubs[2],'.bedpe'),sep='\t',header=T)[,c('old_from','old_to')]

  # Make graphs from both hubs
  hubg1 <- graph_from_data_frame(hub1, directed=FALSE)
  hubg2 <- graph_from_data_frame(hub2, directed=FALSE)
  
  hubg1.edges <- gsize(hubg1)
  
  # Get shared edges
  shared <- intersect(as_data_frame(hubg1),as_data_frame(hubg2))
  len <- nrow(as.data.frame(shared))
  print(len/hubg1.edges)
  #print(length(unlist(shared)))
  return(len/hubg1.edges)
}
louvain_clustering <- function(edges, hubNum, antibody, resolution=1){
  colnames(edges) <- c('from','to','hub.num')
  edges <- subset(edges, hub.num==hubNum)
  t <- graph_communities_sizes(edges[,c(1,2)],
                               sh_out_name=paste0('EBF1HiChIP_',antibody,'_superHub_',hubNum,'_reClustering_res',resolution),
                               cluster_resolution=resolution)
 return(t)
}
#########################################################
# Plotting
#########################################################

# Plot by community size
ggplot(ebf1.smc1.hubSizes, aes(x = rank, y = size))+
  geom_point()+
  xlab("Community Ranking")+
  ylab("Community Size")+#ggplot(sizes_reorder, aes(x = rank, y = size, color = highlight))+
  geom_point(show.legend = F)+
  ggtitle("Communities Ranked by Size")+
  theme_bw() + p0
ggsave("071323_EBF1_HiChIP_overlapSMC1Peak_Granta519_iGraph_communities_rankedBySize.pdf")

# Plot by modularity contribution
ggplot(ebf1.smc1.modByComm, aes(x = rank, y = modularity))+
  geom_point()+
  xlab("Community Ranking")+
  ylab("Community Modularity Contribution")+
  ggtitle("Communities Ranked by Modularity Contribution")+
  theme_bw() + p0
ggsave("071323_EBF1_HiChIP_overlapSMC1Peak_Granta519_iGraph_communities_rankedByModularity.pdf")

# Plot by size with certain hub number highlighted
hub.num <- 881
sizes$highlight <- FALSE
sizes[which(sizes$comm.number == hub.num), 'highlight'] <- TRUE

# Reorder so point of interest won't get covered by other points
sizes_reorder <- sizes[-hub.num,]
hub.row <- sizes[hub.num,]
sizes_reorder[nrow(sizes_reorder),] <- hub.row

ggplot(sizes_reorder, aes(x = rank, y = size, color = highlight))+
  geom_point(show.legend = F)+
  xlab("Community Ranking")+
  ylab("Community Size")+
  ggtitle("Communities Ranked by Size")+
  theme_bw()+
  scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
  geom_hline(yintercept = 52, linetype = 'dashed')+
  geom_vline(xintercept = 1172, linetype = 'dashed')+
  p0 
ggsave("061223_EBF1_HiChIP_Granta519_iGraph_communities_rankedBySize_hub881_sorted.pdf")









