# note: promoters were excluded from *TE.bed
#!/usr/bin/env Rscript

rm(list=ls())
library('GenomicRanges')
library('gplots')
library('pheatmap')
library("reshape")
library("plyr")
library('edgeR')
library('ggplot2')

#****************************************************************
# Passing argument
args = commandArgs(trailingOnly=TRUE)
 if (length(args)==0) {
   stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#*****************************************************************
# NOTE: this has to be changed for revision
#working_folder="/Volumes/D/solexaSuperEnhancer/idea08_GSE50466_KaiGe_eLife/SE_new/Babak_test/"
working_folder=args[1]

#*****************************************************************
# folder for all related scripts
PIPELINE="/mnt/data0/yeqiao/analysis/pipeline/SuperEnhancer"

#*****************************************************************
# To define genes in Granges
#*****************************************************************
TSS_minus = 5000
TSS_plus = 5000
# Promoter Coordinates in RefSeq
RefSeq = read.table(file.path(PIPELINE,"hg19_RefSeq_SE.bed"),
                    header=F,stringsAsFactors=F,sep='\t')

RefSeq_gr <- GRanges(seqnames= Rle(RefSeq[,1]),ranges=IRanges(RefSeq[,2],RefSeq[,3]),strand=RefSeq[,6],intens=RefSeq[,4])

RefSeqPos1_gr <- RefSeq_gr[strand(RefSeq_gr)=='+',]
RefSeqPos2_gr <- GRanges(seqnames= seqnames(RefSeqPos1_gr),
                         ranges=IRanges(start(RefSeqPos1_gr)-TSS_minus,
                                        start(RefSeqPos1_gr)+TSS_plus),strand=strand(RefSeqPos1_gr),
                         intens=elementMetadata(RefSeqPos1_gr)[,"intens"])
# 
RefSeqNeg1_gr <- RefSeq_gr[strand(RefSeq_gr)=='-',]
RefSeqNeg2_gr <- GRanges(seqnames= seqnames(RefSeqNeg1_gr),
                         ranges=IRanges(end(RefSeqNeg1_gr)-TSS_plus,end(RefSeqNeg1_gr)+TSS_minus),
                         strand=strand(RefSeqNeg1_gr),
                         intens=elementMetadata(RefSeqNeg1_gr)[,"intens"])

RefSeq_excludePromoter <- c(RefSeqPos2_gr,RefSeqNeg2_gr)

#******************************************************************************************************************
#This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
#******************************************************************************************************************
calculate_cutoff <- function(inputVector, drawPlot,save_file,...){
  inputVector <- sort(inputVector)
  inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum)
  #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
  
  if(drawPlot){  #if TRUE, draw the plot
    pdf(save_file, useDingbats = F)
    plot(1:length(inputVector), inputVector,type="p",...)
    b <- y_cutoff-(slope* xPt)
    abline(v= xPt,h= y_cutoff,lty=2,col=8)
    #lines(xPt,y_cutoff,col=1)
    points(xPt,y_cutoff,pch=16,cex=0.1,col=2)
    abline(coef=c(b,slope),col=2)
    title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
    axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
    dev.off()
  }
  return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#******************************************************************************************************************
#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
#******************************************************************************************************************
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}

#******************************************************************************************************************
# major function
#******************************************************************************************************************
find_SEs_fromMerged <- function(merged_file)
#function(working_folder, filedir, merged_file)
{#for outputs
  SEdir=file.path(working_folder,"SEs")
  TEdir=file.path(working_folder,"TEs")
  
  merged_data = read.table(file.path(filedir,merged_file),
                           header=F,stringsAsFactors=F,sep='\t')
  
  All_Enhancers_gr <- sort(GRanges(seqnames= Rle(merged_data[,1]),
                                    ranges=IRanges(merged_data[,2],merged_data[,3])))
  
  
  output_filename=strsplit(merged_file,"merged.bedgraph")
 # merged_data[,4]=log10(merged_data[,4])
  merged_data[,4] = (merged_data[,4])/max(merged_data[,4])
  sorted_merged_data = merged_data[(sort(merged_data[,4],index.return = T))$ix,]
  
  inputVector = sorted_merged_data[,4]
  cutoff_options <- calculate_cutoff(inputVector, 
                                     drawPlot='TRUE',
                                     paste(output_filename,'SEplot.pdf',sep=''),
                                    # xlab=paste(rankBy_factor,'_enhancers'),
                                     #ylab=paste(rankByFactor,' Signal','- ',wceName),
                                     lwd=2,col=4)
  superEnhancerRows <- which(inputVector > cutoff_options$absolute)
  # find a1 is the question for super enhancers
  
  OBS <- sorted_merged_data
  SE_merged_data = OBS[OBS[,4]> cutoff_options$absolute,]
  SE_merged_data_gr <- sort(GRanges(seqnames= Rle(SE_merged_data[,1]),
                                    ranges=IRanges(SE_merged_data[,2],SE_merged_data[,3]),
                                    intens=SE_merged_data[,4]))
  
  temp = countOverlaps(All_Enhancers_gr,SE_merged_data_gr)
  TE_gr = All_Enhancers_gr[temp==0,]
  
  ListofPeaks <- TE_gr
  Intermediate <- countOverlaps(ListofPeaks,RefSeq_excludePromoter)
  ListofPeaks_nopromoters <- ListofPeaks[Intermediate == 0, ]
  mytable = data.frame(as.vector(seqnames(ListofPeaks_nopromoters)),
                       start(ListofPeaks_nopromoters), end(ListofPeaks_nopromoters))
                           
  write.table(mytable, file.path(TEdir, paste(output_filename,"TE.bed",sep="")),
                                       row.names=F,quote=F,sep='\t',col.names=F)
                           
  mytable = data.frame(as.vector(seqnames(SE_merged_data_gr)),
                       start(SE_merged_data_gr), end(SE_merged_data_gr),
                       elementMetadata(SE_merged_data_gr)[,"intens"])
  
  write.table(mytable, file.path(SEdir, paste(output_filename,"SE.bedgraph",sep="")),
              row.names=F,quote=F,sep='\t',col.names=F)
  
  
  return((SE_merged_data_gr))
}

#-----------------------------------------------------------------------------
# find all TEs and SEs
#-----------------------------------------------------------------------------

mydir=file.path(working_folder,"SEs")
filedir=file.path(working_folder,"TempFiles")
setwd(mydir)

myfiles=dir(filedir)
merged_files=myfiles[grep("_merged.bedgraph",myfiles)]
out_GR = lapply(merged_files,find_SEs_fromMerged)

