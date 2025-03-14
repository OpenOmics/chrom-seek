#!/usr/bin/env Rscript
####################
#
# Name: DiffBind_PCA_normalizedCounts.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: December 7, 2022
# Updated: November 6, 2024
# 
####################
#
# Purpose: To get PCA data and TMM/RPKM values for chrom-seek data.
#          Uses code from DiffBind version 2. 
# 
# For use on biowulf:
#    sinteractive --mem=50G --gres=lscratch:200 --cpus-per-task=4
#    These are estimates and do not work for all projects. If you have errors, the most common issues beyond not setting up
#       the singularity correctly tend to involve not giving it enough memory.
# 
#  This is to load DiffBind v2. Make sure you list all directories you will need access to directly after the -B flag.
#    module load singularity.
#    singularity exec -B /data/NHLBI_IDSS,/data/OpenOmics,$PWD /data/OpenOmics/SIFs/cfchip_toolkit_v0.5.0.sif /bin/bash
#  
#  To use this function, open R within the singularity, source this script, and run one of the two functions within.
#
# Functions:
#    PCApeaks(csvfile, outroot)
#                To get PCA data using only peak positions
#
#    PCAcounts(csvfile, outroot, counts=TRUE, macsNarrow=FALSE)
#                To get PCA on consensus peaks with TMM-normalized
#                counts with ability to also save per sample
#                TMM- and RPKM-normalized counts
#
# Variables:
#    csvfile: the name of the csvfile that DiffBind uses to the load all the data
#    outroot: a string to add to the start of all output file names
#    counts:  whether or not to make files of TMM- and RPKM-normalized
#             data for all consensus peaks across all samples
#    macsNarrow: if the peak caller used was macsNarrow, set this variable to TRUE
# 
####################
# Want to plot the PCA data? Here's an example using ggplot.
# 
# library(ggplot2)
# plotData <- read.delim("S16S3_vs_Swk26_PCA_peaks_only.txt")
#
# c1p <- as.numeric(gsub("PC1 %: ","",rownames(plotData)[nrow(plotData)-2]))
# c2p <- as.numeric(gsub("PC2 %: ","",rownames(plotData)[nrow(plotData)-1]))
#
# plotData <- plotData[grep("%",rownames(plotData),invert=T),]
#
# p <- ggplot(data=plotData,aes(x=PC1,y=PC2))
# p + geom_point() +
#     xlab(sprintf('PC1 [%2.0f%%]',c1p)) +
#     ylab(sprintf('PC2 [%2.0f%%]',c2p))
#
####################

suppressMessages(library(DiffBind))
library(parallel)

pv.pcmask <- function(pv,numSites, mask, sites,removeComps,cor=F,bLog=T){
   
   if(missing(numSites)) numSites <- nrow(pv$binding)
   if(is.null(numSites)) numSites <- nrow(pv$binding)  
   numSites <- min(numSites,nrow(pv$binding))
   
   if(missing(sites)) sites <- 1:numSites
   if(is.null(sites)) sites <- 1:numSites
   
   if(missing(mask)) mask <- rep(T,ncol(pv$class))
   for(i in which(mask)) {
      if(nrow(pv$peaks[[i]])==0) {
         mask[i]=F
      }
   }
   if(sum(mask)<2) {
      stop('Need at least two samples for PCA.')
   }
   
   res <- NULL   
   res$class <- pv$class
   pv$values <- pv$binding[sites,c(F,F,F,mask)]
   active   <- apply(pv$values,1,pv.activefun)
   numSites <- min(numSites,sum(active))
   
   pv$values <- pv$values[active,][1:numSites,]
   
   if(!missing(removeComps)) {
      pv$values <- pv.removeComp(pv$values,numRemove=removeComps)
   }
   
   if(bLog) {
      if(max(pv$values)>1) {
         pv$values[pv$values<=1] <- 1
         pv$values <- log2(pv$values)
      }
   }
   
   if(nrow(pv$values) >= sum(mask)) {
      res$pc <- prcomp(pv$values) #,cor=cor)
   }
   res$mask <- mask
   
   return(res)
}

pv.activefun <- function(x){
   if(sum(x>0)>0){
      return(TRUE)
   } else {
      return(FALSE)
   }	
}

PCAbasic <- function(pv, outfile) {
  numSites <- nrow(pv$binding)
  mask <- rep(T,ncol(pv$class))
  pv <- pv.pcmask(pv, numSites, mask, sites=NULL, cor=F, bLog=T)

  pc <- pv$pc

  vr <- rep(0,length(pc$sdev))
  for (i in 1:length(vr)) {
      vr[i] <- pc$sdev[i] ^ 2
  }

  comps=1:3
  c1 <- comps[1]
  c2 <- comps[2]
  c3 <- comps[3]

  c1p <- vr[c1] / sum(vr) * 100
  c2p <- vr[c2] / sum(vr) * 100
  c3p <- vr[c3] / sum(vr) * 100

  plotData <- as.data.frame(pc$rotation[,c(c1,c2,c3)])
  colnames(plotData) <- c("PC1","PC2","PC3")
  write.table(plotData, outfile,quote=F,append=T,sep="\t")
  write.table(paste0("PC1 %: ", c1p), outfile,quote=F,row.names=F,col.names=F,append=T)
  write.table(paste0("PC2 %: ", c2p), outfile,quote=F,row.names=F,col.names=F,append=T)
  write.table(paste0("PC3 %: ", c3p), outfile,quote=F,row.names=F,col.names=F,append=T)  
}


PCApeaks <- function(csvfile, outroot) {
  samples <- dba(sampleSheet=csvfile)
  pv <- samples
  outfile <- paste0(outroot, "_PCA_peaks_only.txt")
  PCAbasic(pv,outfile)
}

PCAcounts <- function(csvfile, outroot, counts=TRUE, macsNarrow=FALSE) {
  samples <- dba(sampleSheet=csvfile)
#  if ( grepl("narrow",samples$samples$Peaks[1]) ) {
  if (macsNarrow == TRUE) {
    DBdataCounts <- dba.count(samples, summits=250)
  } else {
    DBdataCounts <- dba.count(samples)
  }
  pv <- DBdataCounts
  outfile <- paste0(outroot, "_PCA_peaks_and_counts.txt")
  PCAbasic(pv,outfile)
  
  if (counts == TRUE) {
    # first TMM
    countsGR <- dba.peakset(DBdataCounts, bRetrieve=TRUE)
    outfile <- paste0(outroot, "_TMM.txt")
    write.table(as.data.frame(countsGR), outfile, sep="\t", quote=F)
    
    # then RPKM
    DBdataCounts <- dba.count(DBdataCounts,peaks=NULL,score=DBA_SCORE_RPKM_FOLD)
    countsGR <- dba.peakset(DBdataCounts, bRetrieve=TRUE)
    outfile <- paste0(outroot, "_RPKM.txt")
    write.table(as.data.frame(countsGR), outfile, sep="\t", quote=F) 
  }
}
