####################
#
# Name: promoterAnnotationByGene.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: August 8, 2022
# Updated: October 26, 2022 to work with uropa 4.0.2
# Updated: November 3, 2022 to fit with pipeline
# 
####################
#
# Purpose: To take UROPA allhits output files using "TSSprot" conditions and
#          create a table of which genes have annotations overlapping their 
#          promoters and how many times. Output format: dataframe
#
# Details: Promoters will be defined as 3kb upstream to 1 kb downstream of the 
#          TSS. Allhits files were chosen to capture information from "peaks" 
#          overlappingmultiple promoters. Finalhits files can also be processed 
#          with this pipeline. This script can handle multiple allhits files as 
#          long as there are equal numbers of sampleNames to go with them. Also, 
#          giving a matching DiffBind txt file will allow the allhits file to be 
#          filtered to only include the significant differential peaks or to 
#          split the data by the direction of log fold-change.
#
# Requires: GenomicRanges, tidyr
#
# Function: promoterAnnotationByGene(allhitsFiles, sampleNames, diffbindFiles=NA, direction=NA)
# 
# Variables:
#     allhitsFiles:  [required] a vector of allhits files to process
#     sampleNames:   [required] a vector of short names for each allhits file 
#                               to use as column headers
#     diffbindFiles: [optional] a vector of diffbind files to use to filter each 
#                               allhits file
#     direction:     [optional] when filtering using diffbindFiles, define how to
#                               filter using log fold change. "Both" is default
#                               when not defined by user.
#                               Options: "both", "pos", "neg", "separate"
#
# Example usage:
#     source("promoterAnnotation_by_Gene.R")
#     out1 <- promoterAnnotationByGene(allhitsA.txt, "A")
#     out2 <- promoterAnnotationByGene(allhitsA.txt, "A", diffbindA.txt, "both")
#     out3 <- promoterAnnotationByGene(allhitsFiles= c(allhitsA.txt, allhitsB.txt), 
#                                      sampleNames=c("A","B"), 
#                                      diffbindFiles=c(diffbindA.txt,diffbindB.txt), 
#                                      direction="pos")
#     out4 <- promoterAnnotationByGene(allhitsFiles= c(allhitsA.txt, allhitsA.txt), 
#                                      sampleNames=c("Deseq2","EdgeR"), 
#                                      diffbindFiles=c(Deseq2.txt,EdgeR.txt), 
#                                      direction="separate")
# 
####################


allhits2promoter <- function(allhitsFile) {
  # cleaning up the allhits file to only keep information about peaks 
  # overlapping promoters
  inData <- read.delim(allhitsFile)
  tmp <- which(inData$name == "query_1")
  if (length(tmp) == 0) {
    print (paste0("Supplied file ", allhitsFile, " has no peaks overlapping promoters."))
  } else {
    promoterData <- inData[tmp,]
    promoterData <- promoterData[,c("peak_chr", "peak_start", "peak_end", "gene_id", "gene_name")]
    return(promoterData)
  }
}

filterPromoter <- function(Diffbind, promoterData, sampleName) {
  # used by DiffbindFilterPromoter
  promoterData2 <- GenomicRanges::makeGRangesFromDataFrame(promoterData, seqnames.field="peak_chr",
                                  start.field="peak_start", end.field="peak_end",
                                  starts.in.df.are.0based=F)
  Diffbind2 <- GenomicRanges::makeGRangesFromDataFrame(Diffbind)
  ov <- GenomicRanges::countOverlaps(promoterData2,Diffbind2,type = "equal",maxgap=1)
  promoterData3 <- promoterData[which(ov != 0),]
  promoterData3$sample_id <- sampleName
  return(promoterData3)
}

DiffbindFilterPromoter <- function(DiffbindFile, promoterData, sampleName, direction) {
  # filters the promoter data based upon whether it matches a different peak and what direction the fold-change is
  # direction can be: "both", "pos", "neg", "separate". If direction is NA, use "both".
    Diffbind <- read.delim(DiffbindFile)
    Diffbind <- Diffbind[which(Diffbind$FDR < 0.05),]
    if ((direction == "both") | is.na(direction)) {
      promoterData2 <- filterPromoter(Diffbind, promoterData, sampleName)
    } else if (direction == "pos") {
      sampleName <- paste0(sampleName, "_pos")
      Diffbind <- Diffbind[which(Diffbind$Fold > 0),]
      promoterData2 <- filterPromoter(Diffbind, promoterData, sampleName)
    } else if (direction == "neg") {
      sampleName <- paste0(sampleName, "_neg")
      Diffbind <- Diffbind[which(Diffbind$Fold < 0),]
      promoterData2 <- filterPromoter(Diffbind, promoterData, sampleName)
    } else {
      sampleNameP <- paste0(sampleName, "_pos")
      DiffbindP <- Diffbind[which(Diffbind$Fold > 0),]
      promoterDataP <- filterPromoter(DiffbindP, promoterData, sampleNameP)
      sampleNameN <- paste0(sampleName, "_neg")
      DiffbindN <- Diffbind[which(Diffbind$Fold < 0),]
      promoterDataN <- filterPromoter(DiffbindN, promoterData, sampleNameN)
      promoterData2 <- rbind(promoterDataP, promoterDataN)
    }
return(promoterData2)
}

createPromoterTable <- function(promoterData) {
  # making final output table
  PromoterTable <- data.frame( table(promoterData[,c("gene_id", "sample_id")] ) )
  PromoterTable2 <- merge( unique(promoterData[,c("gene_id", "gene_name")] ), PromoterTable)
  PromoterTable3 <- tidyr::pivot_wider(PromoterTable2, names_from="sample_id", values_from="Freq")
  return(PromoterTable3)
}

promoterAnnotationByGene <- function(allhitsFiles, sampleNames, diffbindFiles=NA, direction=NA) {
  # the main function
  if ( length(allhitsFiles) != length(sampleNames) ) {
    print("Number of allhits files and sample names don't match.")
  } else {
    if ( (length(allhitsFiles) != length(diffbindFiles)) & (sum(is.na(diffbindFiles)) != 1) ) {
      print("Number of allhits files and diffbind files don't match.")
    } else {
      if ( length(allhitsFiles) == 1 ) {
        promoterData <- allhits2promoter(allhitsFiles)
        if (is.na(diffbindFiles)) {
          promoterData$sample_id <- sampleNames
        } else {
          promoterData <- DiffbindFilterPromoter(diffbindFiles, promoterData, sampleNames, direction)
        } 
      } else {
        for ( a in 1:length(allhitsFiles) ) {
          print(a)
          tmpA <- allhits2promoter(allhitsFiles[a])
          if (sum(is.na(diffbindFiles)) ==1) {
            tmpA$sample_id <- sampleNames[a]
          } else {
            tmpA <- DiffbindFilterPromoter(diffbindFiles[a], tmpA, sampleNames[a], direction)
          }
          if (a == 1) {
            promoterData <- tmpA
          } else {
            promoterData <- rbind(promoterData, tmpA)
          }
        }
      }
    }
    promoterTable <- createPromoterTable(promoterData)
    return(promoterTable)
  }   
}     

peakcallVersion <- function(inFolder,outFile) {
# currently only works for macs outputs
# inFolder here is the folder where the uropa output files are located
  filesA <- list.files(path=inFolder,pattern="allhits.txt")
  samples <- matrix(unlist(strsplit(filesA,"_macs")),ncol=2,byrow=T)[,1]
  filesA <- list.files(path=inFolder,pattern="allhits.txt",full.names = T)
  promoterInfo <- promoterAnnotationByGene(allhitsFiles=filesA, sampleNames=samples)
  write.table(promoterInfo, outFile, quote=F,sep="\t",row.names=F)
}

diffbindVersion <- function(inFolder,outFile) {
# currently designed for macs peaks, analyzed by deseq2
# analyzing both positive and negative together for now
# inFolder here is the root working directory for the project
  uropaFolder <- paste0(inFolder, "/UROPA_annotations/DiffBind")
  diffbindFolder <- paste0(inFolder, "/DiffBind")
  filesU <- list.files(path=uropaFolder, pattern="DiffbindDeseq2_uropa_protTSS_allhits.txt")
  samples <- matrix(unlist(strsplit(filesU,"-macs")),ncol=2,byrow=T)[,1]
  filesU <- list.files(path=uropaFolder, pattern="DiffbindDeseq2_uropa_protTSS_allhits.txt",full.names=T)
  filesD <- list.files(path=diffbindFolder, pattern="Deseq2.txt",full.names=T,recursive=T)
  promoterInfo <- promoterAnnotationByGene(allhitsFiles=filesU, 
                     sampleNames=samples, diffbindFiles=filesD, direction="both")
  write.table(promoterInfo, outFile, quote=F,sep="\t",row.names=F)
}
