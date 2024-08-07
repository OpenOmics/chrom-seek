####################
#
# Name: cfChIP_signatures.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: August 9, 2022
# 
####################
#
# Purpose: To take the individual cfChIP signature tables, combine them, 
#          and create the preferred output plot
#
# Functions: mergeSignatures and plotSignatures
#
# Requires: ggplot2 and ggprism (for plotting)
#
# Details: mergeSignatures will take a folder of signatures and combine them 
#         into one long table. plotSignatures can directly take the the output 
#         of mergeSignatures, but you can also load the data into R and filter
#         to only include a subset of samples before running the function. Also,
#         add a column called "Condition" either to the input txt file or the R
#         object to group columns in the plot using that additional information.
#
# Function1: mergeSignatures(folder, outFile)
#   folder:  [required] the path to the folder containing the individual signature
#            files direct from the cfChIP tool
#   outFile: [required] the name of the output txt file to save the data
#
# Function2: plotSignatures(inTXT, outPDF)
#   inTXT: [required] either the name of the file from mergeSignatures or an
#          an R object containing the data. Column names must match that of the
#          output of mergeSignatures, but column order doesn't matter
#   inPDF: [required] the name of the output pdf file to create
#
# Example usage:
#   source("cfChIP_signatures.R")
#   mergeSignatures("cfChIPtool/Output/H3K4me3/Signatures/", out.txt)
#   plotSignatures(out.txt, out.pdf)
#   plotSignatures(signatureDataFrame, out.pdf)
# 
####################

mergeSignatures <- function(folder, outFile) {
  files <- list.files(folder,full.names = T)
  sigList <- lapply(files, read.csv)
  samples <- gsub(".csv","",grep("csv",unlist(strsplit(files,"/")),value=T))
  for (i in 1:length(samples)) {
    sigList[[i]] <- data.frame(sigList[[i]],Sample=samples[i])
  }
  sigData <- do.call("rbind",sigList)
  write.table(sigData, outFile, quote=F, sep="\t", row.names=F)
}

plotSignatures <- function(inTXT, outPDF) {
  library(ggplot2)
  library(ggprism)

  if (mode(inTXT) == "character") { # if using a file name
    sigData <- read.delim(inTXT)
  } else { # if starting with an object in R
   sigData <- as.data.frame(inTXT) 
  }
  sigData$NormalizedCounts[which(sigData$NormalizedCounts > 3)] <- 3
  sigData$NormalizedCounts[which(sigData$NormalizedCounts < 0.15)] <- 0.15
  sigData$qValue[which(sigData$qValue > 300)] <- 300
  sigData$qValue[which(sigData$qValue < 5)] <- NA
  names(sigData)[1] <- "cellType"

  cellTypes <- data.frame(cellType=c("Neutrophils","Monocytes","Megakaryocyte",
                                   "Erythroblast","T-Cells","B-Cells","NK",
                                   "Vasculary","Adipose","Skin","Sk. Muscle",
                                   "Brain","Heart","Lung","Breast","Digestive",
                                   "Pancreas"),
                        class=c(rep("Blood",7),rep("Global",4),rep("Other",6)) )

  sigData2 <- merge(sigData,cellTypes)
  sigData2$cellType <- factor(sigData2$cellType,levels=rev(cellTypes$cellType))

  pdf(outPDF)
  p <- ggplot(data=sigData2,aes(x=Sample,y=cellType,color=NormalizedCounts,size=qValue))
  p <- p + geom_point() + 
    scale_size(limits=c(5,300),breaks=c(5,50,100,150,200,250,300),
               labels=paste0("e-",c(5,50,100,150,200,250,300))) +
    scale_color_viridis_c(direction = -1, option="A") +
    theme_prism() +
    theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45,vjust=0.9,hjust=1))
  if (sum(names(sigData) == "Condition") == 1) {
    p <- p + facet_grid(rows=vars(class),cols=vars(Condition),scales="free",space="free")
  } else {
    p <- p + facet_grid(rows=vars(class),scales="free",space="free")
  }
  print(p)
  dev.off()
}