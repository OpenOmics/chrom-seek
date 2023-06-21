####################
#
# Name: significantPathways.R
# Created by: Tovah Markowitz, PhD
# Bioinformatics (NCBR)/ Integrated Data Sciences Section (IDSS)
# Research Technologies Branch/DIR/NIAID
#
# Created: August 9, 2022
# Updated: October 28, 2022 to make reactomePA optional
# Updated: November 4, 2022 to accept a txt file or a gtf for the background genes
#                           also to accept a promoter annotation table and analyze every column
# 
####################
#
# Purpose: To take a list of genes and find the significant KEGG or Reactome
#          pathways using overenrichment analysis. See details for specialized functionality.
#
# Requires: clusterProfiler, ReactomePA, enrichplot, org.Hs.eg.db, rtracklayer, ggplot2, and ggprism
#
# Details: Takes input gene lists as Ensembl gene IDs or gene symbols, converts to
#          Entrez gene IDs, and runs ORA against KEGG or Reactome database. Requires
#          a background gene list as cfChIP currently ignores chrs X, Y, and M.
#          Outputs a dataframe of significant pathways, a pdf of the top most
#          significant pathways, or a pdf of just the pathways of interest (if significant).
#
# Function: significantPathways(Genes, bkgGeneTXT, database="KEGG", PDFfile=NA, pathwayVector=NA)
#
# Variables:
#   Genes:         [Required] a vector of the genes to be analyzed through ORA
#   bkgGeneFILE:   [Required] a txt file containing a column of Ensembl IDs listing 
#                  the appropriate background gene set or the gtf file used for the uropa
#                  annotations 
#                  For example: hg19.ensembl.prot_coding.with_annotations.txt
#   database:      [Optional] whether to compare to the KEGG or Reactome database
#                  default: KEGG
#   PDFfile:       [Optional] name of the PDF file to create, if empty no PDF will be made
#   pathwayVector: [Optional] a vector of pathways (descriptions or IDs) to plot in the pdf.
#                  If PDFfile is empty, it is ignored. If this is empty and PDFfile is not,
#                  pdf plot will be of the top 30 most significant pathways instead.
#
# Example usage:
#   source("significantPathways.R")
#   out <- significantPathways(Genes= c("GeneA","GeneB"), 
#                              bkgGeneFILE= "hg19.ensembl.prot_coding.with_annotations.txt",
#                              database="KEGG", PDFfile="a.pdf", 
#                              pathwayVector=c("pathwayA", "pathwayB"))
# 
####################

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggprism)

makeBarplotTop <- function(inData, titleName, PDFfile) {
  inDataCount <- sum(inData@result$p.adjust < 0.1)
  if (inDataCount > 30) { inDataCount = 30 } 
  if (inDataCount > 0) {
    pdf(PDFfile)
    print(barplot(inData, showCategory = inDataCount,
                  label_format=70, title=titleName, x="GeneRatio") + 
            theme_prism(base_size =8) + theme(legend.title = element_text()) )
  }
  dev.off()
}

makeBarPlotSelect <- function(inData, titleName, PDFfile, categories) {
  pdf(PDFfile)
  print(barplot(inData, showCategory = categories,
                label_format=70, title=titleName, x="GeneRatio") + 
          theme_prism(base_size = 8) + theme(legend.title = element_text()) )
  dev.off()
  }

processGenes <- function(geneIDs) {
  if (grepl("^ENSG", geneIDs[1])) {
    ensIDs <- gsub("\\.[0-9]+", "", geneIDs, perl=T)
    entrezIDs <- bitr(ensIDs, from= "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  } else {
    entrezIDs <- bitr(ensIDs, from= "SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  }
  entrezIDs <- entrezIDs$ENTREZID
  return(entrezIDs)
}

significantPathways <- function(Genes, bkgGeneFILE, database="KEGG", PDFfile=NA, pathwayVector=NA) {
  sigGenes <- processGenes(Genes)
  if (grepl("gtf",bkgGeneFILE)) {
    bkgGenesData <- rtracklayer::import(bkgGeneFILE)
    bkgGenes <- unique(bkgGenesData$gene_id)
  } else {
    bkgGenesData <- read.delim(bkgGeneFILE)
    bkgGenes <- bkgGenesData[,grep("^ENSG", bkgGenesData[1,])]
  }
  backgroundGenes <- processGenes(bkgGenes)
  if (database == "KEGG") {
    pathwaySig <- enrichKEGG(sigGenes, organism= "hsa", keyType="kegg", universe=backgroundGenes)
    pathwayData <- pathwaySig@result[which(pathwaySig@result$p.adjust < 0.1),]
  } else {
    library(ReactomePA)
    pathwaySig <- enrichPathway(sigGenes, readable=T, universe=backgroundGenes)
    pathwayData <- pathwaySig@result[which(pathwaySig@result$p.adjust < 0.1),]
  }
  if (!is.na(PDFfile)) {
    if(length(pathwayVector) != 1) {
      if (length(grep("HSA", pathwayVector, ignore.case=T)) != 0) {
        pathwayVector <- pathwayData$Description[which(pathwayData$ID %in% pathwayVector)]
      }
      makeBarPlotSelect(inData=pathwaySig, titleName=database, PDFfile=PDFfile, categories=pathwayVector)
    } else {
      makeBarplotTop(inData=pathwaySig, titleName=database, PDFfile=PDFfile)
    }
  }
  return(pathwayData)
}

promoterAnnotationWrapper <- function(promoterFile, bkgGeneFILE, database="KEGG") {
   promoterData <- read.delim(promoterFile)
   outFolder <- dirname(promoterFile)
   for (i in 3:ncol(promoterData)) {
      colName <- names(promoterData)[i]
      Genes <- promoterData$gene_id[which(promoterData[,i] > 0)]
      outData <- significantPathways(Genes, bkgGeneFILE, database)
      outFileName <- paste0(outFolder,"/",colName,"_",database,".txt")
      write.table(outData, outFileName, quote=F, row.names=F, sep="\t")
   }
}