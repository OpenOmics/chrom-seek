## FRIP_plot.R
## Created by Tovah Markowitz
## June 19, 2020
## Most recent update: October 17, 2024

args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]

library(ggplot2)
library(rjson)
# these are used but called directly by the function
#library(ComplexHeatmap)
#library(tidyr)
#library(circlize)

merge_files <- function(folder) {
  files <- list.files(path=paste0(folder,"/PeakQC"), pattern="FRiP_table.txt", 
  	   	full.names=T)
  allList <- lapply(files,read.table,header=T)
  allData <- do.call(rbind.data.frame, allList)
  write.table(allData, paste0(folder, "/PeakQC/FRiP_All_table.txt"), quote=F, 
  		row.names=F, sep="\t")
  return(allData)
}

plot_barplots <- function(inData, groupName, folder) {
  p <- ggplot(inData,aes(x=bamsample, y=FRiP, fill=bedsample))
  p <- p + geom_bar(position="dodge",stat = "identity") +
    facet_wrap(.~bedtool) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = -15, hjust = 0)) +
    labs(title=groupName, x="bam file", y ="Fraction of Reads in Peaks (FRiP)", 
         fill ="peak file")
  pdf(paste0(folder, "/PeakQC/", groupName,".FRiP_barplot.pdf"))
  print(p)
  dev.off()
}

plot_scatterplots <- function(inData, groupName, folder) {
  p <- ggplot(inData,aes(x=n_basesM, y=FRiP, shape=bedsample, color=bedtool))
  p <- p + geom_point(size=2.5) +
    facet_wrap(.~bamsample) +
    theme_bw() + 
    scale_x_continuous(trans = "log10") +
    labs(title=groupName, x="Number of Bases in Peaks (M)", 
         y="Fraction of Reads in Peaks (FRiP)",
         shape="peak file", color="peak calling tool")
  q <- p + annotation_logticks(sides="b")
  pdf(paste0(folder, "/PeakQC/", groupName,".FRiP_scatterplot.pdf"))
  tryCatch(print(q), error = function(e) {print(p)})
  dev.off()
}

plot_barplots_self <- function(inData2, folder) {
  bedtools <- unique(inData2$bedtool)
  for (i in 1:length(bedtools)) {
    inData3 <- inData2[which(inData2$bedtool == bedtools[i]),]
    p <- ggplot(inData3,aes(x=bamsample, y=FRiP, fill=groupInfo))
      p <- p + geom_bar(position="dodge",stat = "identity") +
        theme_bw() +
        theme(axis.text.x=element_text(angle = -15, hjust = 0)) +
        labs(title="All Samples",x="bam file", y ="Fraction of Reads in Peaks (FRiP)", 
           fill ="Group")
    pdf(paste0(folder, "/PeakQC/",bedtools[i],".FRiP_barplot.pdf"))
    print(p)
    dev.off()
  }
}

plot_scatterplots_self <- function(inData2, folder) {
  p <- ggplot(inData2,aes(x=n_basesM, y=FRiP, shape=bedtool, color=groupInfo))
  p <- p + geom_point(size=2.5) +
    theme_bw() + 
    scale_x_continuous(trans = "log10") +
    annotation_logticks(sides="b") +
    labs(title="All samples", x="Number of Bases in Peaks (M)", 
         y="Fraction of Reads in Peaks (FRiP)",
         shape="peak file", color="peak calling tool")
  pdf(paste0(folder, "/PeakQC/FRiP_scatterplot.pdf"))
  print(p)
  dev.off()
}

PeakQCHeatmap<-function(PeakQC.dir, peakcaller, plot_data){
# From Subrata Paul and adapted for this script
#  files<-grep(peakcaller, list.files(PeakQC.dir, full.names = T), value = T)
#  plot_data = lapply(files, function(file) read.table(file, sep = '\t', header = T))
#  plot_data = do.call('rbind', plot_data)
  plot_data = plot_data[which(plot_data$bedtool == peakcaller),]
  plot_data = plot_data[,c('bedsample', 'bamsample', 'FRiP')]
  plot_data = tidyr::pivot_wider(plot_data,names_from = bamsample, values_from=FRiP)
  plot_data = data.frame(plot_data, check.names = F)
  rownames(plot_data)<-plot_data$bedsample
  plot_data = as.matrix(plot_data[, -1])

  ComplexHeatmap::Heatmap(plot_data, na_col="grey",
                          col = circlize::colorRamp2(c(0, 0.1, max(plot_data)), c('red', 'orange', 'blue')), 
                          heatmap_legend_param = list(title = 'FRiP'), 
                          row_title = 'Sample of peaks', 
                          column_title = 'Sample of reads')  
}

process_json <- function(injson) {
# to get the identities of the groups and the list of samples (ChIP and input)
# associated with it
  json  <- fromJSON(file = injson)
  groupsInfo <- json$project$groups
  inputs <- unlist(json$project$peaks$inputs)
  for (i in 1:length(groupsInfo)) {
    tmp <- unique(unlist(inputs[names(inputs) %in% groupsInfo[[i]]]))
    if (length(tmp) > 1) {
       groupsInfo[[i]] <- c(groupsInfo[[i]],as.character(tmp))
    } else if (tmp != "" ) {
       groupsInfo[[i]] <- c(groupsInfo[[i]],as.character(tmp))
    }
  }
  return(groupsInfo)
}

allData <- merge_files(folder)

bedtools <- unique(allData$bedtool)

for (i in 1:length(bedtools) ) {
    ht <- PeakQCHeatmap(PeakQC.dir=folder, peakcaller=bedtools[i], plot_data=allData)
    pdf(paste0(folder, "/PeakQC/",bedtools[i],".FRiP_heatmap.pdf"))
    print(ht)
    dev.off()
}

groupList <- process_json(paste0(folder,"/config.json"))

for (i in 1:length(groupList)) {
  group <- groupList[[i]]
  groupName <- names(groupList)[i]
  inData <- allData[which((allData$bedsample %in% group) & 
  	    	          (allData$bamsample %in% group)),]
  plot_scatterplots(inData, groupName, folder)
}

selfData <- allData[which(allData$bedsample == allData$bamsample),]
groupInfo <- reshape2::melt(groupList)
names(groupInfo) <- c("bamsample","groupInfo")
selfData2 <- merge(selfData,groupInfo)
plot_barplots_self(selfData2, folder)
write.table(selfData2, paste0(folder, "/PeakQC/FRiP_table.txt"), quote=F,
                row.names=F, sep="\t")
