## FRIP_plot.R
## Created by Tovah Markowitz
## June 19, 2020
## Updated: Jan 19, 2022
## Updated: Novemeber 4, 2022

args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]

library(ggplot2)
library(rjson)

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
  p <- ggplot(inData2,aes(x=bamsample, y=FRiP, fill=groupInfo))
  p <- p + geom_bar(position="dodge",stat = "identity") +
    facet_wrap(.~bedtool) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = -15, hjust = 0)) +
    labs(title="All Samples",x="bam file", y ="Fraction of Reads in Peaks (FRiP)", 
         fill ="Group")
  pdf(paste0(folder, "/PeakQC/FRiP_barplot.pdf"))
  print(p)
  dev.off()
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

process_json <- function(injson) {
# to get the identities of the groups and the list of samples (ChIP and input)
# associated with it
  json  <- fromJSON(file = injson)
  groupsInfo <- json$project$groups
  inputs <- as.data.frame(json$project$peaks$inputs)
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
groupList <- process_json(paste0(folder,"/config.json"))

for (i in 1:length(groupList)) {
  group <- groupList[[i]]
  groupName <- names(groupList)[i]
  inData <- allData[which((allData$bedsample %in% group) & 
  	    	          (allData$bamsample %in% group)),]
  plot_barplots(inData, groupName, folder)
  plot_scatterplots(inData, groupName, folder)
}

selfData <- allData[which(allData$bedsample == allData$bamsample),]
groupInfo <- reshape2::melt(groupList)
names(groupInfo) <- c("bamsample","groupInfo")
selfData2 <- merge(selfData,groupInfo)
plot_barplots_self(selfData2, folder)
plot_scatterplots_self(selfData2, folder)
