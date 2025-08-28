suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicAlignments))

#inFile <- "K4Me3_S16_HC3_200.bam"
#outFile <- "K4Me3_S16_HC3_200.tagAlign"

bam2fragment <- function(inFile, outFile) {
  print(paste0("loading: ",inFile))
  inData <- readGAlignmentPairs(inFile)
  print("converting to fragments")
  fragments <- granges(inData)
  rtracklayer::export(fragments, outFile, format="bed")
  print("done")
}

args<-commandArgs(TRUE)

bam2fragment(inFile=args[1], outFile=args[2])