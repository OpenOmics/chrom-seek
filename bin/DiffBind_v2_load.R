#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description= 'This progrom does stuff')
parser$add_argument('--csvfile', '-c', help= 'I am the input file')
parser$add_argument('--contrasts', '-r', help= 'I am the output file')
parser$add_argument('--peakcaller', '-p', help= 'Some filtering cutoff')
parser$add_argument('--threads', '-t', type="integer", help= 'Some filtering cutoff')
xargs <- parser$parse_args()

csvfile <- xargs$csvfile
contrasts <- xargs$contrasts
contrasts <- gsub('"', "", contrasts)
contrasts <- gsub("'", "", contrasts)
peakcaller <- xargs$peakcaller
peakcaller <- gsub('"', "", peakcaller)
peakcaller <- gsub("'", "", peakcaller)
threads <- xargs$threads
threads <- as.numeric(threads)

suppressMessages(library(DiffBind))
suppressMessages(library(parallel))

samples <- dba(sampleSheet=csvfile)

if ( peakcaller == "macsNarrow" ) {
  summits_arg <- 250
  print ("Ran macsNarrow.")
  print ("Differential peaks are 250bp upstream and downstream of the summits.")
} else {
  summits_arg <- FALSE
  print ("Assuming broad peak calling tool.")
  print ("Differential peaks are consensus peaks.")
}

# count
DBdataCounts <- dba.count(samples, summits=summits_arg, bParallel=threads)

# save counts
counts <- paste0(contrasts, "-", peakcaller, "_Diffbind_counts.RDS")
saveRDS(DBdataCounts, counts)

# save peaklist
peaklist <- paste0(contrasts, "-", peakcaller, "_Diffbind_fullList.bed")
consensus <- dba.peakset(DBdataCounts, bRetrieve=T)
consensus$name <- paste0("Peak", 1:length(consensus))
rtracklayer::export(consensus, peaklist)
