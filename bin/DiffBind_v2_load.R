#!/usr/bin/env Rscript
library(argparse)
library(stringi)
library(DiffBind)
library(parallel)

cleanup_arg <- function(arg) {
    arg <- stri_replace_all_fixed(arg, '"', '')                    # remote double quotes
    arg <- stri_replace_all_fixed(arg, "'", "")                    # remove single quotes
    arg <- stri_replace_all_charclass(arg, "\\p{WHITE_SPACE}", "") # remove all whitespace
    return(arg)
}

remove_negative_coordinates <- function(dbaOBJ) {
  # This version is to remove negative coordinates caused by the DiffBind summits calculation
  negCoordIdx <- unique(which(dbaOBJ$merged[,2] <= 0))
  dbaOBJ$merged[negCoordIdx,2] <- 1
  dbaOBJ$binding[negCoordIdx,2] <- 1
  for ( i in 1:length(dbaOBJ$peaks) ) {
      dbaOBJ$peaks[[i]]$Start[negCoordIdx] <- 1
  }
  return(dbaOBJ)
}

remove_negative_coordinates2 <- function(dbaOBJ) {
  # This version is to remove negative coordinates deriving from the peak callers
  negCoordIdx <- unique(which(dbaOBJ$merged[,2] <= 0))
  dbaOBJ$merged[negCoordIdx,2] <- 1
  return(dbaOBJ)
}


parser <- ArgumentParser(description= 'Load diffbind csv process with R::dba return RDS and BED')
parser$add_argument('--csvfile', '-c', help='CSV input file from `diffbind_csv`')
parser$add_argument('--counts', '-n', help='Peak count RDS output file', default=file.path(getwd(), "peak_counts.rds"))
parser$add_argument('--list', '-l', help='Peak list TXT output file', default=file.path(getwd(), "peak_list.bed"))
parser$add_argument('--peakcaller', '-p', help='String with peakcaller name')
xargs <- parser$parse_args()

csvfile <- cleanup_arg(xargs$csvfile)
counts <- cleanup_arg(xargs$counts)
list <- cleanup_arg(xargs$list)
threads <- as.numeric(cleanup_arg(xargs$threads))
peakcaller <- cleanup_arg(xargs$peakcaller)
samples <- dba(sampleSheet=csvfile)
samples <- remove_negative_coordinates2(samples)

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
DBdataCounts <- dba.count(samples, summits=summits_arg, bParallel=T)

# remove negative coordinates when summits_arg is not FALSE
if ( summits_arg != FALSE ) {
    DBdataCounts <- remove_negative_coordinates(DBdataCounts)
}

# save counts
saveRDS(DBdataCounts, counts)

# save peaklist
consensus <- dba.peakset(DBdataCounts, bRetrieve=T)
consensus$name <- paste0("Peak", 1:length(consensus))
rtracklayer::export(consensus, list)
