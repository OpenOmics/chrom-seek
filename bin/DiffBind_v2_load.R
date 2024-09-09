#!/usr/bin/env Rscript
library(argparse)
library(stringi)

cleanup_arg <- function(arg) {
    arg <- stri_replace_all_fixed(arg, '"', '')                    # remote double quotes
    arg <- stri_replace_all_fixed(arg, "'", "")                    # remove single quotes
    arg <- stri_replace_all_charclass(arg, "\\p{WHITE_SPACE}", "") # remove all whitespace
    return(arg)
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
DBdataCounts <- dba.count(samples, summits=summits_arg, bParallel=T)

# save counts
saveRDS(DBdataCounts, counts)

# save peaklist
consensus <- dba.peakset(DBdataCounts, bRetrieve=T)
consensus$name <- paste0("Peak", 1:length(consensus))
rtracklayer::export(consensus, list)
