# Example usage: Rscript runDiffBind.R 'directory' 'outfilename.html' 'input.csv' 'WT_vs_KO' 'macs_narrow' 
args <- commandArgs(trailingOnly = TRUE)

DIR <- args[1]
setwd(DIR)
outHtml <- args[2]

Sys.setenv(RSTUDIO_PANDOC="/usr/local/apps/rstudio/rstudio-1.1.447/bin/pandoc/")

rmarkdown::render("DiffBind_v2_ChIPseq.Rmd",output_file=outHtml, params = list(
	csvfile = args[3],
	contrasts = args[4],
	peakcaller = args[5]
))