#!/usr/bin/env python3

"""
Name: frip.py
Created by: Tovah Markowitz
Date: 06/18/20

Purpose: To calculate FRiP scores, one bam file and as many bedfiles as wanted as inputs
Currently only works with python/3.5
"""

##########################################
# Modules
import optparse
from pybedtools import BedTool
import pysam
import pandas as pd

##########################################
# Functions

def split_infiles(infiles):
    """ 
    breaks the infile string with space-delimited file names and 
    creates a list
    """
    infileList = infiles.strip("\'").strip('\"').split(" ")
    if len(infileList) == 1:
        infileList = infileList[0].split(";")
    return(infileList)

def count_reads_in_bed(bam, bedfile, genomefile):
    """
    some of this comes directly from the pybedtools site; read in 
    bed (or bed-like) file, sort it, and then count the number of 
    reads within the regions
    """
    bedinfo = BedTool(bedfile)
    bedinfo.sort(g=genomefile)
    return (
        BedTool(bam).intersect( bedinfo, bed=True, stream=True, )
        ).count()

def count_reads_in_bam(bam):
    """ count the number of reads in a given bam file """
    return( pysam.AlignmentFile(bam).mapped )

def calculate_frip(nreads, noverlaps):
    """ calculate FRiP score from nreads and noverlaps """
    return( float(noverlaps) / nreads )

def measure_bedfile_coverage(bedfile, genomefile):
    """ calculate the number of bases covered by a given bed file """
    bedinfo = BedTool(bedfile)
    return( bedinfo.sort(g=genomefile).total_coverage() )

def clip_bamfile_name(bamfile):
    """ 
    clip bam file name for table/plotting purposes; assumes file 
    naming system matches that of Pipeliner
    """
    sample = bamfile.split("/")[-1].split(".")[0]
    condition =  ".".join(bamfile.split("/")[-1].split(".")[1:-1])
    return( sample, condition )

def clip_bedfile_name(bedfile,filetype):
    """
    clip bed file name for table/plotting purposes; assumes file 
    naming system matches that of Pipeliner
    """
    if filetype == "":
        toolused = bedfile.split("/")[-3]
        sample = bedfile.split("/")[-2]
    else:
        toolused = filetype
        sample = bedfile.split("/")[-1].split(".")[0].strip("_peaks").strip("_broadpeaks")
    return( toolused, sample )

def process_files(bamfile, bedfiles, genome, filetypes):
    """ 
    this is the main function to take in list of input files and 
    put out an array containing key file name information, read 
    counts, and FRiP scores
    """
    bedfileL = split_infiles(bedfiles)
    filetypesL = split_infiles(filetypes)
    out = [[ "bedtool", "bedsample", "bamsample", "bamcondition", 
    "n_reads", "n_overlap_reads", "FRiP", "n_basesM" ]]
    nreads = count_reads_in_bam(bamfile)
    (bamsample, condition) = clip_bamfile_name(bamfile)
    for i in range(len(bedfileL)):
        bed = bedfileL[i]
        if len(filetypesL) > 1:
            filetype = filetypesL[i]
        else:
            filetype = filetypesL[0]
        (bedtool, bedsample) = clip_bedfile_name(bed,filetype)
        noverlaps = count_reads_in_bed(bamfile, bed, genome)
        frip = calculate_frip(nreads, noverlaps)
        nbases = measure_bedfile_coverage(bed, genome) / 1000000
        out.append( [bedtool, bedsample, bamsample, condition, 
                nreads, noverlaps, frip, nbases] )
    out2 = pd.DataFrame(out[1:], columns=out[0])
    return(out2)

def create_outfile_name(bamfile, outroot):
    """ uses outroot to create the output file name """
    (bamsample, condition) = clip_bamfile_name(bamfile)
    outtable = bamsample + "." + condition + "." + "FRiP_table.txt"
    if outroot != "":
        outtable = outroot + "." + outtable
    return(outtable)

def write_table(out2, outtable):
    out2.to_csv(outtable,sep='\t',index=False)


###############################################
# Main

def main():
    desc="""
This function takes a space-delimited or semi-colon delimited list
of bed-like files (extensions must be recognizable by bedtools)
and a single bam file. It will then calculate the FRiP score for
all possible combinations of files and save the information in a
txt file. It will also calculate the number of bases covered by 
each bed-like file. Note: this function assumes that the file 
naming system of the input files matches that of Pipeliner.
    """

    parser = optparse.OptionParser(description=desc)

    parser.add_option('-p', dest='peakfiles', default='', 
           help='A space- or semicolon-delimited list of peakfiles \
(or bed-like files).')
    parser.add_option('-b', dest='bamfile', default='', 
           help='The name of a bamfile to analyze.')
    parser.add_option('-g', dest='genomefile', default='', 
           help='The name of the .genome file so bedtools knows the \
size of every chromosome.')
    parser.add_option('-o', dest='outroot', default='', 
           help='The root name of the multiple output files. Default: ""')
    parser.add_option('-t', dest='filetypes', default='', help='A space- \
or semicolon-delimited list of input file sources/types. Only needed when \
source of bed file is not built into the script. Default: ""')

    (options,args) = parser.parse_args()
    bedfiles = options.peakfiles
    bamfile = options.bamfile
    genomefile = options.genomefile
    outroot = options.outroot
    filetypes = options.filetypes

    out2 = process_files(bamfile, bedfiles, genomefile, filetypes)
    outtable = create_outfile_name(bamfile, outroot)
    write_table(out2, outtable)

if __name__ == '__main__':
    main()

###############################################
# example cases

#bedfiles = "macs_broad/mWT_HCF1_mm_i81/mWT_HCF1_mm_i81_peaks.broadPeak macs_broad/mWT_HCF1_mm_i89/mWT_HCF1_mm_i89_peaks.broadPeak"
#bamfiles = "bam/Input_mm_i95.sorted.Q5DD.bam bam/mWT_HCF1_mm_i81.sorted.Q5DD.bam bam/mWT_HCF1_mm_i89.sorted.Q5DD.bam"
#genomefile = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/indexes/mm10.fa.sizes"
#out2 = pd.read_csv("FRIP_test.txt",sep="\t")
