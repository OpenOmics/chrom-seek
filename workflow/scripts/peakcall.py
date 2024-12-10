#!/usr/bin/env python3
from os.path import join

def get_input_bam(input_sample, bam_dir):
    """
    Returns a ChIP samples input BAM file,
    see chip2input for ChIP, Input pairs.
    """
    if input_sample:
        # Runs in a ChIP, input mode
        return join(bam_dir, "{0}.Q5DD.bam".format(input_sample))
    # Runs in ChIP-only mode
    return []


def get_control_input(ext, paired_end, bam_dir):
    i = []
    if paired_end and ext != "":
        i = [join(bam_dir, "{0}.Q5DD.bam".format(ext))]
    elif not paired_end and ext != "":
        i = [join(bam_dir, "{0}.Q5DD_tagAlign.gz".format(ext))]
    return i


def outputIDR(groupswreps, groupdata, chip2input, tools):
    """
    Produces the correct output files for IDR. All supposed replicates
    should be directly compared when possible using IDR. IDR malfunctions
    with bed files and GEM so it will not run with either of those.
    Because there is no q-value calculated for SICER when there is no 
    input file, those samples are also ignored.
    """
    IDRgroup, IDRsample1, IDRsample2, IDRpeaktool = [], [], [], []
    for group in groupswreps:
        nsamples = len(groupdata[group])
        for i in range(nsamples):
            ctrlTF = chip2input[groupdata[group][i]] != ""
            for j in range(i+1,nsamples):
                if ctrlTF == (chip2input[groupdata[group][j]] != ""):
                    if ctrlTF == False:
                        tooltmp = [ tool for tool in tools if tool != "sicer" ]
                    else:
                        tooltmp = tools			           
                    IDRgroup.extend([group] * len(tooltmp))
                    IDRsample1.extend([groupdata[group][i]] * len(tooltmp))
                    IDRsample2.extend([groupdata[group][j]] * len(tooltmp))
                    IDRpeaktool.extend(tooltmp)
    return( IDRgroup, IDRsample1, IDRsample2, IDRpeaktool )


def zip_peak_files(chips, PeakTools, PeakExtensions):
    """Making input file names for FRiP"""
    zipSample, zipTool, zipExt = [], [], []
    for chip in chips:
        for PeakTool in PeakTools:
            zipSample.append(chip)
            zipTool.append(PeakTool)
            zipExt.append(PeakExtensions[PeakTool])
    return(zipSample, zipTool, zipExt)


def calc_effective_genome_fraction(effectivesize, genomefile):
    """
    calculate the effective genome fraction by calculating the
    actual genome size from a .genome-like file and then dividing
    the effective genome size by that number
    """
    lines=list(map(lambda x:x.strip().split("\t"),open(genomefile).readlines()))
    genomelen=0
    for chrom,l in lines:
        if not "_" in chrom and chrom!="chrX" and chrom!="chrM" and chrom!="chrY":
            genomelen+=int(l)
    return(str(float(effectivesize)/ genomelen))


def getSicerChips(bam_dir, name, paired_end):
    if paired_end:
        chip = join(bam_dir, name + ".Q5DD.bam")
    else:
        chip = join(bam_dir, name + ".Q5DD_tagAlign.gz")
    return chip


def getSicerFragLen(ppqt_dir, qc_dir, name, paired_end):
    if paired_end:
        fragLen = join(qc_dir, name + ".Q5DD.insert_size_metrics.txt")
    else:
        fragLen = join(ppqt_dir, name + ".Q5DD_tagAlign.ppqt.txt")
    return fragLen


def get_manorm_sizes(g1, g2, group_data, ppqt_in):
    if not ppqt_in:
        return ""
    file = lambda w, _in: list(map(lambda z: z.strip().split(), open(ppqt_in, 'r').readlines()))
    extsize1 = [ppqt[1] for ppqt in file if ppqt[0] == group_data[g1]][0]
    extsize2 = [ppqt[1] for ppqt in file if ppqt[0] == group_data[g2]][0]
    return f"--s1 {extsize1} --s2 {extsize2}"