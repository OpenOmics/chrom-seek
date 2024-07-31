#!/usr/bin/env python3
# ~~~ Common helper functions for grouping of outputs
from os.path import join

# common functions related to sample grouping or group meta-information
def group_samples_by_reps(groupdata, samples, chip2input):
    groupdatawinput = {}
    groupswreps = []
    for group, chipsamples in groupdata.items() :
        tmp = [ ]
        if len(chipsamples) > 1:
            groupswreps.append(group)
        for chip in chipsamples :
            if chip in samples:
                tmp.append(chip)
                input = chip2input[chip]
                if input != 'NA' and input != '':
                    tmp.append(input)
        if len(tmp) != 0:
            groupdatawinput[group]=set(tmp)
    return groupdatawinput, groupswreps


def group_output_files(extensions, groupslist, inputnorm):
    """
    Produces correct output filenames based on group information.
    Names will be:
    Inputnorm.Q5DD.RPGC.metagene_heatmap.pdf
    {groupName}.Q5DD.RPGC.metagene_heatmap.pdf
    {groupName}.sorted.RPGC.metagene_heatmap.pdf
    Note: Inputnorm will only be included when there are input samples.
    """
    dtoolgroups, dtoolext = [], []
    
    if len(inputnorm) == 2:
            dtoolgroups.extend(["InputNorm"])
            dtoolext.extend([extensions[1]])
    
    for group in groupslist:
            dtoolgroups.extend([group] * 2)
            dtoolext.extend([extensions[1], extensions[0]])
    
    if len(inputnorm) == 2:
            dtoolgroups.extend(["InputNorm.prot"])
            dtoolext.extend([extensions[1]])
    
    for group in groupslist:
            dtoolgroups.extend([group + ".prot"] * 2)
            dtoolext.extend([extensions[1], extensions[0]])
    
    return dtoolgroups, dtoolext


def zip_contrasts(contrast, PeakTools):
    """
        making output file names for differential binding analyses
    """
    zipGroup1, zipGroup2, zipTool, contrasts = [], [], [], []
    for g1, g2 in contrast:
        for PeakTool in PeakTools:
            zipGroup1.append(g1)
            zipGroup2.append(g2)
            zipTool.append(PeakTool)
            contrasts.append( g1 + "_vs_" + g2 + "-" + PeakTool )
    return(zipGroup1, zipGroup2, zipTool, contrasts)


def get_peaktools(assay_type):
    tools = ["macsNarrow"]
    if assay_type == "atac": 
        tools.append("Genrich") 
    elif assay_type == "chip":
        tools.extend(["macsBroad"])
        # turn sicer off for now
        # tools.extend(["macsBroad", "sicer"])
    return tools


def dedup_out7(input, assay, paired_end):
    dd = []
    if assay == "cfchip":
        dd.append(input + ".Q5DD_tagAlign")
    elif paired_end == False and assay == "chip":
        dd.append(input + ".Q5DD_tagAlign.gz")
    return dd


def get_ppqt_input(ppqt_dir, wildcards, paired_end):
    ppqt = []
    if paired_end:
        ppqt.append(join(ppqt_dir, "{0}.{1}.ppqt.txt".format(wildcards.name, wildcards.ext)))
    else:
        if wildcards.ext == "Q5DD":
            ppqt.append(join(ppqt_dir, "{0}.Q5DD_tagAlign.ppqt.txt".format(wildcards.name)))
        elif wildcards.ext == "sorted":
            ppqt.append(join(ppqt_dir, "{0}.sorted.ppqt.txt".format(wildcards.name)))
        else:
            raise ValueError(f'Unknown alignment file extension, name: {wildcards.name}, ext: {wildcards.ext}.')
    return ppqt


def get_bam_input(bam_dir, wildcards, paired_end):
    bams = []
    if paired_end:
        bams.append(join(bam_dir, "{0}.{1}.bam".format(wildcards.name, wildcards.ext)))
    else:
        if wildcards.ext == "Q5DD":
            bams.append(join(bam_dir, "{0}.Q5DD.bam".format(wildcards.name)))
        elif wildcards.ext == "sorted":
            bams.append(join(bam_dir, "{0}.sorted.bam".format(wildcards.name)))
    return bams

