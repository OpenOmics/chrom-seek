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


def get_peaktools(assay_type):
    tools = ["macsNarrow"]
    if assay_type == "atac": 
        tools.append("Genrich") 
    elif assay_type == "chip":
        tools.append("macsBroad")
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


def test_for_block(groupdata, contrast, blocks):
    """ only want to run blocking on contrasts where all
    individuals are on both sides of the contrast """
    contrastBlock = [ ]
    for con in contrast:
        group1 = con[0]
        group2 = con[1]
        block1 = [ blocks[sample] for sample in groupdata[group1] ]
        block2 = [ blocks[sample] for sample in groupdata[group2] ]
        if len(block1) == len(block2):
            if len(set(block1).intersection(block2)) == len(block1):
                contrastBlock.append(con)
    return contrastBlock


def ctrl_test(ctrl_dict, input_name, in_dir, mode=None):
    sample = join(in_dir, f"{input_name}.Q5DD.RPGC.bw")
    assert mode in ('chip', 'ctrl'), 'Unrecognized input file mode.'
    
    if input_name in ctrl_dict:
        norm = join(in_dir, ctrl_dict[input_name] + ".Q5DD.RPGC.bw")
    else:
        raise ValueError(f'ChIP sample {input_name} missing from input lookup: \n{str(ctrl_dict)}')
    outs = {'chip': sample, 'ctrl': norm}
    return outs[mode]