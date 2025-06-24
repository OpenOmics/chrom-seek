#!/usr/bin/env python3
# ~~~ Common helper functions for grouping of outputs
from os.path import join


# common functions related to sample grouping or group meta-information
def group_samples_by_reps(groupdata, samples, chip2input):
    groupdatawinput = {}
    groupswreps = []
    for group, chipsamples in groupdata.items():
        tmp = []
        if len(chipsamples) > 1:
            groupswreps.append(group)
        for chip in chipsamples:
            if chip in samples:
                tmp.append(chip)
                input = chip2input.get(chip, "")
                if input != "NA" and input != "":
                    tmp.append(input)
        if len(tmp) != 0:
            groupdatawinput[group] = set(tmp)
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
    elif assay_type == "cutnrun":
        tools.append("SEACR")
    return tools


def test_for_block(groupdata, contrast, blocks):
    """only want to run blocking on contrasts where all
    individuals are on both sides of the contrast"""
    contrastBlock = []
    for con in contrast:
        group1 = con[0]
        group2 = con[1]
        block1 = [blocks[sample] for sample in groupdata[group1]]
        block2 = [blocks[sample] for sample in groupdata[group2]]
        if len(block1) == len(block2):
            if len(set(block1).intersection(block2)) == len(block1):
                contrastBlock.append(con)
    return contrastBlock
