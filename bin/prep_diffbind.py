#!/usr/bin/env python3
import json
import argparse
from csv import DictWriter
from os.path import join


def main(condition, peaktool, peakext, peakcaller, csvfile, wp, bam_dir):
    config = json.load(open(join(wp, "config.json")))
    chip2input = config['project']['peaks']['inputs']
    groupdata = config['project']['groups']
    blocks = config['project']['blocks']
    blocking = False if set(blocks.values()) in ({None}, {''}) else True
    
    if blocking:
        cols = ["SampleID", "Condition", "Treatment", "Replicate", "bamReads", 
                "ControlID", "bamControl", "Peaks", "PeakCaller"]
    else:
        cols = ["SampleID", "Condition", "Replicate", "bamReads", "ControlID", 
                "bamControl", "Peaks", "PeakCaller"]
    
    samplesheet = []
    # {group1}_vs_{group2} == condition
    g1, g2 = condition.split('_')[0], condition.split('_')[2]
    for group in [g1, g2]:
        for chip in groupdata[group]:
            replicate = str([ i + 1 for i in range(len(groupdata[group])) if groupdata[group][i]== chip ][0])
            bamReads = join(bam_dir, chip + ".Q5DD.bam")
            controlID = chip2input[chip]
            if controlID != "":
                bamControl = join(bam_dir, controlID + ".Q5DD.bam")
            else:
                bamControl = ""
            peaks = join(wp, peaktool, chip, chip + peakext)
            if blocking:
                block = blocks[chip]
                this_row = dict(zip(cols, [chip, group, block, replicate, bamReads, 
                    controlID, bamControl, peaks, peakcaller]))
            else:
                this_row = dict(zip(cols, [chip, group, replicate, bamReads, 
                    controlID, bamControl, peaks, peakcaller]))
            samplesheet.append(this_row)
                

    with open(csvfile, 'w') as f:
        wrt = DictWriter(f, dialect='excel', fieldnames=cols)
        wrt.writeheader()
        wrt.writerows(samplesheet)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to prepare the DiffBind input csv')
    parser.add_argument('--con', dest='contrast', required=True, help='Contrast string in [GROUP1]_vs_[GROUP2] format')
    parser.add_argument('--wp', dest='wp', required=True, help='Full path of the working directory')
    parser.add_argument('--pt', dest='peaktool', required=True, 
                        help='Name of the the peak calling tool, also the directory where the peak file will be located')
    parser.add_argument('--pe', dest='peakext', required=True, help='The file extension of the peakcall output')
    parser.add_argument('--pc', dest='peakcaller', required=True, 
                        help='Value for the PeakCaller column of the DiffBind csv')
    parser.add_argument('--bd', dest='bam_dir', required=True, 
                        help='Name of the directory where the bam files are located')
    parser.add_argument('--csv', dest='csvfile', required=True, help='Name of the output csv file')
    args = parser.parse_args()
    main(args.contrast, args.peaktool, args.peakext, args.peakcaller, args.csvfile, args.wp, args.bam_dir)
