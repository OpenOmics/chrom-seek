#!/usr/bin/env python3

import json
import argparse

parser = argparse.ArgumentParser(description='Script to prepare the DiffBind input csv')
parser.add_argument('--wp',dest='workpath',required=True,help='Full path of the working directory')
parser.add_argument('--pt',dest='peaktool',required=True,help='Name of the the peak calling tool, also the directory where the peak file will be located')
parser.add_argument('--pe',dest='peakextension',required=True,help='The file extension of the peakcall output')
parser.add_argument('--pc',dest='peakcaller',required=True,help='Value for the PeakCaller column of the DiffBind csv')
parser.add_argument('--bd',dest='bamdir',required=True,help='Name of the directory where the bam files are located')
parser.add_argument('--csv',dest='csvfile',required=True,help='Name of the output csv file')

args = parser.parse_args()

with open("config.json","r") as read_file:
   config=json.load(read_file)
   
chip2input = config['project']['peaks']['inputs']
groupdata = config['project']['groups']

tmpIDs = [x for xs in groupdata.values() for x in xs]
Ncounts = [tmpIDs.count(tmp) for tmp in set(tmpIDs)]

samplesheet = [",".join(["SampleID","Condition", "Replicate", "bamReads", 
         "ControlID", "bamControl", "Peaks", "PeakCaller"])]

count = 1
for chip in chip2input.keys():
   if set(Ncounts) == {1}: # if all samples only in one group
      for key in groupdata.keys():
          if chip in groupdata[key]:
             condition = key
      replicate = str([ i + 1 for i in range(len(groupdata[condition])) if groupdata[condition][i]== chip ][0])
   else:
      condition = ""
      replicate = str(count)
      count = count +1
   bamReads = args.bamdir + "/" + chip + ".Q5DD.bam"
   controlID = chip2input[chip]
   if controlID != "":
      bamControl = args.bamdir + "/" +  controlID + ".Q5DD.bam"
   else:
      bamControl = ""
   peaks = args.workpath + "/" + args.peaktool + "/" + chip + "/" + chip + args.peakextension
   samplesheet.append(",".join([chip, condition, replicate, bamReads, 
                   controlID, bamControl, peaks, args.peakcaller]))            

f = open(args.csvfile, 'w')
f.write ("\n".join(samplesheet))
f.close()
