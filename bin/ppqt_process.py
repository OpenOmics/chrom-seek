#!/usr/bin/env python3

#Purpose: To grab the estimated fragment length from the ppqt output and a small txt with that information. For input files, adding an extra value of 200bp as an alternative.
import argparse
parser = argparse.ArgumentParser(description='Script to extract the the estimated fragment length from the ppqt output.')
parser.add_argument('-i', required=True,help='Name of the ppqt txt file')
parser.add_argument('-o', required=True,help='Name of the output file')
args = parser.parse_args()

output = args.o
inppqt = args.i

o=open(output,'w')

file = list(map(lambda z:z.strip().split(),open(inppqt,'r').readlines()))


ppqt_values = file[0][2].split(",")
extenders = []
for ppqt_value in ppqt_values:
    if int(ppqt_value) > 150:
        extenders.append(ppqt_value)
if len(extenders) > 0:
    o.write(extenders[0])
else:
    o.write("200")      
o.close()