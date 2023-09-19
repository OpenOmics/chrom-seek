#!/usr/bin/env python3

Purpose: To grab the estimated fragment length from the ppqt output and a small txt with that information. For input files, adding an extra value of 200bp as an alternative.

parser = argparse.ArgumentParser(description='Script to extract the the estimated fragment length from the ppqt output.')
parser.add_argument('-i',dest='inppqt',required=True,help='Name of the ppqt txt file')
parser.add_argument('-o',dest='output',required=True,help='Name of the output file')
parser.add_argument('--input',dest='inputSample',required=True,help='True or False; whether the sample is considered an input or not')
args = parser.parse_args()

o=open(output,'w')
file = list(map(lambda z:z.strip().split(),open(inppqt,'r').readlines()))
ppqt_values = file[0][2].split(",")
extenders = []
for ppqt_value in ppqt_values:
  if int(ppqt_value) > 150:
     extenders.append(ppqt_value)
if len(extenders) > 0:
  o.write(extenders[0] + "\n")
else:
    print("All estimated fragments lengths were less than 150bp which will may cause the pipeline to fail.")
    print("Potential causes include: wrong ref genome selected or low starting DNA.")
    print("Assuming default estimated fragment length of 200bp.\n")
    o.write("200" + "\n")
#if inputSample:
    o.write("200" + "\n")
o.close()
