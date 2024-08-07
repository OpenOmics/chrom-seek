#!/usr/bin/env python
"""
Author:		Skyler Kuhn
Date created:	08/18/2018
Last modified:	09/10/2019 by Tovah Markowitz markowitzte@nih.gov
Email:		kuhnsa@nih.gov

About: This program takes in a parsed data from Samtools flagstat, atac_nrf.py, ngsqc,
       PhantomPeakQuailty tools, or atac_frip.py and parses it further to append it to
       a nested dictionary (which is store a json object). The program is designed to 
       work with standarad input (grep 'pattern' filename.txt | filerMetrics sampleName filterType).
       
       Example Usage:
       --------------  
	1.) Find total number of reads
          # grep 'in total' H3k4me3_gran_1.sorted.bam.flagstat | awk '{print $1,$3}' | ./filterMetrics H3k4me3_gran_1 tnreads	
	
	2.) Find total number of mapped reads
	  # grep 'mapped (' H3k4me3_gran_1.sorted.bam.flagstat | awk '{print $1,$3}' | ./filterMetrics H3k4me3_gran_1 mnreads  
	
	3.) Find total number of uniquely mapped reads
	  # grep 'mapped (' H3k4me3_gran_1.sorted.Q5DD.bam.flagstat | awk '{print $1,$3}' | ./filterMetrics H3k4me3_gran_1 unreads
	
	4.) Find NRF, PCB1, PCB2 
	  # cat H3k4me3_gran_1.nrf | ./filterMetrics H3k4me3_gran_1 nrf
	
	7.) Find FRiP (in the second half of ChIP-seq pipeline)
	  # TO-DO
	
	8.) Find NGSQC statistics (detla RCIs)
	  # grep '<' NGSQC_report.txt | awk '{print $(NF)}' | xargs | ./filterMetrics H3k4me3_gran_1 ngsqc
	
	9.) Find NSC, RSC, Qtag
	  # awk '{print $(NF-2),$(NF-1),$(NF)}' H3k4me3_gran_1.sorted.Q5DD.ppqt | ./filterMetrics H3k4me3_gran_1 ppqt

	10.) Find the Fragment Length
	  # awk -F '\t' '{print $3}' H3k4me3_gran_1.sorted.Q5DD.ppqt | sed -e 's/,/ /g' | ../Scripts/filterMetrics H3k4me3_gran_1 fragLen

Python version(s):	2.7 or 3.X
"""

from __future__ import print_function, division
import sys


def add(i1, i2):
	return int(i1) + int(i2)

def divide(num):
	return int(num/4.0)

def getmetadata(type):
	metadata = []
	if type == 'ppqt':
		metadata = ["NSC", "RSC", "Qtag"]
	elif type == 'ngsqc':
		metadata = ["NGSQC_dRCI_2.5", "NGSQC_dRCI_5.0", "NGSQC_dRCI_10.0"]
	elif type == 'nrf':
		metadata = ["NRF", "PBC1", "PBC2"]
	elif type == 'tnreads':
		metadata = 'NReads'
	elif type == 'mnreads':
		metadata = 'NMappedReads'
	elif type == 'unreads':
		metadata = 'NUniqMappedReads'
	elif type == 'fragLen':
		metadata = 'FragmentLength'
	return metadata
	

def filteredData(sample, ftype):
	""" 
	Data grabbed by the awk or grep commands in the above example cases becomes 
	variable 'line' and gets split by space to make 'linelist'
	"""
	for line in sys.stdin:
		linelist = line.strip().split()
		if 'reads' in ftype:
			mtypes = getmetadata(ftype)
			if ftype == 'tnreads':
				v1 = int(linelist[0])
				print("{}\t{}\t{}".format(sample, mtypes, divide(v1)))
			else:
				v1, v2 = linelist[0], linelist[1]
				print("{}\t{}\t{}".format(sample, mtypes, add(v1,v2)))
		elif ftype == 'fragLen':
			mtypes = getmetadata(ftype)
			extenders = []
			for ppqt_value in linelist:
				if int(ppqt_value) > 150:
					extenders.append(ppqt_value)
			if len(extenders) > 0:
				print("{}\t{}\t{}".format(sample, mtypes, extenders[0]))
			else:
				print("{}\t{}\t{}".format(sample, mtypes, linelist[0]))
		elif ftype == 'ppqt' or ftype == 'ngsqc' or ftype == 'nrf':
			mtypes = getmetadata(ftype)
			for i in range(len(linelist)):
				print("{}\t{}\t{}".format(sample, mtypes[i], linelist[i]))

def main():
	# Grab Arguements
	arguements = sys.argv[1:]
	# Sample Name
	samplename = arguements[0]

	# Filter Type (extremely important: determines the type of filtering to complete)
	filter = arguements[1]
	filteredData(samplename, filter)


if __name__ == '__main__':
	main()
