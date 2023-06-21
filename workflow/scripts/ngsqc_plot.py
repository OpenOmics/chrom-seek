#!/usr/bin/env python3

"""
Name: ngsqc_plot.py
Created by: Tovah Markowitz
Date: 1/3/19
Updated: 4/18/19

Purpose: To take a number of NGSQC_report.txt files and convert into a figure.
Currently designed to work in case where all files have been moved from their original
subdirectories and the filename includes information of about the source data. All input
files must be in the noted directory (not recursive), end with ".NGSQC_report.txt", and
contain 'ext'. Output figures are named based upon 'ext'. Now also creates a single output
table that can be read into the multiqc report.
"""

import matplotlib
matplotlib.use('Agg')
import optparse
import os
import re
import matplotlib.pyplot as plt


# Functions
def read_ngsqc(ngsqcFile):
    """Purpose: read in ngsqc file and extract similarity QCi data
    """
    f = open(ngsqcFile, 'r')
    ngsqc = f.readlines()
    f.close()
    ngsqc2 = [[],[]]
    for i in range( len(ngsqc) ):
        if re.search( "<", ngsqc[i] ):
            tmp = ngsqc[i].strip()
            tmp2 = re.search( "RCI < (.*)%: (.*)", tmp )
            ngsqc2[0].append(float(tmp2.group(1)))
            ngsqc2[1].append(float(tmp2.group(2)))
    return ngsqc2


def create_plot(ngsqcAll, ext, outfileName):
	""" save data as a figure.
	s90 is subsampling 90% and S50 is subsampling 50%.
	RCI (Read count intensity)
	"""
	for ind in ngsqcAll:
		plt.plot(ind[1][0],ind[1][1],'.-',label=ind[0])	
	plt.xlabel('percent RCI difference')
	plt.ylabel('Similarity QCi (s90/s50)')
	plt.title(ext)
	plt.legend()
	# plt.show()
	plt.savefig(outfileName)
	plt.close("all")


def write_ngsqc(ngsqcAll,directory,ext):
        """to write a set of tables of the ngsqc data used in the image for incorporation
        into multiqc
        """
        for sample in ngsqcAll:
            if ext == "":
                outfile = directory + "/" + sample[0] + ".NGSQC.txt"
            else:
                outfile = directory + "/" + sample[0] + "." + ext + ".NGSQC.txt"
            f = open(outfile, 'w')
            f.write( "\n".join( str(sample[1][0][i]) + "," + str(sample[1][1][i]) for i in range(len(sample[1][0])) ) )
            f.close()


def all_ngsqc(directory, ext):
        """ To get the ngsqc data for all NGSQC_report.txt files within a single folder
        assumes that all files have been moved from their original subdirectories and 
        filename include information about the source of the data
        """
        files = os.listdir(directory)
        files2 = [ file for file in files if re.search("NGSQC_report.txt",file) ]
        files3 = [ file for file in files2 if re.search(ext,file) ]
        ngsqcAll = [ ]
        if ext == "":
            regex = r"(.*).NGSQC_report.txt"
        else:
            regex = r"(.*)." + re.escape(ext) + r".NGSQC_report.txt"
        for file in files3:
            tmp = re.search(regex,file)
            ngsqcAll.append([ tmp.group(1), read_ngsqc(directory + "/" + file) ])
        return ngsqcAll


def name_outimage(directory, ext, group):
	"""this version uses ext to create the output file name
	"""
	if group == "":
		outfileName =  directory + "/NGSQC." + ext + ".png"
	else:
		outfileName =  directory + "/" + group + ".NGSQC." + ext + ".png"
	return outfileName


# Main
def main():
	desc="""
	To take a number of NGSQC_report.txt files and convert into a figure.
	Currently designed to work in case where all files have been moved from their original
	subdirectories and the filename includes information of about the source data. All input
	files must be in the noted directory (not recursive), end with ".NGSQC_report.txt", and
	contain 'ext'. Output figures are named based upon 'ext'.
	"""

	parser = optparse.OptionParser(description=desc)

	parser.add_option('-d', dest='directory', default='', help='The name of the directory \
with the NGSQC_report.txt files.')
	parser.add_option('-e', dest='ext', default='', help='The string found in all \
NGSQC_report.txt files. This string must be directly before ".NGSQC_report.txt".')
	parser.add_option('-g', dest='group', default='', help='Any additional naming to add \
the output file name.')

	(options,args) = parser.parse_args()
	directory = options.directory
	ext = options.ext
	group = options.group

	outfileName = name_outimage(directory, ext, group)
	ngsqcAll = all_ngsqc(directory, ext)
	write_ngsqc(ngsqcAll, directory, ext)
	create_plot(ngsqcAll, ext, outfileName)


if __name__ == '__main__':
    main()
