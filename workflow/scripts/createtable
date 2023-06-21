#!/usr/bin/env python

"""
Author:		Skyler Kuhn
Date created:	08/18/2018
Last modified:	04/26/2019 by Tovah Markowitz: markowitzte@nih.gov
Email:		kuhnsa@nih.gov

About: This program takes standard input from the concatenation of each *qc.metric file.The results 
       are then stored into a nested dictionary (where: dict[sample][qcmetric] = value). This data 
       structure is then converted into a pandas dataframe. 
       Example Usage:
       --------------  
	   cat sample1.qcmetrics sample2.qcmetrics sampleNth.qcmetrics | ./createQCTable > ChIPseq_QC_Table.txt
Python version(s):	2.7 or 3.X
"""

from __future__ import print_function
import pandas as pd
import sys

def file2table():
    tabledict = {}
    for line in sys.stdin:
        linelist = line.strip().split()
        #print(linelist)
        if linelist[1] not in tabledict:
            tabledict[linelist[1]] = {}
        try:
            tabledict[linelist[1]][linelist[0]] = linelist[2]
        except IndexError:
            pass
    
    df = pd.DataFrame(tabledict)
    df.index.name = 'SampleName'
    df.reset_index(inplace=True)
    #print(df[['NSC', 'FRiP', 'PCB1', 'PCB2', 'RSC']])  #re-order columns
    #cols = df.columns.tolist() # view df columns names
    #orderedcols = ordercolumns(cols)
    #print(df.to_string())
    print(df[['SampleName', 'NReads', 'NMappedReads', 'NUniqMappedReads', 'NRF', 'PBC1', 'PBC2', 'FragmentLength', 'NSC', 'RSC', 'Qtag']].to_string(index=False,justify="left"))


if __name__ == '__main__':
    file2table()

