#!/usr/bin/env python3

"""
Name: jaccard_score.py
Created by: Tovah Markowitz
Date: 1/23/19
Updated: 8/5/19 to compare multiple tools and create plots

Purpose: To do all pairwise comparisons of bed/peak files given. Uses bedtools
to calculate a jaccard score for every comparison. All data is saved in a 
single tab-delimited file.
"""

##########################################
# Modules
import optparse
from pybedtools import BedTool
import pandas as pd
from sklearn.decomposition import PCA as sklearnPCA
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

##########################################
# Functions

def split_infiles(infiles):
    """ breaks the infile string with space-delimited file names and creates a list.
    also works for infile types
    """
    infileList = infiles.strip("\'").strip('\"').split(" ")
    if len(infileList) == 1:
        infileList = infileList[0].split(";")
    return(infileList)

def loop_jaccard(infileList, genomefile, filetypeList):
    """ Uses two loops to do all possible pairwise comparisons of files 
    in a list. Returns a writeable output and a pandas object
    """
    nfiles = len(infileList)
    (colnames, snames) = get_colnames(infileList, filetypeList)
    out = [[1.000] * nfiles for i in range(nfiles)]
    outTable = []
    for z in range(nfiles):
        fileA = infileList[z]
        print("fileA is: " + fileA) 
        for y in range(z+1,nfiles):
            fileB = infileList[y]
            (data, keylist) = run_jaccard(fileA, fileB, genomefile)
            out[z][y] = data[3]
            out[y][z] = data[3]
            if filetypeList != [""]:
                keylist.insert(1, "toolA")
                keylist.insert(3, "toolB")
                data.insert(1, filetypeList[z])
                data.insert(3, filetypeList[y])
            if len(outTable) == 0:
                outTable.append( "\t".join(keylist) )
            outTable.append( "\t".join(data) )
        out2 = pd.DataFrame(out, columns=colnames, index=colnames,dtype="float")
    return(outTable, out2, snames)

def run_jaccard(fileA, fileB, genomefile):
    """ Running bedtools. Reads in two bedtools approved file types, sorts the files, 
    and calculates a jaccard score.
    """
    a = BedTool(fileA)
    a = a.sort(g=genomefile)
    b = BedTool(fileB)
    b = b.sort(g=genomefile)
    j = a.jaccard(b,g=genomefile)
    j["fileA"] = fileA.split("/")[-1]
    j["fileB"] = fileB.split("/")[-1]
    keylist = list(j.keys())
    keylist.sort()
    data = [ str(j[key]) for key in keylist ]
    return(data, keylist)

def get_colnames(infileList, filetypeList):
    snames = [ i.split("/")[-1].split(".")[0].strip("_peaks").strip("_broadpeaks") for i in infileList ]
    if filetypeList == [""]:
        colnames = snames
    else:
        colnames = [ snames[i] + "_" + filetypeList[i] for i in range(len(snames)) ]
    return(colnames, snames)

def create_outfile_names(outroot):
    """ uses outroot to create the output file names """
    outTableFile = "jaccard.txt"
    outPCAFile = "jaccard_PCA.pdf"
    outHeatmapFile = "jaccard_heatmap.pdf"
    if outroot != "":
        if outroot[-1] == "/":
            outTableFile= outroot + outTableFile
            outPCAFile = outroot + outPCAFile
            outHeatmapFile = outroot + outHeatmapFile
        else:
            outTableFile= outroot + "_" + outTableFile
            outPCAFile = outroot + "." + outPCAFile
            outHeatmapFile = outroot + "." + outHeatmapFile
    return(outTableFile, outPCAFile, outHeatmapFile)

def pca_plot(out, filetypeList, snames, outPCAFile):
    """ creates a 2D PCA plot comparing the files based upon jaccard scores
    """
    sklearn_pca = sklearnPCA(n_components=2)
    Y_sklearn = sklearn_pca.fit_transform(out)
    PCAdata = pd.DataFrame(Y_sklearn,columns=["PC1","PC2"])
    PCAdata.insert(0,"sample name",snames)
    fig, ax =plt.subplots()
    snames_pal = sns.hls_palette(len(set(snames)),s=.8)
    sns.set_palette(snames_pal)
    if filetypeList != [""]:
        PCAdata.insert(1,"tool",filetypeList)
        ax = sns.scatterplot(x="PC1",y="PC2",hue="sample name",style="tool",data=PCAdata,s=100)
    else:
        ax = sns.scatterplot(x="PC1",y="PC2",hue="sample name",data=PCAdata,s=100)
    ax.axhline(y=0, color='grey', linewidth=1,linestyle="--")
    ax.axvline(x=0, color='grey', linewidth=1,linestyle="--")
    ax.set(xlabel= "PC1 (" + str(round(100*sklearn_pca.explained_variance_[0],2)) + "%)",
           ylabel= "PC2 (" + str(round(100*sklearn_pca.explained_variance_[1],2)) + "%)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    #plt.show()
    plt.savefig(outPCAFile, bbox_inches='tight')
    plt.close("all")

def plot_heatmap(out, outHeatmapFile, snames, filetypeList):
    snames_pal = sns.hls_palette(len(set(snames)),s=.8)
    snames_lut = dict(zip(set(snames), snames_pal))
    snames_cols = pd.Series(snames,index=out.index).map(snames_lut)
    if filetypeList != [""]:
       tool_pal = sns.cubehelix_palette(len(set(filetypeList)))
       tool_lut = dict(zip(set(filetypeList), tool_pal))
       tool_cols = pd.Series(filetypeList,index=out.index).map(tool_lut)
       g = sns.clustermap(out,cmap="YlGnBu",col_cluster=False,
                    row_colors=[snames_cols,tool_cols])
       for label in set(snames):
            g.ax_col_dendrogram.bar(0, 0, color=snames_lut[label],
                            label=label, linewidth=0)
       for label in set(filetypeList):
            g.ax_col_dendrogram.bar(0, 0, color=tool_lut[label],
                            label=label, linewidth=0)
       g.ax_col_dendrogram.legend(loc="center", ncol=3, 
                                bbox_to_anchor=(0.4, 0.8))
    else:
       g = sns.clustermap(out,cmap="YlGnBu",col_cluster=False,
                    row_colors=snames_cols)
       for label in set(snames):
            g.ax_col_dendrogram.bar(0, 0, color=snames_lut[label],
                            label=label, linewidth=0)
       g.ax_col_dendrogram.legend(loc="center", ncol=3, 
                                bbox_to_anchor=(0.5, 0.8))
    #plt.show()
    plt.savefig(outHeatmapFile, bbox_inches='tight')
    plt.close("all")

def write_out(out, outFile):
    f = open(outFile, 'w')
    f.write( "\n".join(out) )
    f.close()

##########################################
# Main

def main():
    desc="""
    This function takes a space-delimited list of files (bed, bedgraph, gff, gtf, etc.)
    and calculates all possible pairwise jaccard scores. From bedtools: 'Jaccard is the 
    length of the intersection over the union. Values range from 0 (no intersection) to 
    1 (self intersection)'. The columns of the output file are: fileA, fileB, 
    intersection, jaccard, n_intersections, and union-intersection.
    """

    parser = optparse.OptionParser(description=desc)

    parser.add_option('-i', dest='infiles', default='', help='A space- or semicolon-delimited list of \
input files for analysis.')
    parser.add_option('-t', dest='filetypes', default='', help='A space- or semicolon-delimited list \
of input file sources/types.')
    parser.add_option('-o', dest='outroot', default='', help='The root name of the output files \
where all the jaccard score information will be saved.')
    parser.add_option('-g', dest='genomefile', default='', help='The name of the .genome file.')

    (options,args) = parser.parse_args()
    infiles = options.infiles
    filetypes = options.filetypes
    outroot = options.outroot
    genomefile = options.genomefile

    infileList = split_infiles(infiles)
    filetypeList = split_infiles(filetypes)
    (outTable, out, snames) = loop_jaccard(infileList, genomefile, filetypeList)
    (outTableFile, outPCAFile, outHeatmapFile) = create_outfile_names(outroot)
    write_out(outTable, outTableFile)
    pca_plot(out, filetypeList, snames, outPCAFile)
    plot_heatmap(out, outHeatmapFile, snames, filetypeList)

if __name__ == '__main__':
    main()


