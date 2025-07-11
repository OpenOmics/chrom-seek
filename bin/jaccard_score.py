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
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from textwrap import dedent
from pybedtools import BedTool
from sklearn.decomposition import PCA as sklearnPCA


# matplotlib
mpl.use('Agg')


##########################################
# Functions
def split_infiles(infiles):
    """
    breaks the infile string with space-delimited file names and creates a list.
    also works for infile types
    """
    infileList = infiles.strip("\'").strip('\"').split(" ")
    if len(infileList) == 1:
        infileList = infileList[0].split(";")
    return(infileList)


def loop_jaccard(infileList, genomefile):
    """
    Uses two loops to do all possible pairwise comparisons of files 
    in a list. Returns a writeable output and a pandas object
    """
    nfiles = len(infileList)
    (colnames, snames) = get_colnames(infileList)
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
            if len(outTable) == 0:
                outTable.append( "\t".join(keylist) )
            outTable.append( "\t".join(data) )
        out2 = pd.DataFrame(out, columns=colnames, index=colnames,dtype="float")
    return (outTable, out2, snames)


def run_jaccard(fileA, fileB, genomefile):
    """
    Running bedtools. Reads in two bedtools approved file types, sorts the files, 
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
    return (data, keylist)


def get_colnames(infileList):
    snames = [ i.split("/")[-1].split(".")[0].strip("_peaks").strip("_broadpeaks") for i in infileList ]
    colnames = snames
    return (colnames, snames)


def pca_plot(out, snames, peakcaller, pcatabout, outPCAFile):
    """
    creates a 2D PCA plot comparing the files based upon jaccard scores
    """
    sklearn_pca = sklearnPCA(n_components=2)
    Y_sklearn = sklearn_pca.fit_transform(out)
    PCAdata = pd.DataFrame(Y_sklearn, columns=["PC1", "PC2"])
    PCAdata["sample_name"] = snames
    PCAdata["peak_caller"] = peakcaller
    PCAdata.to_csv(pcatabout, sep='\t', index=False)

    fig, ax = plt.subplots()
    snames_pal = sns.hls_palette(len(set(snames)),s=.8)
    sns.set_palette(snames_pal)
    ax = sns.scatterplot(x="PC1", y="PC2", hue="sample_name", data=PCAdata, s=100)
    ax.axhline(y=0, color='grey', linewidth=1,linestyle="--")
    ax.axvline(x=0, color='grey', linewidth=1,linestyle="--")
    ax.set(
        xlabel= "PC1 (" + str(round(100*sklearn_pca.explained_variance_[0],2)) + "%)",
        ylabel= "PC2 (" + str(round(100*sklearn_pca.explained_variance_[1],2)) + "%)"
    )
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    plt.savefig(outPCAFile, bbox_inches='tight')
    plt.close("all")
    return


def plot_heatmap(out, outHeatmapFile, snames):
    snames_pal = sns.hls_palette(len(set(snames)),s=.8)
    snames_lut = dict(zip(set(snames), snames_pal))
    snames_cols = pd.Series(snames,index=out.index).map(snames_lut)
    g = sns.clustermap(out,cmap="YlGnBu",col_cluster=False,
                row_colors=snames_cols)
    for label in set(snames):
        g.ax_col_dendrogram.bar(0, 0, color=snames_lut[label],
                        label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc="center", ncol=3, 
                            bbox_to_anchor=(0.5, 0.8))
    plt.savefig(outHeatmapFile, bbox_inches='tight')
    plt.close("all")
    return


def write_out(out, outFile):
    with open(outFile, 'w') as f:
        f.write( "\n".join(out) )
        f.close()
    return


def main():
    desc = \
    dedent("""
    This function takes a space-delimited list of files (bed, bedgraph, gff, gtf, etc.)
    and calculates all possible pairwise jaccard scores. From bedtools: 'Jaccard is the 
    length of the intersection over the union. Values range from 0 (no intersection) to 
    1 (self intersection)'. The columns of the output file are: fileA, fileB, 
    intersection, jaccard, n_intersections, and union-intersection.
    """)

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        '-i', 
        dest='infiles', 
        required=True,
        help="A space or semi-colon delimited list of peak call input files for jaccard analysis"
    )
    parser.add_argument(
        '--caller', 
        dest='peakcaller', 
        required=True,
        help="Name of the peak caller used"
    )
    parser.add_argument(
        '--outtable', 
        dest='table', 
        required=True, 
        help='jaccard tabular output file name'
    )
    parser.add_argument(
        '--pcaplot', 
        dest='pcaplot', 
        required=True, 
        help='jaccard pca plot output file name'
    )
    parser.add_argument(
        '--pcatab', 
        dest='pcatab', 
        required=True, 
        help='jaccard pca tabular output file name'
    )
    parser.add_argument(
        '--outheatmap',
        required=True, 
        dest='heatmap', 
        help='jaccard heatmap output plot file name'
    )
    parser.add_argument(
        '-g', 
        dest='genomefile', 
        required=True,
        help='The genome contig sizes reference file'
    )

    args = parser.parse_args()

    # incoming arguments
    infiles = args.infiles
    genomefile = args.genomefile
    outTableFile = args.table
    outPCAplot = args.pcaplot
    outPCAtab = args.pcatab
    outHeatmapFile = args.heatmap
    pkcaller = args.peakcaller

    # downstream processing
    infileList = split_infiles(infiles)

    outTable, out, snames = loop_jaccard(infileList, genomefile)
    write_out(
        outTable, 
        outTableFile
    )
    pca_plot(
        out,
        snames,
        pkcaller,
        outPCAtab, 
        outPCAplot
    )
    plot_heatmap(
        out, 
        outHeatmapFile, 
        snames
    )

if __name__ == '__main__':
    main()


