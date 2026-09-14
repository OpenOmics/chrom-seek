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
import math
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from textwrap import dedent
from pybedtools import BedTool
from sklearn.decomposition import PCA as sklearnPCA


# matplotlib
mpl.use('Agg')


# How an incalculable jaccard score is reported, and the columns of the
# tabular output: the bedtools keys of a jaccard record, sorted, plus the
# two file names.
NA_VALUE = "NA"
TABLE_COLUMNS = (
    "fileA", "fileB", "intersection", "jaccard", "n_intersections", "union-intersection"
)


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


def peak_file_has_intervals(path):
    """True when a peak file has at least one non-comment interval line."""
    if not path or (not os.path.isfile(path)):
        return False
    with open(path, "r") as peak_file:
        for line in peak_file:
            line = line.strip()
            if line and not line.startswith("#") and not line.startswith("track") and not line.startswith("browser"):
                return True
    return False


def write_placeholder_plot(outfile, title, message):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.axis("off")
    ax.set_title(title)
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    plt.savefig(outfile, bbox_inches='tight')
    plt.close("all")


def parse_score(score):
    """
    Converts a jaccard score from bedtools into a float. Scores that are not a
    number stay NaN, so that they are reported as NA rather than as a real
    score: the "nan" string bedtools prints for an empty union, numpy.nan,
    numpy float NaNs, pandas.NA and None are all caught here.
    """
    try:
        if pd.isna(score):
            return float("nan")
    except (TypeError, ValueError):
        # non-scalar or otherwise untestable value, fall through to float()
        pass
    try:
        return float(score)
    except (TypeError, ValueError):
        return float("nan")


def na_record(fileA, fileB):
    """
    Builds the tabular output row of a comparison that cannot be calculated,
    that is, one where a peak file has no peaks. Every score is NA.
    """
    record = dict.fromkeys(TABLE_COLUMNS, NA_VALUE)
    record["fileA"] = fileA.split("/")[-1]
    record["fileB"] = fileB.split("/")[-1]
    keylist = list(TABLE_COLUMNS)
    return ([record[key] for key in keylist], keylist)


def drop_na_samples(out, snames):
    """
    Keeps only the samples that have a jaccard score for every remaining
    comparison. Samples with an NA score, that is, peak files with no peaks,
    cannot be plotted or clustered. Note that an NA sample puts an NA in every
    row of the matrix, so the samples themselves have to be dropped, not the
    rows that mention them.
    """
    # a peak file with no peaks has no score against itself either
    keep = [i for i in range(out.shape[0]) if not pd.isna(out.iloc[i, i])]
    # any NA left in the block belongs to a single pair of samples, drop
    # whichever of them accounts for the most NAs and look again
    while keep:
        na_counts = out.iloc[keep, keep].isna().sum(axis=1).tolist()
        if max(na_counts) == 0:
            break
        keep.pop(na_counts.index(max(na_counts)))
    return (out.iloc[keep, keep], [snames[i] for i in keep])


def loop_jaccard(infileList, genomefile):
    """
    Uses two loops to do all possible pairwise comparisons of files
    in a list. Returns a writeable output and a pandas object. Any comparison
    involving a peak file with no peaks is reported as NA instead of being
    handed to bedtools, which cannot score it.
    """
    nfiles = len(infileList)
    (colnames, snames) = get_colnames(infileList)
    has_peaks = [peak_file_has_intervals(infile) for infile in infileList]
    out = [[float("nan")] * nfiles for i in range(nfiles)]
    for i in range(nfiles):
        if has_peaks[i]:
            out[i][i] = 1.0
    outTable = [ "\t".join(TABLE_COLUMNS) ]
    for z in range(nfiles):
        fileA = infileList[z]
        print("fileA is: " + fileA)
        for y in range(z+1,nfiles):
            fileB = infileList[y]
            if has_peaks[z] and has_peaks[y]:
                (data, keylist) = run_jaccard(fileA, fileB, genomefile)
                score = parse_score(data[keylist.index("jaccard")])
                if math.isnan(score):
                    data[keylist.index("jaccard")] = NA_VALUE
            else:
                (data, keylist) = na_record(fileA, fileB)
                score = float("nan")
            out[z][y] = score
            out[y][z] = score
            outTable.append( "\t".join(data) )
    out2 = pd.DataFrame(out, columns=colnames, index=colnames, dtype="float")
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


def strip_suffix(text, *suffixes):
    """Remove trailing suffixes from a string, one full suffix at a time.

    Unlike str.strip(), which removes any leading/trailing characters found in
    its argument, this only removes the exact suffix(es) given.
    """
    for suffix in suffixes:
        if suffix and text.endswith(suffix):
            text = text[: -len(suffix)]
    return text


def get_colnames(infileList):
    snames = [ strip_suffix(i.split("/")[-1].split(".")[0], "_peaks", "_broadpeaks") for i in infileList ]
    colnames = snames
    return (colnames, snames)


def load_sample_groups(config_file):
    """Load sample-to-group(s) mapping from pipeline config, if available.

    Returns a dict mapping each sample to the list of groups it belongs to.
    Samples can appear in more than one group, so memberships are accumulated
    rather than overwritten.
    """
    if not config_file:
        return {}

    try:
        with open(config_file, "r") as handle:
            cfg = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}

    groups = cfg.get("project", {}).get("groups", {})
    sample2groups = {}
    for grp, sample_list in groups.items():
        for sample in sample_list:
            sample2groups.setdefault(sample, []).append(grp)
    return sample2groups


def resolve_group_coloring(snames, sample2groups):
    """Decide how to colorize plots based on group membership.

    If every sample belongs to exactly one group, colorize by group label and
    return that per-sample group list. Otherwise (any sample has zero or
    multiple group memberships), fall back to colorizing by sample name.

    Returns (use_group_coloring, labels), where labels is the per-sample list of
    labels to color by.
    """
    if sample2groups:
        memberships = [sample2groups.get(s, []) for s in snames]
        if all(len(m) == 1 for m in memberships):
            return True, [m[0] for m in memberships]
    return False, list(snames)


def pca_plot(out, snames, peakcaller, pcatabout, outPCAFile, sample2groups=None):
    """
    creates a 2D PCA plot comparing the files based upon jaccard scores.
    Only samples with a score for every comparison are plotted, NA samples
    are left out.
    """
    (out, snames) = drop_na_samples(out, snames)

    if out.shape[0] < 2 or out.shape[1] < 2:
        PCAdata = pd.DataFrame({
            "PC1": [0.0] * len(snames),
            "PC2": [0.0] * len(snames),
            "sample_name": snames,
            "peak_caller": [peakcaller] * len(snames),
        })
        PCAdata.to_csv(pcatabout, sep='\t', index=False)
        write_placeholder_plot(
            outPCAFile,
            f"{peakcaller} Jaccard PCA",
            "Insufficient valid (non-NA) samples for PCA."
        )
        return

    sklearn_pca = sklearnPCA(n_components=2)
    Y_sklearn = sklearn_pca.fit_transform(out)
    PCAdata = pd.DataFrame(Y_sklearn, columns=["PC1", "PC2"])
    PCAdata["sample_name"] = snames
    PCAdata["peak_caller"] = peakcaller
    if sample2groups is None:
        sample2groups = {}
    PCAdata["group"] = [";".join(sample2groups.get(s, [])) for s in snames]
    use_group_coloring, _ = resolve_group_coloring(snames, sample2groups)
    PCAdata.to_csv(pcatabout, sep='\t', index=False)

    fig, ax = plt.subplots()
    if use_group_coloring:
        group_pal = sns.hls_palette(len(set(PCAdata["group"])), s=.8)
        sns.set_palette(group_pal)
        ax = sns.scatterplot(x="PC1", y="PC2", hue="group", style="sample_name", data=PCAdata, s=100)
    else:
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


def plot_heatmap(out, outHeatmapFile, peakcaller, heatmap_tab, snames):
    """
    clusters and plots the jaccard score matrix. Every sample is written to
    the tabular output, including the NA ones, but only samples with a score
    for every comparison can be clustered.
    """
    # the full matrix, NA samples included, keeps the columns of this table
    # aligned across peak callers for jaccard_summary.py
    out_hm = out.copy()
    out_hm['peakcaller'] = peakcaller
    out_hm.to_csv(heatmap_tab, sep='\t', index=False, na_rep=NA_VALUE)

    (out, snames) = drop_na_samples(out, snames)

    if out.shape[0] < 2 or out.shape[1] < 2:
        write_placeholder_plot(
            outHeatmapFile,
            f"{peakcaller} Jaccard Heatmap",
            "Insufficient valid (non-NA) samples for heatmap clustering."
        )
        return

    snames_pal = sns.hls_palette(len(set(snames)),s=.8)
    snames_lut = dict(zip(set(snames), snames_pal))
    snames_cols = pd.Series(snames, index=out.index).map(snames_lut)
    g = sns.clustermap(out, cmap="YlGnBu", col_cluster=False, row_colors=snames_cols)
    for label in set(row_labels):
        g.ax_col_dendrogram.bar(0, 0, color=row_lut[label],
                        label=label, linewidth=0)
    g.ax_col_dendrogram.legend(title=legend_title, loc="center", ncol=3, 
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
    intersection, jaccard, n_intersections, and union-intersection. Peak files
    with no peaks cannot be scored, their comparisons are reported as NA and
    they are left out of the PCA and heatmap plots.
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
        '--tabheatmap',
        required=True, 
        dest='heatmap_tab', 
        help='jaccard heatmap output plot file name'
    )
    parser.add_argument(
        '-g', 
        dest='genomefile', 
        required=True,
        help='The genome contig sizes reference file'
    )
    parser.add_argument(
        '-c',
        '--config',
        dest='config_file',
        required=False,
        default=None,
        help='Pipeline config.json file for group annotations'
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
    hm_tsv = args.heatmap_tab
    config_file = args.config_file

    # downstream processing
    infileList = split_infiles(infiles)
    sample2groups = load_sample_groups(config_file)

    # Inputs without any peaks cannot be scored, their comparisons are
    # reported as NA and they are left out of the plots.
    no_peak_files = [f for f in infileList if not peak_file_has_intervals(f)]

    if no_peak_files:
        print(
            "WARNING: Peak files with no peaks, scored as " + NA_VALUE + ": "
            + ", ".join(no_peak_files)
        )

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
        outPCAplot,
        sample2groups=sample2groups
    )
    plot_heatmap(
        out, 
        outHeatmapFile,
        pkcaller,
        hm_tsv,
        snames,
        sample2groups=sample2groups
    )

if __name__ == '__main__':
    main()


