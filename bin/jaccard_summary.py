#!/usr/bin/env python3
"""
Name: jaccard_summary.py
Date: 7/10/2025

Purpose: Combine all individual jaccard PCA coordinates into one figure
         colorized by peak caller.
"""
import argparse
import os
import pandas as pd
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from textwrap import dedent


def write_placeholder_plot(outfile, title, message):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.axis("off")
    ax.set_title(title)
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    plt.savefig(outfile, bbox_inches='tight')
    plt.close("all")


def drop_na_samples(hm_full, peakcallers):
    """
    Samples with an NA jaccard score, that is, samples with no peaks, cannot
    be clustered. The per-peak-caller matrices are square and written without
    an index, so a sample's row sits at the position of its column and every
    peak caller contributes the same number of rows. An NA sample puts an NA
    in every row of its peak caller's block, so samples are spotted by their
    own self comparison on the diagonal rather than by scanning columns.
    """
    snames = hm_full.columns.tolist()
    nsamples = len(snames)
    nblocks, remainder = divmod(len(hm_full), nsamples)
    if remainder != 0:
        print("WARNING: Unexpected number of jaccard heatmap rows, keeping all samples")
        return (hm_full, peakcallers)

    keep = [
        i for i in range(nsamples)
        if not any(
            pd.isna(hm_full.iloc[block * nsamples + i, i]) for block in range(nblocks)
        )
    ]
    # what is left can still hold NAs belonging to a single pair of samples,
    # drop whichever of them accounts for the most NAs and look again
    while keep:
        rows = [block * nsamples + i for block in range(nblocks) for i in keep]
        na_counts = hm_full.iloc[rows, keep].isna().sum(axis=0).tolist()
        if max(na_counts) == 0:
            break
        keep.pop(na_counts.index(max(na_counts)))

    if len(keep) == nsamples:
        return (hm_full, peakcallers)

    dropped = [sname for i, sname in enumerate(snames) if i not in keep]
    print("WARNING: Samples without a jaccard score for every comparison, "
          "left out of the summary heatmap: " + ", ".join(dropped))

    keep_rows = [block * nsamples + i for block in range(nblocks) for i in keep]
    return (hm_full.iloc[keep_rows, keep], [peakcallers[i] for i in keep_rows])


def main(args):
    outdir = os.path.dirname(args.pca_files[0])

    # pca summary
    pca_full = pd.DataFrame()
    for pca_file in args.pca_files:
        pca_full = pd.concat([pca_full, pd.read_csv(pca_file, sep='\t')])
    pca_full = pca_full.dropna(subset=['PC1', 'PC2'])
    if len(pca_full['sample_name'].unique()) >= 60:
        fig = px.scatter(pca_full, x='PC1', y='PC2', color='peak_caller')
    else:
        fig = px.scatter(pca_full, x='PC1', y='PC2', color='peak_caller', symbol="sample_name")
    fig.write_image(os.path.join(outdir, 'jaccard_summary_pca.pdf'))

    # heatmap summary
    hm_full = pd.DataFrame()
    for hm_file in args.heatmap_files:
        hm_full = pd.concat([hm_full, pd.read_csv(hm_file, sep='\t')])

    # make color map for peak callers
    peakcallers = hm_full['peakcaller'].tolist()
    hm_full = hm_full.drop(columns='peakcaller')

    # only plot the samples that have a score for every comparison
    (hm_full, peakcallers) = drop_na_samples(hm_full, peakcallers)
    if hm_full.shape[1] < 2:
        write_placeholder_plot(
            os.path.join(outdir, 'jaccard_summary_heatmap.pdf'),
            "Jaccard Heatmap Summary",
            "Insufficient valid (non-NA) samples for heatmap clustering."
        )
        return

    colors = sns.color_palette("husl", len(set(peakcallers)))
    color_map = dict(zip(list(set(peakcallers)), colors))

    # set up col indexes for cluster map
    col_labels = hm_full.columns.tolist() * len(set(peakcallers))
    hm_full.index = col_labels

    # heatmap plot

    g = sns.clustermap(hm_full,
                       cmap="YlGnBu", 
                       figsize=(8.5, 11), 
                       col_cluster=False, 
                       row_colors=[color_map[label] for label in peakcallers])
    
    # Add first legend
    legend_elements = []
    for _label, _color in color_map.items():
        legend_elements.append(mpatches.Patch(color=_color, label=_label))
    legend1 = g.ax_heatmap.legend(handles=legend_elements, 
                             title='Peak callers',
                             bbox_to_anchor=(0.3, 1.2),
                             loc='upper left')

    plt.savefig(os.path.join(outdir, 'jaccard_summary_heatmap.pdf'))
    return
    

if __name__ == "__main__":
    desc = \
    dedent("""
    A script to combine multiple PCA coordinates from different
    peak callers, and form one colorized and symboled (if number 
    of samples < 60) PCA plot.
    """)

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        '--pca',
        dest='pca_files',
        nargs="+",
        help="A space delimited list of pca coordinate files for jaccard summary analysis"
    )
    parser.add_argument(
        '--hm',
        nargs="+",
        dest='heatmap_files',
        help="A space delimited list of heatmap coordinate files for jaccard summary analysis"
    )
    main(parser.parse_args())