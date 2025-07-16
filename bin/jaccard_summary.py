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


def main(args):
    outdir = os.path.dirname(args.pca_files[0])

    # pca summary
    pca_full = pd.DataFrame()
    for pca_file in args.pca_files:
        pca_full = pd.concat([pca_full, pd.read_csv(pca_file, sep='\t')])
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