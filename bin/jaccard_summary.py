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
from textwrap import dedent


def main(args):
    pca_full = pd.DataFrame()
    outdir = os.path.dirname(args.pca_files[0])
    for pca_file in args.pca_files:
        pca_full = pd.concat([pca_full, pd.read_csv(pca_file, sep='\t')])
    if len(pca_full['sample_name'].unique()) >= 60:
        fig = px.scatter(pca_full, x='PC1', y='PC2', color='peak_caller')
    else:
        fig = px.scatter(pca_full, x='PC1', y='PC2', color='peak_caller', symbol="sample_name")
    fig.write_image(os.path.join(outdir, 'jaccard_summary_pca.pdf'))
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
        'pca_files',
        nargs="+",
        help="A space delimited list of pca coordinate files for jaccard summary analysis"
    )
    main(parser.parse_args())