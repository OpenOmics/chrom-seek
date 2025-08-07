#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def ceildiv(a, b):
    """
    ceiling division, opposite of floor division (//)
    used for calculating the number of rows/columns needed in a grid layout
    """
    return -(a // -b)


def get_command_args():
    """Parse command line arguments for FRiP analysis"""
    parser = argparse.ArgumentParser(
        description="Generate FRiP (Fraction of Reads in Peaks) analysis plots and tables from a single input table",
        epilog="Example: python FRIP_plot.py -t FRiP_data.txt -c config.json -b barplot.png --hm heatmap.png -s scatter.png",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "-t", "--table",
        type=str,
        required=True,
        nargs='+',
        help="Path to FRiP table file",
        metavar="file.txt"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Print verbose output"
    )
    
    parser.add_argument(
        "-b", "--barplot",
        type=str,
        required=True,
        help="Bar plot output file",
        metavar="barplot"
    )
    
    parser.add_argument(
        "--hm",
        type=str,
        required=True,
        help="Heatmap output file",
        metavar="heatmap"
    )
    
    # parser.add_argument(
    #     "-s", "--scatter",
    #     type=str,
    #     required=True,
    #     help="Scatter output file",
    #     metavar="scatter"
    # )
    
    parser.add_argument(
        "-c", "--config",
        type=str,
        required=True,
        help="Config input file",
        metavar="config"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate file existence
    for tbl in args.table:
        if not os.path.exists(tbl):
            print(f"Error: Table file does not exist: {tbl}", file=sys.stderr)
            sys.exit(1)
    
    if not os.path.exists(args.config):
        print(f"Error: Config file does not exist: {args.config}", file=sys.stderr)
        sys.exit(1)
    
    return args

def main(args):
    """Main function to execute FRiP analysis"""
    if args.verbose:
        print(f"Processing FRiP table files: {args.table}")
        print(f"Output files - Bar plot: {args.barplot}, Heatmap: {args.hm}")

    full_table = pd.DataFrame()
    for tbl in args.table:
        if args.verbose:
            print(f"Loading table: {tbl}")
        this_data = pd.read_csv(tbl, sep="\t", header=0)
        full_table = pd.concat([full_table, this_data], ignore_index=True)
    
    #### heatmap ####
    # > make data frame of frip scores with bam id on rows and bed id on columns
    heatmap_data = full_table.pivot_table(index='bamsample', columns='bedsample', values='FRiP', aggfunc='mean')
    # > figure
    sns.clustermap(heatmap_data, cmap='viridis', figsize=(10, 10), annot=False, cbar_kws={'label': 'FRiP Score'})
    # > save figure
    plt.savefig(args.hm, dpi=300)

    #### barplot ####
    # > convert inches to pixels (using 96 pixels per inch as a common standard)
    pixels_per_inch = 96
    width_pixels = 10 * pixels_per_inch
    height_pixels = 10 * pixels_per_inch
    # > figure
    all_peak_tools = full_table['bedtool'].value_counts().to_dict()
    n_peaktools = len(all_peak_tools)
    n_rows = 2 # default number of rows for barplot
    n_columns = ceildiv(n_peaktools, n_rows) 
    
    # > create subplots, scalable in case number of peak tools changes
    subplot_titles = list(all_peak_tools.keys())
    if len(subplot_titles) < (n_rows * n_columns):
        subplot_titles += [''] * ((n_rows * n_columns) - len(subplot_titles))
    fig = make_subplots(rows=n_rows, cols=n_columns, subplot_titles=subplot_titles)
    # > create bar plot for each peak tool
    aligned_tbl = full_table[full_table['bamsample'] == full_table['bedsample']]
    for peak_tool, data in aligned_tbl.groupby('bedtool'):
        if args.verbose:
            print(f"Bar plot - Processing peak tool: {peak_tool}")
        fig.add_trace(
            go.Bar(
                x=data['bamsample'],
                y=data['FRiP'],
                name=peak_tool
            ),
            row=(list(all_peak_tools.keys()).index(peak_tool) // n_columns) + 1,
            col=(list(all_peak_tools.keys()).index(peak_tool) % n_columns) + 1
        )
        fig.update_yaxes(
            title_text="Fraction of Reads in Peaks (FRiP)", 
            row=(list(all_peak_tools.keys()).index(peak_tool) // n_columns) + 1, 
            col=(list(all_peak_tools.keys()).index(peak_tool) % n_columns) + 1
        )
    fig.update_layout(showlegend=False)
    # > save figure
    fig.write_image(args.barplot, width=width_pixels, height=height_pixels, scale=1)

    #### done ####
    print("FRiP analysis completed successfully.")

# Example usage
if __name__ == "__main__":
    args = get_command_args()
    main(args)