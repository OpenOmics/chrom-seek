#!/usr/bin/env python
import pandas as pd
import argparse
import os


"""
Script for joining diffbind peak set to uropa gene annotations

Caveat: this won't work with the uropa allhits file because of the 1:Many
        relationship between peaks and genes.

Example uropa entry

    peak_chr                          chr1
    peak_start                      713768
    peak_end                        714581
    peak_id                          Peak1
    peak_score                           0
    peak_strand                          .
    feature                           gene
    feat_start                    621058.0
    feat_end                      622053.0
    feat_strand                          -
    feat_anchor                      start
    distance                       92121.0
    relative_location             Upstream
    feat_ovl_peak                      0.0
    peak_ovl_feat                      0.0
    gene_id              ENSG00000185097.2
    gene_name                       OR4F16
    gene_type               protein_coding
    name                           query_3

Example diffbind entry

    seqnames           chr1
    start            713769
    end              714581
    width               813
    strand                *
    Conc           5.864859
    Conc_IFN0h     5.622903
    Conc_IFN24h    6.072004
    Fold          -0.449101
    p.value        0.398368
    FDR            0.633526
    Called1               1
    Called2               2

Joined entry

    chr                               chr1
    start                           713768
    end                             714581
    peak_id                          Peak1
    feature                           gene
    feat_start                    621058.0
    feat_end                      622053.0
    feat_strand                          -
    feat_anchor                      start
    distance                       92121.0
    relative_location             Upstream
    feat_ovl_peak                      0.0
    peak_ovl_feat                      0.0
    gene_id              ENSG00000185097.2
    gene_name                       OR4F16
    gene_type               protein_coding
    name                           query_3
    width                              813
    Conc_IFN0h                    5.622903
    Conc_IFN24h                   6.072004
    Fold                         -0.449101
    p.value                       0.398368
    FDR                           0.633526
    Called1                              1
    Called2                              2

"""


def valid_path(path_str):
    if not os.path.exists(path_str):
        raise argparse.ArgumentTypeError(f"Path '{path_str}' does not exist.")
    return path_str


def main(args):
    diffbind = pd.read_csv(args.diffbind, sep="\t")
    diffbind = diffbind.rename(columns={"seqnames": "chr"})  # start & end exist
    diffbind = diffbind[
        (diffbind["Fold"] >= args.fold) & (diffbind["FDR"] <= args.fdr)
    ]  # filter
    diffbind = diffbind.drop(columns=["strand", "Conc"])
    diffbind = diffbind.reset_index(drop=True)

    uropa = pd.read_csv(args.uropa, sep="\t")
    uropa = uropa.rename(
        columns={"peak_chr": "chr", "peak_start": "start", "peak_end": "end"}
    )
    uropa = uropa.drop(columns=["peak_score", "peak_strand"])
    uropa = uropa.reset_index(drop=True)

    merged = uropa.merge(diffbind, on=["chr", "start", "end"])
    merged = merged.reset_index(drop=True)
    merged.to_csv(args.output, sep="\t", index=False)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to join uropa and diffbind outputs"
    )
    parser.add_argument(
        "--diffbind",
        "-d",
        type=valid_path,
        required=True,
        help="CSV input file from `diffbind_prep`",
    )
    parser.add_argument(
        "--uropa",
        "-u",
        type=valid_path,
        required=True,
        help="Allhits input file form `UROPA_diffbind`",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Location to save joined *_finalhits.txt and DiffBind_prep.csv",
    )
    parser.add_argument(
        "--fdr",
        "-f",
        type=float,
        help="FDR cutoff for filtering",
        default=0.05,
    )
    parser.add_argument(
        "--fold",
        "-l",
        type=float,
        help="Fold change cutoff for filtering",
        default=0,
    )
    args = parser.parse_args()
    main(args)
