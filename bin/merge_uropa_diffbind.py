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


def fuzz_merge(db, uropa):
    uropa_small = uropa[["chr", "start", "end"]]
    db_small = db[["chr", "start", "end"]]
    db_uropa_map = {}

    for contig in db_small["chr"].unique():
        this_db = db_small[db_small["chr"] == contig]
        this_uropa = uropa_small[uropa_small["chr"] == contig]
        this_uropa['distance'] = 0
        for i, row in this_db.iterrows():
            check_uropa = this_uropa
            check_uropa['distance'] = check_uropa['start'].apply(lambda x: abs(row['start'] - x))
            check_uropa = check_uropa.dropna(subset=['distance'])
            check_uropa = check_uropa[check_uropa['distance'] <= 2]
            if check_uropa.shape[0] < 1:
                continue
            else:
                db_uropa_map[i] = check_uropa.index.tolist()

    merge_df = pd.DataFrame()
    for db_i, ur_i in db_uropa_map.items():
        this_db_row = db.loc[db_i].to_dict()
        this_ur_row = uropa.loc[ur_i].fillna('NA').to_dict(orient='records')
        
        for ur_row in this_ur_row:
            del ur_row['chr']
            del ur_row['start']
            del ur_row['end']
            newrow = this_db_row
            newrow.update(ur_row)
            merge_df = pd.concat([merge_df, pd.DataFrame([newrow])], ignore_index=True)
    return merge_df


def main(args):
    diffbind = pd.read_csv(args.diffbind, sep="\t")
    diffbind = diffbind.rename(columns={"seqnames": "chr"})  # start & end exist
    diffbind['chr'] = diffbind['chr'].astype(str)
    # count filter
    fold_filter = diffbind[diffbind["Fold"].abs() < args.fold]
    n_filter_fold = fold_filter.shape[0]
    fdr_filter = diffbind[diffbind["FDR"] > args.fdr]
    n_filter_fdr = fdr_filter.shape[0]
    # filter
    diffbind = diffbind[
        (diffbind["Fold"].abs() >= args.fold) & (diffbind["FDR"] <= args.fdr)
    ]
    diffbind.reset_index(drop=True, inplace=True)
    n_after_filter = diffbind.shape[0]
    # describe filter
    print(f"-- {str(n_filter_fdr)} peaks filtered for being > {str(args.fdr)} FDR --\n")
    print(str(fold_filter) + "\n\n")
    print(f"-- {str(n_filter_fold)} peaks filtered for being < abs({str(args.fold)}) fold-change --\n")
    print(str(fdr_filter) + "\n\n")
    print(f"-- {str(n_after_filter)} total peaks removed by FDR and fold-change filters")
    # format
    diffbind = diffbind.drop(columns=["strand", "Conc"])
    diffbind = diffbind.reset_index(drop=True)
    uropa = pd.read_csv(args.uropa, sep="\t")
    uropa = uropa.rename(
        columns={"peak_chr": "chr", "peak_start": "start", "peak_end": "end"}
    )
    uropa['chr'] = uropa['chr'].astype(str)
    uropa = uropa.drop(columns=["peak_score", "peak_strand"])
    uropa = uropa.reset_index(drop=True)
    # join tables
    n_uropa = uropa.shape[0]
    n_diffbind = diffbind.shape[0]
    merged = fuzz_merge(diffbind, uropa)
    merged = merged.reset_index(drop=True)
    n_merged = merged.shape[0]
    print(f"-- {n_uropa} peaks from uropa annotation, {n_diffbind} peaks from diffbind consensus -- \n\n")
    print(f"-- {n_merged} peaks merged from uropa + diffbind sources -- \n\n")
    merged = merged.sort_values(by='FDR', ascending=True)
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
        help="Finalhits input file form `UROPA_diffbind`",
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
