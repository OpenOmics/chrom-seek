#!/usr/bin/env python
import argparse
from csv import DictWriter
from os.path import basename, dirname, exists
from os import makedirs
from itertools import repeat

##
## Objective : gather all Q5DD bams, their respective controls (if they exist),
##              and their peaksets together in diffbind-esque csv
##              see : https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf


def main(args):
    extract_sid = lambda fn: basename(fn).replace(".Q5DD.bam", "")
    n = len(args.sample)
    columns = [
        "SampleID",
        "Condition",
        "Replicate",
        "bamReads",
        "ControlID",
        "bamControl",
        "Peaks",
        "PeakCaller",
    ]
    tbl = {}
    tbl["SampleID"] = list(map(extract_sid, args.sample))
    tbl["Condition"] = list(repeat("", n))
    tbl["Replicate"] = list(repeat("1", n))
    tbl["bamReads"] = args.sample
    if args.control:
        tbl["ControlID"] = list(map(extract_sid, args.control))
        tbl["bamControl"] = args.control
    else:
        tbl["ControlID"] = list(repeat("", n))
        tbl["bamControl"] = list(repeat("", n))
    tbl["Peaks"] = args.peaks
    tbl["PeakCaller"] = list(repeat(args.pktool, n))
    csv = []
    for i in range(n):
        this_row = {}
        for col in columns:
            this_row[col] = tbl[col][i]
        csv.append(this_row)

    out_dir = dirname(args.output)
    if not exists(out_dir):
        makedirs(out_dir, 0o755, exist_ok=True)

    with open(args.output, "w") as csv_out:
        wrtr = DictWriter(csv_out, columns, delimiter=",")
        wrtr.writeheader()
        for row in csv:
            wrtr.writerow(row)
    print(f"\t> File {args.output} written.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script for creating diffbind CSVs for QC purposes"
    )
    parser.add_argument(
        "-t",
        "--tool",
        dest="pktool",
        type=str,
        help="Single string, identify peak tool to diffbind\n"
        + "(see: https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        help="Path to output AllSamples-* csvfile",
    )
    parser.add_argument(
        "-s",
        "--samplebams",
        dest="sample",
        nargs="+",
        help="List of the sample BAM files",
    )
    parser.add_argument(
        "-c",
        "--controlbams",
        dest="control",
        nargs="?",
        default=None,
        help="List of the control BAM files",
    )
    parser.add_argument(
        "-p", "--peaks", dest="peaks", nargs="+", help="List of sample PEAKSETs"
    )
    main(parser.parse_args())
