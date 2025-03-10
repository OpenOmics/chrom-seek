#!/usr/bin/env python
import argparse
from csv import DictWriter
from os.path import basename
from itertools import repeat


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
    tbl["Replicate"] = list(map(str, range(1, n + 1)))
    tbl["bamReads"] = args.sample
    tbl["ControlID"] = list(map(extract_sid, args.control))
    tbl["bamControl"] = args.control
    tbl["Peaks"] = args.peaks
    tbl["PeakCaller"] = list(repeat(args.pktool, n))
    csv = zip(*[tbl[col] for col in columns])
    with open(args.output) as csv_out:
        wrtr = DictWriter(csv_out, columns, delimter=",")
        wrtr.writeheader()
        for row in csv:
            wrtr.writerow(row)
    print(f"\t> File {args.output} written.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-t",
        "--tool",
        dest="pktool",
        type="string",
        help="Single string, identify peak tool to diffbind\n"
        + "(see: https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type="string",
        help="Path to output AllSamples-* csvfile",
    )
    parser.add_argument("-s", "--samplebams", dest="sample", nargs="+", help="")
    parser.add_argument("-c", "--controlbams", dest="control", nargs="+", help="")
    parser.add_argument("-p", "--peaks", dest="peaks", nargs="+", help="")
    main(parser.parse_args())
