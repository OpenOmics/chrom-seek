#!/usr/bin/env python
import argparse
import json
from csv import DictWriter
from os.path import basename, dirname, exists, isfile, abspath
from os import makedirs
from itertools import repeat

##
## Objective : gather all Q5DD bams, their respective controls (if they exist),
##              and their peaksets together in diffbind-esque csv
##              see : https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf
##


def valid_json(path):
    path = abspath(path)
    if not isfile(path):
        raise argparse.ArgumentTypeError(f"'{path}' is not a valid file path")
    with open(path, "r") as file:
        data = json.load(file)
        return data


def main(args):
    extract_sid = lambda fn: basename(fn).replace(".Q5DD.bam", "")
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

    control_map = args.cfg['project']['peaks']['inputs']
    samples = args.cfg['project']['peaks']['chips']
    n = len(samples)
    grp2sample = args.cfg['project']['groups']
    sample2grp = dict.fromkeys(samples, "")
    bam_map = {extract_sid(b): b for b in args.bams}
    peak_map = {extract_sid(p): p for p in args.peaks}
    for grp, grp_sample in grp2sample.items():
        for s in grp_sample:
            sample2grp[s] = grp

    tbl["SampleID"] = samples
    tbl["Condition"] = list(map(sample2grp.get, samples))
    tbl["Replicate"] = list(repeat("1", n))
    tbl["bamReads"] = list(map(bam_map.get, map(extract_sid, samples)))
    tbl["ControlID"] = list(map(control_map.get, samples))
    tbl["bamControl"] = list(map(bam_map.get, map(control_map.get, samples)))
    tbl["Peaks"] = list(map(peak_map.get, map(extract_sid, args.peaks)))
    tbl["PeakCaller"] = list(repeat(args.pktool, n))

    csv = []
    for i in range(n):
        this_row = {}
        for col in columns:
            this_row[col] = list(tbl[col])[i]
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
        type=lambda p: abspath(p),
        help="Path to output AllSamples-* csvfile",
    )
    parser.add_argument(
        "-p", "--peaks", dest="peaks", nargs="+", help="List of sample PEAKSETs"
    )
    parser.add_argument(
        "-b", "--bams", dest="bams", nargs="+", help="List of sample and control BAMs"
    )
    parser.add_argument(
        "-c",
        "--config",
        dest="cfg",
        type=valid_json,
        help="JSON pipeline configuration file",
        required=True,
    )
    main(parser.parse_args())
