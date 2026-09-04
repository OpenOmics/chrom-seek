#!/usr/bin/env python
import argparse
import json
from csv import DictWriter
from collections import defaultdict
from os.path import basename, dirname, exists, isfile, abspath
from os import makedirs

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


def nan_like(value):
    if value is None:
        return True
    v = str(value).strip().lower()
    return v in {"", "na", "nan", "none", "null"}


def peak_has_intervals(path):
    if nan_like(path) or (not isfile(path)):
        return False
    with open(path, "r") as peak_file:
        for line in peak_file:
            line = line.strip()
            if line and not line.startswith("#") and not line.startswith("track") and not line.startswith("browser"):
                return True
    return False


def extract_sid_from_peak(peak_path):
    peak_name = basename(peak_path)
    suffixes = (
        "_peaks.narrowPeak",
        "_peaks.broadPeak",
        ".narrowPeak",
        ".broadPeak",
        ".stringent.bed",
        ".bed",
    )
    for suffix in suffixes:
        if peak_name.endswith(suffix):
            return peak_name[: -len(suffix)]
    return peak_name.split(".")[0]


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
    control_map = args.cfg['project']['peaks']['inputs']
    samples = args.cfg['project']['peaks']['chips']
    grp2sample = args.cfg['project']['groups']
    sample2grp = dict.fromkeys(samples, "")
    bam_map = {extract_sid(b): b for b in args.bams}
    peak_map = {extract_sid_from_peak(p): p for p in args.peaks}
    for grp, grp_sample in grp2sample.items():
        for s in grp_sample:
            sample2grp[s] = grp

    dropped = []
    condition_replicates = defaultdict(int)
    csv = []
    for sample in samples:
        condition = sample2grp.get(sample, "")
        bam_reads = bam_map.get(sample, "")
        control_id = control_map.get(sample, "")
        peaks = peak_map.get(sample, "")

        # Drop NaN/missing/empty peaks rows and empty BAM rows.
        if nan_like(sample) or nan_like(bam_reads) or nan_like(peaks) or (not peak_has_intervals(peaks)):
            dropped.append(sample)
            continue

        condition_replicates[condition] += 1
        csv.append({
            "SampleID": sample,
            "Condition": condition,
            "Replicate": str(condition_replicates[condition]),
            "bamReads": bam_reads,
            "ControlID": "" if nan_like(control_id) else control_id,
            "bamControl": "" if nan_like(control_id) else bam_map.get(control_id, ""),
            "Peaks": peaks,
            "PeakCaller": args.pktool,
        })

    if dropped:
        print("\t> WARNING: Dropped samples with NaN/empty peak inputs: " + ", ".join(dropped))

    out_dir = dirname(args.output)
    if not exists(out_dir):
        makedirs(out_dir, 0o755, exist_ok=True)

    with open(args.output, "w") as csv_out:
        wrtr = DictWriter(csv_out, columns, delimiter=",")
        wrtr.writeheader()
        for row in csv:
            wrtr.writerow(row)
    print(f"\t> File {args.output} written with {len(csv)} sample(s).")


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
