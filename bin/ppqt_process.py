#!/usr/bin/env python3

# Purpose:
#   To grab the estimated fragment length from the ppqt output and a
#   small txt with that information. For input files, adding an extra
#   value of 200bp as an alternative.
import argparse


def main(args):
    cln_file = lambda z: z.strip().split()
    this_file = list(map(cln_file, open(args.ppqt, "r").readlines()))
    ppqt_values = this_file[0][2].split(",")
    frag_len = None
    for ppqt_value in ppqt_values:
        if int(ppqt_value) > 150:
            frag_len = ppqt_value
            break

    if frag_len is None:
        frag_len = 200
    print(frag_len)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to extract the the estimated fragment length from the ppqt output."
    )
    parser.add_argument("ppqt", help="Name of the ppqt txt file")
    args = parser.parse_args()
    main(args)
