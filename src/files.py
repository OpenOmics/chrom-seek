#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
import csv
import os, sys, re, json
import mimetypes

# Local imports
from utils import Colors, err, fatal


def clean(s, remove=['"', "'"]):
    """Cleans a string to remove any defined leading or trailing characters.
    @param s <str>:
        String to clean.
    @param remove list[<str>]:
        List of characters to remove from beginning or end of string 's'.
    @return s <str>:
        Cleaned string
    """
    for c in remove:
        s = s.strip(c)
    return s


def index(
    file, delim="\t", required=["sample", "group"], optional=["inputcontrol", "block"]
):
    """Return the index of expected columns in provided file. If an optional
    column is not provided, then None is returned its column index. The peakcall
    file is expected to have the following required and optional columns.
    @Required columns:
        - chip, group
    @Optional columns:
        - block, inputcontrol
    If an optional column is not provided, then None is returned its index.
    @param file <str>:
        Path to peakcall TSV file.
    @return tuple(indices <dict[int/None]>, hasHeader <boolean>):
        [0] Dictionary containing information the index of each required/optional column
        [1] Boolean to indicate whether file has a header
    """
    c = Colors()
    indices = {}
    has_header = True

    # Check to see if the file is empty
    fh = open(file, "r")
    try:
        header = [clean(col.lower().strip()) for col in next(fh).strip().split(delim)]
    except StopIteration:
        err(
            "{}{}Error: peakcall file, {}, is empty!{}".format(
                c.bg_red, c.white, file, c.end
            )
        )
        fatal(
            "{}{}Please add ChIP-Input pairs and group information to the file and try again.{}".format(
                c.bg_red, c.white, c.end
            )
        )
    finally:
        fh.close()

    # Parse the header to get the index of required fields
    try:
        # Get index of ChIP, Input, Group
        # columns for parsing the file
        for col in required:
            indices[col] = header.index(col)
    except ValueError:
        # Missing column names or header in peakcall file
        # This can also occur if the file is not actually
        # a tab delimited file.
        # TODO: Add a check to see if the file is actually
        # a tab delimited file, i.e. a TSV file.
        has_header = False
        err(
            "{}{}Warning: {} is missing at least one of the following column names: Sample, InputControl, Group {}".format(
                c.bg_yellow, c.black, file, c.end
            )
        )
        err(
            "{}{}\t  └── Making assumptions about columns in the peakcall file... 1=Sample, 2=InputControl, 3=Group, 4=Block {}".format(
                c.bg_yellow, c.black, c.end
            )
        )
        # Setting column indexes to the following defaults:
        # 0 = ChIP sample column
        # 1 = Input sample column
        # 2 = Group information column
        for i in range(len(required)):
            indices[required[i]] = i

    # Parse the header to get the index of optional fields
    # Get index of optional column, default to None if DNE
    for i, col in enumerate(optional):
        if has_header:
            # Contains header, if DNE set to None
            try:
                indices[col] = header.index(col)
            except ValueError:
                indices[col] = None
        else:
            # No header, check if in expected column
            # 0 = ChIP sample column
            # 1 = Input sample column
            # 2 = Group information column
            # 3 = Block information column <- Optional
            try:
                offset = len(required)
                default_index = offset + i
                tmp = header[default_index]
                indices[col] = default_index
            except IndexError:
                # DNE, set to None
                indices[col] = None

    return indices, has_header


def peakcalls(file, delim="\t"):
    """
    Reads and parses a sample sheet, peakcall.tsv, into a dictionary.
    This file acts as a sample sheet to gather sample metadata and define
    relationship between groups of samples. This file is used to pair a
    ChIP sample with its input sample. This tab delimited file contains
    three columns. One column for the basename of the ChIP sample, one
    column for the basename of the paired sample, and lastly, one column
    for the name of the sample's group. It is worth noting that a sample
    can belong to more than one group. A 1:M sample to group relationship
    can be denoted by seperating muliptle groups with commas (i.e. ',').
    This group information is used downstream in the pipeline for DBA.
    Comparisons between groups can be made with a constrast.tsv file.
    This function returns a tuple containing the ChIP-input dictionary
    and a second dictionary containing group to sample lists.

    @Example: peakcall.tsv
        Sample    InputControl   Group   Block
        cfChIP_001	Input_001	G1,G4   B1
        cfChIP_002	Input_002	G1,G4   B1
        cfChIP_003	Input_003	G1,G4   B2
        cfChIP_004	Input_004	G2,G5   B2
        cfChIP_005	Input_005	G2,G5   B3
        cfChIP_006	Input_006	G2,G5   B3
        cfChIP_007	Input_007	G3,G5   B4
        cfChIP_008	Input_008	G3,G5   B5
        cfChIP_009	Input_009	G3  B6
        cfChIP_000	Input_000	G3  B7

    >> chip2input, groups = peakcalls('peakcall.tsv')
    >> chip2input
    {
        'cfChIP_002': 'Input_002',
        'cfChIP_003': 'Input_003',
        'cfChIP_000': 'Input_000',
        'cfChIP_001': 'Input_001',
        'cfChIP_006': 'Input_006',
        'cfChIP_007': 'Input_007',
        'cfChIP_004': 'Input_004',
        'cfChIP_005': 'Input_005',
        'cfChIP_008': 'Input_008',
        'cfChIP_009': 'Input_009'
    }
    >> groups
    {
        'G1': ['cfChIP_001', 'cfChIP_002', 'cfChIP_003'],
        'G2': ['cfChIP_004', 'cfChIP_005', 'cfChIP_006'],
        'G3': ['cfChIP_007', 'cfChIP_008', 'cfChIP_009', 'cfChIP_000'],
        'G4': ['cfChIP_001', 'cfChIP_002', 'cfChIP_003'],
        'G5': ['cfChIP_004', 'cfChIP_005', 'cfChIP_006', 'cfChIP_007', 'cfChIP_008'],
    }
    >> blocks
    {
        'cfChIP_002': 'B1',
        'cfChIP_003': 'B1',
        'cfChIP_000': 'B2',
        'cfChIP_001': 'B2',
        'cfChIP_006': 'B3',
        'cfChIP_007': 'B3',
        'cfChIP_004': 'B4',
        'cfChIP_005': 'B5',
        'cfChIP_008': 'B6',
        'cfChIP_009': 'B7'
    }

    @param file <str>:
        Path to peakcall TSV file.

    @return pairs <dict[str]>:
        Dictionary containing ChIP-input pairs, where each key is ChIP
        sample and its value is its matched input sample
    @return groups <dict[str]>:
        Dictionary containing group to samples, where each key is group
        and its value is a list of samples belonging to that group
    @return block <dict[str]>:
        Dictionary containing samples to blocking information, where each
        key is a sample and each value is blocking information for building
        a linear model
    """
    SAMPLE_COL = 'Sample'
    INPUT_COL = 'InputControl'
    GROUP_COL = 'Group'
    BLOCK_COL = 'Block'

    with open(file) as fo:
        rdr = csv.DictReader(fo, delimiter=delim)
        inputs_exist = INPUT_COL in rdr.fieldnames
        blocks_exist = BLOCK_COL in rdr.fieldnames
        dont_exist = []
        for col in (SAMPLE_COL, GROUP_COL):
            if col not in rdr.fieldnames:
                dont_exist.append(col)

        if dont_exist:
            _c = ', '.join(dont_exist)
            raise ValueError(f'peakcall file missing columns {_c}')

        all_groups = []
        for row in rdr:
            if ',' in row[GROUP_COL]:
                all_groups.extend(row[GROUP_COL].split(','))
            else:
                all_groups.append(row[GROUP_COL])
        groups = {k: [] for k in all_groups}

        fo.seek(0); next(rdr) # skip header
        pairs = {}
        block = {}
        for row in rdr:
            if not inputs_exist:
                pairs[row[SAMPLE_COL]] = None
            else:
                pairs[row[SAMPLE_COL]] = row[INPUT_COL]
            if ',' in row[GROUP_COL]:
                row[GROUP_COL] = row[GROUP_COL].split(',')
            else:
                row[GROUP_COL] = [row[GROUP_COL]]
            for grp in row[GROUP_COL]:
                groups[grp].append(row[SAMPLE_COL])
            if blocks_exist:
                block[row[SAMPLE_COL]] = row[BLOCK_COL]
            else:
                block[row[SAMPLE_COL]] = ''

    return pairs, groups, block


def contrasts(file, groups, delim="\t"):
    """Reads and parses the group comparison file, contrasts.tsv, into a
    dictionary. This file acts as a config file to setup contrasts between
    two groups, where groups of samples are defined in the peakcalls.tsv file.
    This information is used in differential analysis, like differential binding
    analysis or differential gene expression, etc.
    @Example: contrasts.tsv
        G2  G1
        G4  G3
        G5  G1
    >> contrasts = contrasts('contrasts.tsv', groups = ['G1', 'G2', 'G3', 'G4', 'G5'])
    >> contrasts
    [
        ["G2",  "G1"],
        ["G4",  "G3"],
        ["G5",  "G1"]
    ]
    @param file <str>:
        Path to contrasts TSV file.
    @param groups list[<str>]:
        List of groups defined in the peakcall file, enforces groups exist.
    @return comparisons <list[list[str, str]]>:
        Nested list contain comparsions of interest.
    """

    c = Colors()
    errors = []
    comparsions = []
    line_number = 0
    with open(file) as fh:
        for line in fh:
            line_number += 1
            linelist = [clean(l.strip()) for l in line.split(delim)]
            try:
                g1 = linelist[0]
                g2 = linelist[1]
                if not g1 or not g2:
                    continue  # skip over empty lines
            except IndexError:
                # Missing a group, need two groups to tango
                # This can happen if the file is NOT a TSV file,
                # and it is seperated by white spaces, :(
                err(
                    "{}{}Warning: {} is missing at least one group on line {}: {}{}".format(
                        c.bg_yellow, c.black, file, line_number, line.strip(), c.end
                    )
                )
                err(
                    "{}{}\t  └── Skipping over line, check if line is tab seperated... {}".format(
                        c.bg_yellow, c.black, c.end
                    )
                )
                continue
            # Check to see if groups where defined already,
            # avoids user errors and spelling errors
            for g in [g1, g2]:
                if g not in groups:
                    # Collect all error and report them at end
                    errors.append(g)

            # Add comparsion to list of comparisons
            if [g1, g2] not in comparsions:
                comparsions.append([g1, g2])

    if errors:
        # One of the groups is not defined in peakcalls
        err(
            '{}{}Error: the following group(s) in "{}" are not defined in peakcall file! {}'.format(
                c.bg_red, c.white, file, c.end
            )
        )
        fatal("{}{}\t  └── {} {}".format(c.bg_red, c.white, ",".join(errors), c.end))

    return comparsions


def validate_custom_genome(genome_json):
    if not os.path.exists(genome_json):
        raise FileNotFoundError(
            f"Custom genome definition {genome_json} does not exist!"
        )
    if not mimetypes.guess_type(genome_json)[0] in ("text/plain", "application/json"):
        raise ValueError(
            f"Custom genome definition {genome_json} is not a plain text json file"
        )
    with open(genome_json, "r") as file:
        try:
            this_genome = json.load(file)
        except json.JSONDecodeError:
            raise ValueError(
                f"JSON syntax is broken in custom genome definition {genome_json}"
            )
    required_keys = [
        "ALIAS",
        "SUPPORTED_PIPELINES",
        "BLACKLISTBWAINDEX",
        "BLACKLISTGENRICH",
        "BWA",
        "cfChIP_TOOLS_SRC",
        "EFFECTIVEGENOMESIZE",
        "GENEINFO",
        "GENOME",
        "GENOMECHR",
        "GTFFILE",
        "REFLEN",
        "FRAC",
        "MEME_VERTEBRATES_DB",
        "MEME_EUKARYOTE_DB",
        "MEME_GENOME_DB",
    ]
    bad_columns = []
    genome_alias = list(this_genome["references"].values())[0]["ALIAS"]
    genome_ks = list(this_genome["references"].values())[0].keys()
    for k in required_keys:
        if k not in genome_ks:
            bad_columns.append(k)

    if bad_columns:
        raise ValueError(
            f"Custom genome definition {genome_json} (alias {genome_alias}) is missing keys: {','.join(bad_columns)}"
        )
    return this_genome


if __name__ == "__main__":
    # Testing peakcall TSV parser
    print("Parsing peakcall file...")
    chip2input, groups, blocking = peakcalls(sys.argv[1])
    print(chip2input)
    print(groups)
    print(blocking)
    print("Parsing contrasts file...")
    comparsions = contrasts(sys.argv[2], groups=groups.keys())
    print(comparsions)
