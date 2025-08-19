#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
import csv
import os, sys, re, json
import mimetypes

# Local imports
from utils import Colors, err, fatal


def print_tsv_highlighted(tsv_path):
    """
    Prints TSV data with spaces and tabs highlighted.
    """
    tsv_data = open(tsv_path).read()
    for line in tsv_data.splitlines():
        highlighted_line = ""
        for char in line:
            if char == ' ':
                highlighted_line += Colors.cyan + '•' + Colors.end # Highlight spaces with a middle dot
            elif char == '\t':
                highlighted_line += Colors.green + '→' + Colors.end # Highlight tabs with an arrow
            else:
                highlighted_line += char
        print(highlighted_line)

def check_for_spaces_in_tsv(filepath):
    """
    Checks if any field in a TSV file contains a space character.

    Args:
        filepath (str): The path to the TSV file.

    Returns:
        bool: True if any space is found within a field, False otherwise.
    """
    spaces = []
    spaces_exist = False
    with open(filepath, 'r', newline='', encoding='utf-8') as tsvfile:
        tsv_reader = csv.reader(tsvfile, delimiter='\t')
        for row_num, row in enumerate(tsv_reader, 1):
            for col_num, field in enumerate(row, 1):
                if ' ' in field:
                    spaces.append((row_num, col_num, field))
                    spaces_exist = True
    return spaces_exist, spaces


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


def peakcalls(file, delim="\t"):
    """
    Reads and parses a sample sheet, peakcall.tsv, into a dictionary.
    This file acts as a sample sheet to gather sample metadata and define
    relationship between groups of samples. This file is used to pair a
    ATAC/CHIP/cfCHIP/cutnrun sample with its input sample. This tab 
    delimited file contains two required columns and two optional, all 
    columns are case __in__sensitive.

        - Required: 
            - Sample
            - Group
        - Optionial
            - Blocks
            - InputControl

    It is worth noting that a sample can belong to more than one group. 
    A 1:M sample to group relationship can be denoted by seperating muliptle 
    groups with commas (i.e. ','). This group information is used downstream 
    in the pipeline for DBA.

    Comparisons between groups can be made with a constrast.tsv file.

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

    @return pairs <dict[str or None]>:
        Dictionary containing ChIP-input pairs, where each key is ChIP
        sample and its value is its matched input sample
    @return groups <dict[str]>:
        Dictionary containing group to samples, where each key is group
        and its value is a list of samples belonging to that group
    @return block <dict[str or None]>:
        Dictionary containing samples to blocking information, where each
        key is a sample and each value is blocking information for building
        a linear model
    """
    SAMPLE_COL = 'Sample'.lower()
    INPUT_COL = 'InputControl'.lower()
    GROUP_COL = 'Group'.lower()
    BLOCK_COL = 'Block'.lower()
    tolowerlist = lambda _list: [str(elem).lower() for elem in _list]
    
    spaces_check = check_for_spaces_in_tsv(file)
    if spaces_check[0]:
        print("Spaces detected within peakcall TSV file:\n")
        for row, col, field in spaces_check[1]:
            print(f'\tSpace at row: {row}, column: {col}, field: {field}')
        print('')
        print_tsv_highlighted(file)
        print('')
        raise ValueError('Spaces exist in peakcall file')

    with open(file) as fo:
        rdr = csv.DictReader(fo, delimiter=delim)

        supported_column_names = {"Sample", "InputControl", "Group", "Block"}
        unrecognized_column_names =  set(rdr.fieldnames) - supported_column_names
        if unrecognized_column_names:
            print('')
            print_tsv_highlighted(file)
            print('')
            print("Error: The provided peakcall file contains contains the following unsupported column names: {0}".format(unrecognized_column_names))
            print("Please update the header of your peakcall file. Here is a list of valid column name: {0}".format(supported_column_names))
            raise ValueError('Peakcall file has unsupported headers!')

        rdr.fieldnames = tolowerlist(rdr.fieldnames)
        inputs_exist = INPUT_COL in rdr.fieldnames
        blocks_exist = BLOCK_COL in rdr.fieldnames
        dont_exist = []
        for col in (SAMPLE_COL, GROUP_COL):
            if col not in rdr.fieldnames:
                dont_exist.append(col)

        if dont_exist:
            _c = ', '.join(dont_exist)
            print('')
            print_tsv_highlighted(file)
            print('')
            raise ValueError(f'Peakcall file missing columns {_c}!')
        
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
                pairs[row[SAMPLE_COL]] = ''
            else:
                pairs[row[SAMPLE_COL]] = row[INPUT_COL]
            if ',' in row[GROUP_COL]:
                row[GROUP_COL] = row[GROUP_COL].split(',')
            else:
                row[GROUP_COL] = [row[GROUP_COL]]
            for grp in row[GROUP_COL]:
                groups[grp].append(row[SAMPLE_COL])
            if blocks_exist:
                if row[BLOCK_COL] == '':
                    block[row[SAMPLE_COL]] = ''
                else :
                    block[row[SAMPLE_COL]] = row[BLOCK_COL]
            else:
                block[row[SAMPLE_COL]] = ''

    # group name requirements validation
    bad_labels = []
    ## character black list
    bad_chars = ('*', '_', '-')
    for grp_label in groups.keys():
        for _char in bad_chars:
            if _char in grp_label:
                bad_labels.append(grp_label)

    if bad_labels:
        if len(bad_labels) > 1:
            bad_labels = ', '.join(bad_labels)
        else:
            bad_labels = bad_labels[0]
        raise ValueError('Group(s): ' + bad_labels + '; contain the invalid characters *, -, and/or _ replace and resubmit pipeline')

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
