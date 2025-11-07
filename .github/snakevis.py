#!/usr/bin/env python
import json
import os
import argparse
import subprocess
import re
from glob import glob
from configparser import ConfigParser


REPO_ROOT = os.path.abspath(os.path.join(os.path.basename(__file__), '..'))


def join_jsons(templates):
    """Joins multiple JSON files to into one data structure
    Used to join multiple template JSON files to create a global config dictionary.
    @params templates <list[str]>:
        List of template JSON files to join together
    @return aggregated <dict>:
        Dictionary containing the contents of all the input JSON files
    """
    # Get absolute PATH to templates in git repo
    repo_path = os.path.dirname(os.path.abspath(__file__))
    aggregated = {}

    for file in templates:
        with open(os.path.join(repo_path, file), 'r') as fh:
            aggregated.update(json.load(fh))

    return aggregated


def get_nested_value(data, keys):
    """
    Retrieves a value from a nested dictionary using a list of keys.

    Args:
        data (dict): The nested dictionary.
        keys (list): A list of keys representing the path to the desired value.

    Returns:
        The value at the specified path, or None if any key in the path is not found.
    """
    current_level = data
    for key in keys:
        if isinstance(current_level, dict) and key in current_level:
            current_level = current_level[key]
        else:
            return None  # Key not found at this level
    return current_level


def add_sample_metadata(input_files, config):
    """Adds sample metadata such as sample basename, label, and group information.
    If sample sheet is provided, it will default to using information in that file.
    If no sample sheet is provided, it will only add sample basenames and labels.
    @params input_files list[<str>]:
        List containing pipeline input fastq files
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @params group <str>:
        Sample sheet containing basename, group, and label for each sample
    @return config <dict>:
        Updated config with basenames, labels, and groups (if provided)
    """
    config['samples'] = []
    for file in input_files:
        # Split sample name on file extension
        sample = re.split('(?i)\\.(?:sorted\\.)?(?:bam|sam|cram)', os.path.basename(file))[0]
        if sample not in config['samples']:
            config['samples'].append(sample)

    return config


def main(args):
    if not os.path.exists(args.config):
        print('Config file not present -- not running snakevis')
    cfg = ConfigParser()
    cfg.read(args.config)
    this_cfg = cfg['SNAKEVIS']

    test_files = glob(os.path.join(REPO_ROOT, this_cfg['test_files']))
    snk_file = os.path.join(REPO_ROOT, this_cfg['snk_file'])
    jsons_to_join = [this_cfg['config_include']] if ',' not in this_cfg['config_include'] else json.loads(this_cfg['config_include'])
    jsons_to_join = [os.path.join(REPO_ROOT, this_json) for this_json in jsons_to_join]
    config_output_fn = os.path.join(REPO_ROOT, this_cfg['config_out'])
    output_dir = os.path.dirname(config_output_fn)
    
    pipeline_config = join_jsons(jsons_to_join)
    pipeline_config = add_sample_metadata(test_files, pipeline_config)

    if '.' in this_cfg['config_sample_paths_key']:
        cfg_levels = this_cfg['config_sample_paths_key'].split('.')
        to_change = None
        for i, level in enumerate(cfg_levels[:-1]):
            if i == 0:
                to_change = pipeline_config[level]
            else:
                to_change = to_change[level]
    else:
        to_change = pipeline_config
    
    to_change[this_cfg['config_sample_paths_key'].split('.')[-1]] = test_files

    if 'config_cluster' in this_cfg:
        pipeline_config['cluster'] = json.load(open(os.path.join(REPO_ROOT, this_cfg['config_cluster'])))

    if any('CONFIG' in k for k in cfg.keys()):
        cfg_keys = [k for k in cfg.keys() if 'CONFIG' in k]
        for k in cfg_keys:
            levels = [x for x in k.split('.') if x != 'CONFIG']
            to_change = None
            if levels:
                for i, level in enumerate(levels):
                    if i == 0:
                        if level not in pipeline_config:
                            pipeline_config[level] = {}
                        to_change = pipeline_config[level]
                    else:
                        if level not in to_change:
                            to_change[level] = {}
                        to_change = to_change[level]
            else:
                to_change = pipeline_config
            
            for k, v in cfg[k].items():
                to_change[k] = v

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    with open(config_output_fn, "w") as f:
        json.dump(pipeline_config, f)

    dot_out = os.path.join(output_dir, 'pipeline_rulegraph.dot')
    snakemake_cmd = f'snakemake --configfile={config_output_fn} -s {snk_file} -d {output_dir} --forceall --rulegraph > {dot_out}'
    os.system(snakemake_cmd)

    if not os.path.exists(args.upload_path):
        os.makedirs(args.upload_path, exist_ok=True)
    rulegraph_out = os.path.join(args.upload_path, 'pipeline_rulegraph.svg')
    if 'rule_exclude' in this_cfg:
        exclude_rules = ' '.join(json.loads(this_cfg['rule_exclude']))
        snakevis_cmd = f'snakevision -s {exclude_rules} -o {rulegraph_out} {dot_out}'
    else:
        snakevis_cmd = f'snakevision -o {rulegraph_out} {dot_out}'
    os.system(snakevis_cmd)

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="ARG PARSE DESCRIPTION"
    )
    parser.add_argument(
        "config",
        type=str,
        default=None,
        help="Valid CFG file for snakevis running"
    )
    parser.add_argument(
        "upload_path",
        type=str,
        default=None,
        help="Path to artifact upload directory"
    )
    args = parser.parse_args()
    main(args)
