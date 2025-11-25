#!/usr/bin/env python
import json
import os
import sys
import argparse
import shutil
import tempfile
import jinja2
import ast
import re
import subprocess


def rename(filename):
    """Dynamically renames FastQ file to have one of the following extensions: *.R1.fastq.gz, *.R2.fastq.gz
    To automatically rename the fastq files, a few assumptions are made. If the extension of the
    FastQ file cannot be infered, an exception is raised telling the user to fix the filename
    of the fastq files.
    @param filename <str>:
        Original name of file to be renamed
    @return filename <str>:
        A renamed FastQ filename
    """
    # Covers common extensions from SF, SRA, EBI, TCGA, and external sequencing providers
    # key = regex to match string and value = how it will be renamed
    extensions = {
        # Matches: _R[12]_fastq.gz, _R[12].fastq.gz, _R[12]_fq.gz, etc.
        ".R1.f(ast)?q.gz$": ".R1.fastq.gz",
        ".R2.f(ast)?q.gz$": ".R2.fastq.gz",
        # Matches: _R[12]_001_fastq_gz, _R[12].001.fastq.gz, _R[12]_001.fq.gz, etc.
        # Capture lane information as named group
        ".R1.(?P<lane>...).f(ast)?q.gz$": ".R1.fastq.gz",
        ".R2.(?P<lane>...).f(ast)?q.gz$": ".R2.fastq.gz",
        # Matches: _[12].fastq.gz, _[12].fq.gz, _[12]_fastq_gz, etc.
        "_1.f(ast)?q.gz$": ".R1.fastq.gz",
        "_2.f(ast)?q.gz$": ".R2.fastq.gz",
    }

    if filename.endswith(".R1.fastq.gz") or filename.endswith(".R2.fastq.gz"):
        # Filename is already in the correct format
        return filename

    converted = False
    for regex, new_ext in extensions.items():
        matched = re.search(regex, filename)
        if matched:
            # regex matches with a pattern in extensions
            converted = True
            filename = re.sub(regex, new_ext, filename)
            break  # only rename once

    if not converted:
        raise NameError(
            """\n\tFatal: Failed to rename provided input '{}'!
        Cannot determine the extension of the user provided input file.
        Please rename the file list above before trying again.
        Here is example of acceptable input file extensions:
          sampleName.R1.fastq.gz      sampleName.R2.fastq.gz
          sampleName_R1_001.fastq.gz  sampleName_R2_001.fastq.gz
          sampleName_1.fastq.gz       sampleName_2.fastq.gz
        Please also check that your input files are gzipped?
        If they are not, please gzip them before proceeding again.
        """.format(filename)
        )

    return filename


def validate_json_file(filepath):
    """Validate that a file exists and contains valid JSON.
    @param filepath: Path to the JSON file
    @returns: Parsed JSON data
    @raises argparse.ArgumentTypeError: If validation fails
    """
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError(f"File '{filepath}' does not exist")

    if not os.path.isfile(filepath):
        raise argparse.ArgumentTypeError(f"Path '{filepath}' is not a file")

    try:
        with open(filepath, "r") as f:
            _j = json.load(f)
            _j["__FILEPATH"] = os.path.abspath(filepath)
            return _j
    except json.JSONDecodeError as e:
        raise argparse.ArgumentTypeError(f"File '{filepath}' is not valid JSON: {e}")
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Cannot read file '{filepath}': {e}")


def validate_upload_path(dirpath):
    """Validate or create the upload directory.
    @param dirpath: Path to the upload directory
    @returns: The validated directory path
    @raises argparse.ArgumentTypeError: If validation fails
    """
    if os.path.exists(dirpath):
        if not os.path.isdir(dirpath):
            raise argparse.ArgumentTypeError(
                f"Path '{dirpath}' exists but is not a directory"
            )
    else:
        try:
            os.makedirs(dirpath, exist_ok=True)
            print(f"Created upload directory: {dirpath}")
        except Exception as e:
            raise argparse.ArgumentTypeError(
                f"Cannot create directory '{dirpath}': {e}"
            )

    return dirpath


def get_file_paths(config):
    """Extract file paths from the application configuration.
    @param config: Application configuration dictionary
    @returns: List of file paths
    """
    file_paths = []
    if "files" not in config:
        raise ValueError(
            "Application configuration must contain a 'files' key with a list of file paths relevent to root directory of configuration file."
        )
    input_files = config.get("files", [])
    _exists = []
    for f in input_files:
        abs_path = os.path.abspath(
            os.path.join(os.path.dirname(config.get("__FILEPATH", None)), f)
        )
        file_paths.append(abs_path)
        _exists.append(os.path.exists(abs_path))

    if not all(_exists):
        missing_files = [
            file_paths[i] for i, exists in enumerate(_exists) if not exists
        ]
        raise FileNotFoundError(
            f"The following input files specified in the configuration do not exist: {missing_files}"
        )

    return file_paths


def carry_pipeline_components(config, target_dir):
    """
    Copy pipeline component directories and files to the target directory.

    @param config: Application configuration dictionary
    @param target_dir: Target directory to copy components into
    """
    if "copy_files" in config:
        for f in config["copy_files"]:
            abs_path = os.path.abspath(
                os.path.join(os.path.dirname(config.get("__FILEPATH", None)), f)
            )
            if os.path.isdir(abs_path):
                shutil.copytree(abs_path, os.path.join(target_dir, os.path.basename(f)))
            else:
                shutil.copy(abs_path, os.path.join(target_dir, os.path.basename(f)))
    return


def link_files(files, target_dir):
    """
    Create symbolic links for the specified files in the target directory.

    @param files: List of file paths to link
    @param target_dir: Target directory to create links in
    """
    for f in files:
        abs_path = os.path.abspath(f)
        link_name = os.path.join(target_dir, os.path.basename(rename(f)))
        if not os.path.exists(link_name):
            os.symlink(abs_path, link_name)
    return


def touch_files(files, target_dir):
    """Create an empty file at the specified path.
    @param filepath: Path to the file to be created
    """
    if isinstance(files, str):
        files = [files]

    files = [os.path.join(target_dir, f) for f in files]
    for filepath in files:
        filepath = os.path.abspath(filepath)
        if not os.path.exists(os.path.dirname(filepath)):
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, "a"):
            os.utime(filepath, None)
    return


def main(args):
    SNK_CFG_ROOT = os.path.abspath(
        os.path.join(os.path.dirname(args.app_config.get("__FILEPATH", None)))
    )
    INPUT_FILES = get_file_paths(args.app_config)
    INJECT_TEMPLATE = json.dumps(args.app_config.pop("inject", {}))

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Initialize file paths
        NEW_SNK_FILE = os.path.join(tmpdirname, "Snakefile")
        CFG_FILE = os.path.join(tmpdirname, "config.json")
        BIN_DIR = os.path.join(tmpdirname, "bin")
        DATA_DIR = os.path.dirname(INPUT_FILES[0])
        DOT_GRAPH_FILE = os.path.join(tmpdirname, "pipeline_rulegraph.dot")
        if os.path.isdir(args.upload_path):
            RULE_GRAPH_FILENAME = os.path.splitext(os.path.basename(args.snakemake_config['__FILEPATH']))[0]
            FINAL_RULE_GRAPH = os.path.abspath(os.path.join(args.upload_path, f"{RULE_GRAPH_FILENAME}.svg"))
        else:
            FINAL_RULE_GRAPH = os.path.abspath(args.upload_path)

        # Copy snakefile
        shutil.copy(
            os.path.abspath(
                os.path.join(SNK_CFG_ROOT, args.app_config.get("snakefile", None))
            ),
            NEW_SNK_FILE,
        )

        # Copy pipeline subfolder components
        carry_pipeline_components(args.app_config, tmpdirname)

        # Link data files into pipeline directory
        link_files(INPUT_FILES, tmpdirname)

        # Touch stump files if specified in config
        if "stump" in args.app_config:
            touch_files(args.app_config["stump"], tmpdirname)

        # Template pipeline config with snakevis context
        if INJECT_TEMPLATE:
            template = jinja2.Template(INJECT_TEMPLATE)
            rendered_inject = json.loads(
                template.render(
                    {
                        "output_dir": tmpdirname,
                        "input_files": INPUT_FILES,
                        "data_dir": DATA_DIR,
                        "bin_dir": BIN_DIR,
                    }
                )
            )
            rendered_inject["options"]["input"] = ast.literal_eval(
                rendered_inject["options"]["input"]
            )

        # NOTE: dict.update does not do a deep merge, so we need to update nested dicts manually
        for k, v in rendered_inject.items():
            if k in args.snakemake_config and isinstance(args.snakemake_config[k], dict):
                args.snakemake_config[k].update(v)
            else:
                args.snakemake_config[k] = v

        # Write Snakemake config to pipeline directory
        with open(CFG_FILE, "w") as f:
            json.dump(args.snakemake_config, f, indent=2)
        
        # Call snakemake dry run, make dot rulegraph
        cmd = ''.join(
            "snakemake "
            f"--configfile {CFG_FILE} "
            f"-s {NEW_SNK_FILE} "
            f"-d {tmpdirname} "
            "--rulegraph "
            "--forceall "
            f"> {DOT_GRAPH_FILE}"
        )
        proc = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(
                f"Error: Snakemake dry run failed with exit code {proc.returncode}\n{stderr.decode()}"
            )
        
        # Call snakevis to convert dot to rule graph svg
        cmd = f"snakevision -s all multiqc -o {FINAL_RULE_GRAPH} {DOT_GRAPH_FILE}"
        proc = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(
                f"Error: Snakevis failed with exit code {proc.returncode}\n{stderr.decode()}"
            )
        if os.path.exists(FINAL_RULE_GRAPH):
            print(f"Snakemake pipeline visualization generated at: {FINAL_RULE_GRAPH}")
        else:
            raise FileNotFoundError(
                f"Failed to generate Snakemake pipeline visualization at: {FINAL_RULE_GRAPH}"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Snakemake pipeline visualization using JSON configuration"
    )
    parser.add_argument(
        "-c",
        "--app-config",
        type=validate_json_file,
        required=True,
        metavar="JSON_FILE",
        help="Path to application configuration JSON file",
    )
    parser.add_argument(
        "-s",
        "--snakemake-config",
        type=validate_json_file,
        required=True,
        metavar="JSON_FILE",
        help="Path to Snakemake configuration JSON file",
    )
    parser.add_argument(
        "-u",
        "--upload-path",
        type=validate_upload_path,
        required=True,
        metavar="DIRECTORY",
        help="Path where generated visualization files will be uploaded (will be created if it doesn't exist)",
    )
    args = parser.parse_args()
    main(args)
