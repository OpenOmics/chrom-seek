#!/usr/bin/env python3
import argparse
import pandas as pd
import os
from glob import glob
from typing import Dict, Tuple


files_out = {
    'mapping_summary': 'MappingSummary.txt',
    'encode_qc': 'EncodeQC.txt',
    'peak_metrics_summary': 'PeakSummary.txt'
}


def parse_picard_insert_size_metrics(input_file: str) -> Tuple[Dict, pd.DataFrame]:
    """
    Parse Picard CollectInsertSizeMetrics output file.
    
    The file contains:
    1. Header comments starting with #
    2. Metrics section with column headers and one data row
    3. Empty line or continuation
    4. Histogram section with column headers and multiple data rows
    
    Args:
        input_file: Path to the Picard output file
        
    Returns:
        Tuple of (metrics_dict, histogram_dataframe)
    """
    # Read and parse the file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    # Remove header comments and empty lines
    data_lines = []
    for line in lines:
        if not line.startswith('#') and line.strip():
            data_lines.append(line.strip())
    if len(data_lines) < 3:
        raise ValueError("File does not contain expected metrics and histogram sections")
    # Parse metrics section (first two non-comment lines)
    metrics_headers = data_lines[0].split('\t')
    metrics_values = data_lines[1].split('\t')
    # Create metrics dictionary
    metrics = {}
    for header, value in zip(metrics_headers, metrics_values):
        # Try to convert to appropriate type
        try:
            if '.' in value and value.replace('.', '').replace('-', '').isdigit():
                metrics[header] = float(value)
            elif value.replace('-', '').isdigit():
                metrics[header] = int(value)
            else:
                metrics[header] = value
        except (ValueError, AttributeError):
            metrics[header] = value
    # Find histogram section (starts after empty line or metrics)
    histogram_start = 2
    # Look for histogram header (typically starts with "insert_size")
    while histogram_start < len(data_lines):
        if 'insert_size' in data_lines[histogram_start].lower():
            break
        histogram_start += 1
    if histogram_start >= len(data_lines):
        raise ValueError("Could not find histogram section")
    # Parse histogram data
    histogram_headers = data_lines[histogram_start].split('\t')
    histogram_rows = []
    for i in range(histogram_start + 1, len(data_lines)):
        row_values = data_lines[i].split('\t')
        if len(row_values) == len(histogram_headers):
            histogram_rows.append(row_values)
    # Create histogram DataFrame
    histogram_df = pd.DataFrame(histogram_rows, columns=histogram_headers)
    # Convert numeric columns in histogram
    for col in histogram_df.columns:
        try:
            histogram_df[col] = pd.to_numeric(histogram_df[col])
        except (ValueError, TypeError):
            pass  # Keep as string if conversion fails
    return metrics, histogram_df


def parse_flagstat_text(flagstat_output):
    """
    Parses the standard text output from samtools flagstat.

    Args:
        flagstat_output (str): The raw text output from `samtools flagstat`.

    Returns:
        dict: A dictionary containing the parsed flagstat statistics.
    """
    parsed_stats = {}
    for line in flagstat_output:
        if '\n' in line:
            line = line.strip()
        parts = line.split(' + ')
        if len(parts) == 2:
            qc_pass = int(parts[0].strip().split(' ')[0])
            qc_fail_and_desc = parts[1].strip().split(' ', 1)
            qc_fail = int(qc_fail_and_desc[0])
            description = qc_fail_and_desc[1].strip()
            if description.startswith('mapped ('):
                description = 'mapped'
            parsed_stats[description] = {"QC_Pass": qc_pass, "QC_Fail": qc_fail}
        else:
            # Handle the total reads line (e.g., "1000000 in total")
            if "in total" in line:
                total_reads = int(line.split(' ')[0])
                parsed_stats["Total Reads"] = total_reads
    return parsed_stats


def percent(numerator, denominator):
    if denominator == 0:
        return 0
    return round((numerator / denominator) * 100, 2)


def main(args):
    workpath = args.workpath
    # mapping stats
    flagdata = {}
    flagstats = glob(f"{workpath}/**/*flagstat", recursive=True)
    for flagfile in flagstats:
        name = os.path.basename(flagfile)
        this_data = parse_flagstat_text(open(flagfile).readlines())
        flagdata[name] = this_data
    tagged_flag_data = {}
    for flagkey, this_flagdata in flagdata.items():
        sample_name = flagkey.split(".")[0]
        file_context = flagkey.split(".")[1]
        if sample_name not in tagged_flag_data:
            tagged_flag_data[sample_name] = {}
        if file_context not in tagged_flag_data[sample_name]:
            tagged_flag_data[sample_name][file_context] = this_flagdata
    mapping_summary = pd.DataFrame()
    for sample, contexts in tagged_flag_data.items():
        assert 'Q5' in contexts, f"Sample {sample} missing Q5 flagstat data"
        assert 'Q5DD' in contexts, f"Sample {sample} missing Q5DD flagstat data"
        assert 'sorted' in contexts, f"Sample {sample} missing sorted flagstat data"
        map_summary_row = {
            'SampleID': sample,
            'TrimmedReads': contexts['sorted']['read1']['QC_Pass'],
            'AlignedReads': contexts['sorted']['mapped']['QC_Pass'],
            'QualityReads': contexts['Q5']['read1']['QC_Pass'],
            'DedupReads': contexts['Q5DD']['read1']['QC_Pass'],
            'PercentDuplicated': percent(contexts['Q5']['read1']['QC_Pass'], contexts['Q5DD']['read1']['QC_Pass'])
        }
        mapping_summary = pd.concat([mapping_summary, pd.DataFrame([map_summary_row])], ignore_index=True)
    mapping_summary.reset_index(drop=True, inplace=True)
    mapping_summary.to_csv(os.path.join(workpath, files_out['mapping_summary']), index=False)

    # encode stats
    if args.paired:
        insert_size_data = {}
        q5dd_insert_size = glob(f"{workpath}/**/*Q5DD.insert_size_metrics.txt")
        nrf_data = {}
        nrf = glob(f"{workpath}/**/*nrf")
        for insertfile in q5dd_insert_size:
            name = os.path.basename(insertfile)
            metrics, histogram = parse_picard_insert_size_metrics(insertfile)
            insert_size_data[name] = {
                "metrics": metrics,
                "histogram": histogram
            }
        for nrf_file in nrf:   
            name = os.path.basename(nrf_file).split('.')[0]
            with open(nrf_file, 'r') as f:
                nrf_value = f.readlines()
            nrf_value = [line.strip().split('\t') for line in nrf_value]
            assert len(nrf_value) == 1, f"NRF file {nrf_file} has unexpected format"
            nrf_value = nrf_value[0]
            nrf_data[name] = nrf_value

        encode_qc = pd.DataFrame()
        for sample, data in insert_size_data.items():
            sample_name = sample.split(".")[0]
            metrics, histogram = data['metrics'], data['histogram']
            this_nrf_data = nrf_data.get(sample_name, ['NA', 'NA', 'NA'])
            encode_row = {
                'SampleID': sample_name,
                'MEAN_INSERT_SIZE': metrics.get('MEAN_INSERT_SIZE', 'NA'),
                'NRF': this_nrf_data[0],
                'PBC1': this_nrf_data[1],
                'PBC2': this_nrf_data[2]
            }
            encode_qc = pd.concat([encode_qc, pd.DataFrame([encode_row])], ignore_index=True)
        encode_qc.reset_index(drop=True, inplace=True)
        encode_qc.to_csv(os.path.join(workpath, files_out['encode_qc']), index=False)
    # peak metrics  
    q5dd_peak_metrics = glob(f"{workpath}/PeakQC/FRiP/**/*.FRiP_table.txt")
    peak_data = pd.DataFrame()
    for peakfile in q5dd_peak_metrics:
        df = pd.read_csv(peakfile, sep="\t")
        peak_data = pd.concat([peak_data, df], ignore_index=True)
    peak_callers = {
        "macsNarrow": "*_peaks.narrowPeak",
        "macsBroad": "*_peaks.broadPeak",
        "SEACR": "*.stringent.bed",
        "genrich": "*.narrowPeak",
    }
    peaks_called_data = {}
    for caller in peak_callers:
        if os.path.exists(f"{workpath}/{caller}"):
            get_peak_files = glob(f"{workpath}/{caller}/**/{peak_callers[caller]}", recursive=True)
            for peakfile in get_peak_files:
                sample_name = os.path.basename(peakfile).split(peak_callers[caller].replace('*', ''))[0]
                caller = os.path.normpath(peakfile).split(os.sep)[-3]
                num_lines = sum(1 for _ in open(peakfile)) - 1  # subtract header
                if sample_name not in peaks_called_data:
                    peaks_called_data[sample_name] = {}
                peaks_called_data[sample_name][caller] = num_lines
    peak_metric_summary = pd.DataFrame()
    aligned_peak_data = peak_data[peak_data['bedsample'] == peak_data['bamsample']]
    for sample in peak_data['bedsample'].unique():
        pk_call_row = {'SampleID': sample}
        this_peak_data = aligned_peak_data[aligned_peak_data['bedsample'] == sample]
        pk_caller_count = 1
        for _, row in this_peak_data.iterrows():
            caller = row['bedtool']
            peak_counts = peaks_called_data[sample][caller]
            pk_call_row['PeakCaller' + str(pk_caller_count)] = row['bedtool']
            pk_call_row['PeakCaller' + str(pk_caller_count) + ' Peaks Called'] = str(peak_counts)
            pk_call_row['PeakCaller' + str(pk_caller_count) + ' FRiP'] = row['FRiP']
            pk_caller_count += 1
        peak_metric_summary = pd.concat([peak_metric_summary, pd.DataFrame([pk_call_row])], ignore_index=True)
    peak_metric_summary.reset_index(drop=True, inplace=True)
    peak_metric_summary.to_csv(os.path.join(workpath, files_out['peak_metrics_summary']), index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='A script that takes the working directory and collects summary statistics into three files.'
    )
    parser.add_argument(
        'workpath',
        help='The root directory for all output of the pipeline.'
    )
    parser.add_argument(
        '--paired',
        action='store_true',
        default=False, 
        help='True/False flag for endedness.'
    )
    main(parser.parse_args())