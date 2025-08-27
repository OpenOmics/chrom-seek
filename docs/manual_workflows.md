# Chrom-seek Manual Workflows Documentation

This documentation provides step-by-step instructions for performing post-pipeline analysis workflows outside of Snakemake, using the same working directory and scripts available in the `/bin` directory of chrom-seek.

## Overview

After running the chrom-seek pipeline, you may want to perform additional analysis steps manually, especially if you don't want to retain the raw fastqs. This guide covers four main steps:

1. Creating new DiffBind CSV files based on DiffBindQC CSV file
2. Running DiffBind in QC or differential peak calling mode
3. Annotating peaks with UROPA
4. Merging annotations with differential peak calls (for differential analysis)

## Prerequisites

- Completed chrom-seek pipeline run
- Access to the pipeline working directory
- Access to the docker image containing DiffBind v2
- Python environment (for merge)

## Workflow 1: Creating New DiffBind CSV Files

### Purpose
DiffBind requires a CSV file with important metadata (in a specific
format) including information on the location of important sample
files (bams and peak files). Manually consolidating this information
can be time-consuming and confusing, but the chrom-seek DiffBindQC
steps produce CSVs with information on all samples in the directory
in one place.

1. Grab the relevant CSV from within `<PIPELINE_OUTPUT_DIR>/PeakQC/DB_QC`
2. Remove rows with samples that should not be included in the analysis.  
3. Update information in ONLY the following three columns:  
   1. **CONDITION:**  
      This column contains the group information. This is the column
      DiffBind will use to define the contrast when running the differential binding
      analysis scripts. Do not start the values with numbers.  
   2. **TREATMENT:**
      Feel free to add this column as needed. The differential binding analysis
      scripts use this column for blocking/paired analysis when available. Do not start the values with numbers.  
   3. **REPLICATE:**  
      This is a simple counting device for the tool and nothing more.
4. For differential analysis, you will want to sort your rows so that the top rows are group1 for the contrast group1-group2.

## Workflow 2: Running DiffBind Analysis

### Purpose
Perform either quality control analysis or differential binding analysis using DiffBind.

### Input Files Required
- DiffBind CSV file (from Workflow 1)
- BAM files and peak files referenced in the CSV

### QC Mode Analysis

1. Get an interactive session with 100G mem, 200G lscratch, and 4 threads. See [biowulf user guide](https://hpc.nih.gov/docs/userguide.html) for more help with `sinteractive`.
  ```bash
  sinteractive -N 1 -n 1 --time=1-12:00:00 --mem=100G --gres=lscratch:200 --cpus-per-task=4
  ```
2. Initialize the singularity container for Diffbind v2. See [singularity run documentation](https://docs.sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) and singularity [metadata and environment guide](https://docs.sylabs.io/guides/3.7/user-guide/environment_and_metadata.html) for more help with `singularity run`.  
  **Usage**:
  ```bash
    singularity run \
      -C \
      -e \
      --env TMPDIR=/tmp,TMP=/tmp \
      -B /lscratch/$SLURM_JOBID:/tmp,<PROJECT_WORKING_DIR>:/work:rw \
      --pwd /work \
      docker://skchronicles/cfchip_toolkit:v0.5.0 \
      bash
  ```

  **Example**:
  ```bash
    working_dir=/data/OpenOmics/project1
    singularity run \
      -C \
      -e \
      --env TMPDIR=/tmp,TMP=/tmp \
      -B /lscratch/$SLURM_JOBID:/tmp,$working_dir:/work:rw \
      --pwd /work \
      docker://skchronicles/cfchip_toolkit:v0.5.0 \
      bash
  ```
3. Run `DiffBind_v2_QC.Rmd` using example code. See [Diffbind v3.12.0](https://bioconductor.statistik.tu-dortmund.de/packages/3.18/bioc/html/DiffBind.html) documentation for reference.

  > [!NOTE]
  > Contextual output tokens (`<OUTPUT_*>`) need to point to a writable location (/work)

  > [!NOTE]
  > The `<PEAK_TOOL>` token should be one of macsNarrow, Genrich, macsBroad, SEACR
  
  **Usage**:
  ```bash
  Rscript -e 'rmarkdown::render("/work/bin/DiffBind_v2_QC.Rmd", output_file="/work/CUSTOM_DiffBindQC.html",
              params=list(csvfile="<INPUT_CSV_FILE>", counts_bed="<OUTPUT_BED_FILE>", 
              counts_csv="<OUTPUT_COUNTS_CSV_FILE>", peakcaller="<PEAK_TOOL>"))'
  ```

  **Example**:
  ```bash
  Rscript -e 'rmarkdown::render("/work/bin/DiffBind_v2_QC.Rmd", output_file="/work/CUSTOM_DiffBindQC.html",
              params=list(csvfile="/data/OpenOmics/project1/diffbind_macsNarrow_sample1.csv", counts_bed="/output/diffbind_macsnarrow_counts.bed", 
              counts_csv="/output/diffbind_macsnarrow_counts.csv", peakcaller="macsNarrow"))'
  ```

### Differential Peak Calling Mode ###

1. Get an interactive session with 100G mem, 200G lscratch, and 4 threads.
  ```bash
  sinteractive -N 1 -n 1 --time=1-12:00:00 --mem=100G --gres=lscratch:200 --cpus-per-task=4
  ```
2. Initialize the singularity container for Diffbind v2. See [singularity run documentation](https://docs.sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) and singularity [metadata and environment](https://docs.sylabs.io/guides/3.7/user-guide/environment_and_metadata.html) for more help with `singularity run`. 

  **Usage**:
  ```bash
  singularity run \
    -C \
    -e \
    --env TMPDIR=/tmp,TMP=/tmp \
    -B /lscratch/$SLURM_JOBID:/tmp,<PROJECT_WORKING_DIRECTORY>:/work:rw \
    --pwd /work \
    docker://skchronicles/cfchip_toolkit:v0.5.0 \
    bash
  ``` 

  **Example**:
  ```bash
  working_dir=/data/OpenOmics/project1
  singularity run \
    -C \
    -e \
    --env TMPDIR=/tmp,TMP=/tmp \
    -B /lscratch/$SLURM_JOBID:/tmp,$working_dir:/work:rw \
    --pwd /work \
    docker://skchronicles/cfchip_toolkit:v0.5.0 \
    bash
  ``` 

3. Run `DiffBind_v2_load.R`.
  > [!NOTE]
  > Contextual output tokens (`<OUTPUT_*>`) need to point to a writable location (`/work`)

  > [!NOTE]
  > The `<PEAK_TOOL>` token should be one of macsNarrow, Genrich, macsBroad, SEACR
  

  **Usage**:
  ```bash
  /work/bin/DiffBind_v2_load.R \
    --csvfile <INPUT_CSV_FILE> \
    --counts <INPUT_RDS_FILE> \
    --list <OUTPUT_PEAK_BED_FILE> \
    --peakcaller <PEAK_TOOL>
  ```

  **Example**:
  ```bash
  /work/bin/DiffBind_v2_load.R \
    --csvfile /data/OpenOmics/project1/diffbind_macsNarrow_sample1.csv \
    --counts /data/OpenOmics/project1/diffbind_macsNarrow_sample1.rds \
    --list /ouput/project1/diffbind_macsNarrow_sample1.bed \
    --peakcaller macsNarrow
  ```

4. Execute differential comparisons using `DiffBind v2.15.2`.

   1. Determine which differential application suits your needs: `DeSeq2 or edgeR`
     - large variance in library size
     - library size is confounding
     - prior expectations that there should be an equivalent number of peaks on both sides of the contrast

   2. Verify if your experimential design requires the uses of `blocking` or `no blocking`.

   3. Find relevant script (`<ANALYSIS_SCRIPT>`): 
     - /work/bin/DiffBind_v2_Deseq2.Rmd  
     - /work/bin/DiffBind_v2_Deseq2_block.Rmd  
     - /work/bin/DiffBind_v2_EdgeR.Rmd  
     - /work/bin/DiffBind_v2_EdgeR_block.Rmd

   4. Establish group contrast from experimental setup: "{group1}_vs_{group2}" [string] = `<INPUT_CONTRASTS>`

   5. Establish output locations for tokens: `<OUTPUT_DIFFBIND_REPORT_FILE>`, `<OUTPUT_UP_REGULATED_FILE>`, `<OUTPUT_DOWN_REGULATED_FILE>`, `<OUTPUT_PEAK_BED_LIST>`
      > Contextual output tokens (<OUTPUT_*>) need to point to a writable location (`/work`)

   6. Execute script:
      > The `<PEAK_TOOL>` token should be one of macsNarrow, Genrich, macsBroad, SEACR

      **Usage**:
      ```bash
      Rscript -e 'rmarkdown::render("<ANALYSIS_SCRIPT>", output_file="<OUTPUT_DIFFBIND_REPORT_FILE>",
        params=list(csvfile="<INPUT_CSV_FILE>", peakcaller="<PEAK_TOOL>", list_file="<OUTPUT_PEAK_BED_LIST>",
        up_file="<OUTPUT_UP_REGULATED_FILE>", down_file="<OUTPUT_DOWN_REGULATED_FILE>", contrasts="<INPUT_CONTRASTS>", counts="<INPUT_PEAK_COUNTS_CSV_FILE>"))'
      ```

      **Example**:
      ```bash
      Rscript -e 'rmarkdown::render("/work/bin/DiffBind_v2_Deseq2.Rmd", output_file="/output/project1/diffbind_deseq2_report.html",
        params=list(csvfile="/data/OpenOmics/project1/diffbind_macsNarrow_sample1.csv", peakcaller="macsNarrow", list_file="/ouput/project1/diffbind_macsNarrow_sample1.bed",
        up_file="/output/project1/up_regulated_genes.bed", down_file="/output/project1/down_regulated_genes.bed", contrasts="group1_vs_group2", counts="/data/OpenOmics/project1/diffbind_macsNarrow_sample1.rds"))'
      ```

## Workflow 3: Peak Annotation with UROPA

### Purpose
Annotate genomic peaks with gene information, regulatory elements, and genomic context using UROPA.

### Input Files Required
- Peak files (BED or narrowPeak format) `<PEAK_FILE>`
- GTF annotation file `<GTF_FILE>`
- UROPA configuration file `<UROPA_CONFIG>`

### Steps

1. Create UROPA configuration file:
  
  > [!NOTE]
  > Contextual output tokens (`<OUTPUT_*>`) need to point to a writable location (/work)

  ```json
  {
    "queries": [
        {
            "feature": ["gene"],
            "filter_attribute": "gene_type",
            "attribute_values": ["protein_coding"],
            "feature_anchor": ["start"],
            "relative_location": [
                "PeakInsideFeature",
                "FeatureInsidePeak",
                "Upstream",
                "Downstream",
                "OverlapStart",
                "OverlapEnd"
            ],
            "strand": "ignore",
            "distance": [3000, 1000],
            "name": "query_1"
        },
        {
            "feature": ["gene"],
            "filter_attribute": "gene_type",
            "attribute_values": ["protein_coding"],
            "feature_anchor": ["start"],
            "relative_location": [
                "PeakInsideFeature",
                "FeatureInsidePeak",
                "Upstream",
                "Downstream",
                "OverlapStart",
                "OverlapEnd"
            ],
            "strand": "ignore",
            "distance": [10000],
            "name": "query_2"
        },
        {
            "feature": ["gene"],
            "filter_attribute": "gene_type",
            "attribute_values": ["protein_coding"],
            "feature_anchor": ["start"],
            "relative_location": [
                "PeakInsideFeature",
                "FeatureInsidePeak",
                "Upstream",
                "Downstream",
                "OverlapStart",
                "OverlapEnd"
            ],
            "strand": "ignore",
            "distance": [100000],
            "name": "query_3"
        }
    ],
    "show_attributes": [
        "gene_id",
        "gene_name",
        "gene_type"
    ],
    "priority": "False",
    "gtf": "<GTF_FILE>",
    "bed": "<PEAK_FILE>"
  }
  ```

2. **Run UROPA annotation:**
  ```bash
  uropa -i sample_config.json -t 4
  ```

3. **Expected outputs:**
   - `<OUTPUT_UROPA_DIR>/Sample1_finalhits.txt` - Detailed annotations
   - `<OUTPUT_UROPA_DIR>/Sample1_allhits.txt` - Annotation summary statistics


### Configuration Options
The UROPA config file supports various annotation priorities:

- **Promoter regions**: TSS ± 3kb
- **Gene body**: Exons and introns
- **Intergenic regions**: Regions between genes
- **Custom features**: User-defined genomic regions

See [UROPA manual](https://uropa-manual.readthedocs.io/) for more information.

## Workflow 4: Merging Annotations with Differential Results

### Purpose
Combine differential binding results from DiffBind with gene annotations from UROPA to create comprehensive results tables.

### Input Files Required
- Differential peaks results from Workflow 2
- UROPA annotation results from Workflow 3

### Steps

1. **Merge differential results with annotations:**
  1. Collect diffbind CSV output from **Differential Peak Calling Mode**: `<OUTPUT_PEAK_COUNTS_CSV_FILE>` = `<DIFFBIND_CSV_FILE>`
  2. Collect UROPA *_finalhits.txt output for input to merge script: `<UROPA_FINAL_HITS>`
  3. Collect merge file output file location: `<MERGE_OUTPUT>`
  2. Decide if you want to filter at a higher FDR or fold change than default UROPA settings [Optional]: `<FDR_THRESHOLD>`, `<FC_THRESHOLD>`

  > [!NOTE]
  > Contextual output tokens (`<OUTPUT_*>`) need to point to a writable location (/work)

  **Usage**:
  ```bash
  python /work/bin/merge_diffbind_uropa.py \
    --diffbind <DIFFBIND_CSV_FILE> \
    --uropa <UROPA_FINAL_HITS> \
    --fdr <FDR_THRESHOLD> \
    --fold <FC_THRESHOLD> \
    --output <MERGE_OUTPUT>
  ```

  **Example**:
  ```bash
  python /work/bin/merge_diffbind_uropa.py \
    --diffbind /data/OpenOmics/project1/diffbind_macsNarrow_sample1.csv \
    --uropa /data/OpenOmics/project1/uropa_config.json \
    --fdr 0.05 \
    --fold 2 \
    --output /outputs/uropa_annotation_macsnarrow_diffbind_results.txt
  ```

2. **Expected outputs:**
   - `<MERGE_OUTPUT>` - Combined differential and annotation data

## Troubleshooting

### Common Issues and Solutions

1. **Missing BAM files**
   - Check file paths in the DiffBind CSV
   - Ensure BAM files are indexed (`.bai` files present)

2. **UROPA annotation errors**
   - Verify GTF file format and chromosome naming consistency
   - Check peak file format (BED vs narrowPeak)

3. **Memory issues with large datasets**
   - Reduce the number of cores used
   - Process samples in smaller batches
   - Consider using cluster resources

4. **R package dependencies**
   ```bash
   # Install required R packages
   Rscript -e "
   if (!requireNamespace('BiocManager', quietly = TRUE))
     install.packages('BiocManager')
   BiocManager::install(c('DiffBind', 'ChIPseeker', 'clusterProfiler'))
   "
   ```

### File Format Requirements

- **Peak files**: BED format or narrowPeak (tab-separated)
- **BAM files**: Properly sorted and indexed
- **CSV files**: Comma-separated with proper headers
- **GTF files**: Standard GTF format with gene annotations

## Best Practices

1. **File organization**: Keep input files organized in clear directory structures
2. **Backup results**: Make copies of important intermediate files
3. **Document parameters**: Keep track of analysis parameters used
4. **Quality control**: Always review QC plots before proceeding to differential analysis
5. **Reproducibility**: Use version-controlled scripts and document software versions

## Getting Help

For additional support:
- Check the main chrom-seek documentation: https://openomics.github.io/chrom-seek/
- Review the FAQ section for common issues
- Open an issue on GitHub: https://github.com/OpenOmics/chrom-seek/issues

This documentation assumes familiarity with command-line operations and basic bioinformatics concepts. For beginners, we recommend starting with the QC workflow before attempting differential analysis.
