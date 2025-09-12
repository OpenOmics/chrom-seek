# Chrom-seek Manual DiffBind Analyses and Downstream Steps

This documentation provides step-by-step instructions for performing post-pipeline analysis 
workflows outside of Snakemake, using the same working directory and scripts available 
in the `/bin` directory of chrom-seek.

## Overview

After running the chrom-seek pipeline, you may want to perform additional analysis steps 
manually, especially if you don't want to retain the raw fastqs and sorted bams. This 
guide covers four main stages:

1. Creating new DiffBind CSV files based on DiffBindQC CSV file
2. Running DiffBind in QC or differential peak calling mode
3. Annotating peaks with UROPA
4. Merging annotations with differential peak calls (for differential analysis)

## Prerequisites

- Completed chrom-seek pipeline run
- Access to the pipeline working directory
- Access to the docker image containing DiffBind v2
- Python environment (for merge)

## Stage 1: Creating New DiffBind CSV Files

### Purpose
DiffBind requires a CSV file with important metadata (in a specific
format) including information on the location of important sample
files (Q5DD bams and peak files). Manually consolidating this information
can be time-consuming and confusing, but the chrom-seek DiffBindQC
rules produce CSVs with information on all samples in the directory
in one place.

Steps:  
1. Grab the relevant CSV from within `<PIPELINE_OUTPUT_DIR>/PeakQC/DB_QC`
2. Remove rows with samples that should not be included in the analysis.  
3. Update information in ONLY the following three columns:  
   1. **CONDITION:**  
      This column contains the group information. This is the column
      DiffBind will use to define the contrast when running the differential binding
      analysis scripts. Do not start the values with numbers.  
   2. **TREATMENT:**
      Feel free to add this column as needed. The differential binding analysis
      scripts use this column for blocking/paired analysis when available. Do not start 
      the values with numbers.  
   3. **REPLICATE:**  
      This is a simple counting device for the tool and nothing more.
4. For differential analysis using scripts in `/bin`, you will want to sort your rows 
   so that the top rows are group1 for the contrast group1-group2.

## Stage 2: Running DiffBind Analysis

### Purpose
The main goal of this stage is to use the CSV created in Stage 1 and DiffBind to
identify consensus peaks and normalize the read counts following DiffBind v2 run
parameters. "QC mode analysis" will produce a TMM-normalized counts table for use with an external
differential peak caller of choice or to produce heatmaps. It will also calculate
PCA and UMAP coordinates for figures in a paper. "Differential peak calling analysis"
will walk you through the steps to use DiffBind for identifying differential peaks 
as accomplished in the chrom-seek pipeline.

### Input Files Required
- DiffBind CSV file (from Stage 1)
- BAM (and associated BAI) files and peak files referenced in the CSV

### Prepare Run Conditions

1. Get an interactive session with 100G mem, 200G lscratch, and 4 threads. These values are
an estimate based upon a typical chrom-seek project and may need to be adjusted for your
particular samples. 

  ```bash
  sinteractive -N 1 -n 1 --time=1-12:00:00 --mem=100G --gres=lscratch:200 --cpus-per-task=4
  ```
  
2. Initialize the singularity container for Diffbind v2. 

  > [!NOTE]
  > Exact conditions chosen here are required for lscratch to remain accessible and allow
    for the CSV file from Stage 1 to continue to have full paths.
   
  **Usage**:
  ```bash
    singularity run \
      -C \
      -e \
      --env TMPDIR=/tmp,TMP=/tmp \
      -B /lscratch/$SLURM_JOBID:/tmp,<PROJECT_WORKING_DIR>:<PROJECT_WORKING_DIR>:rw \
      docker://skchronicles/cfchip_toolkit:v0.5.0 \
      bash
  ```

  **Example**:
  ```bash
    working_dir=/data/OpenOmics/dev/datasets/outputs/test_homer
    singularity run \
      -C \
      -e \
      --env TMPDIR=/tmp,TMP=/tmp,working_dir=$working_dir \
      -B /lscratch/$SLURM_JOBID:/tmp,$working_dir:$working_dir:rw \
      docker://skchronicles/cfchip_toolkit:v0.5.0 \
      bash
  ```
  
### QC Mode Analysis

1. Run `DiffBind_v2_QC.Rmd` using example code.

  > [!NOTE]
  > Rmarkdown requires full paths or else it looks for the `<INPUT_CSV_FILE>` in `<WORKING_DIR>/bin`

  > [!NOTE]
  > The `<PEAK_TOOL>` token should be one of macsNarrow, Genrich, macsBroad, SEACR
  
  **Usage**:
  ```bash
  Rscript -e 'rmarkdown::render("<WORKING_DIR>/bin/DiffBind_v2_QC.Rmd", 
              output_file="<OUTPUT_HTML_FILE>",
              params=list(csvfile="<INPUT_CSV_FILE>", counts_bed="<OUTPUT_BED_FILE>", 
              counts_csv="<OUTPUT_COUNTS_CSV_FILE>", peakcaller="<PEAK_TOOL>"))'
  ```

  **Example**:
  ```bash
  Rscript -e 'rmarkdown::render("/data/OpenOmics/dev/datasets/outputs/test_homer/bin/DiffBind_v2_QC.Rmd", output_file="/data/OpenOmics/dev/datasets/outputs/test_homer/PeakQC/DB_QC/AllSamples-macsNarrow/AllSamples-macsNarrow_DiffBindQC.html",
    params=list(csvfile="/data/OpenOmics/dev/datasets/outputs/test_homer/PeakQC/DB_QC/AllSamples-macsNarrow/AllSamples-macsNarrow_DiffBind_prep.csv", counts_bed="/data/OpenOmics/dev/datasets/outputs/test_homer/PeakQC/DB_QC/AllSamples-macsNarrow/AllSamples-macsNarrow_DiffBindQC_TMMcounts.bed", 
    counts_csv="/data/OpenOmics/dev/datasets/outputs/test_homer/PeakQC/DB_QC/AllSamples-macsNarrow/AllSamples-macsNarrow_DiffBindQC_TMMcounts.csv", peakcaller="macsNarrow"))'
  ```

### Differential Peak Calling Analysis ###

1. Run `DiffBind_v2_load.R`.

  > [!NOTE]
  > The `<PEAK_TOOL>` token should be one of macsNarrow, Genrich, macsBroad, SEACR
  
  **Usage**:
  ```bash
  <WORKING_DIR>/bin/DiffBind_v2_load.R \
    --csvfile <INPUT_CSV_FILE> \
    --counts <INPUT_RDS_FILE> \
    --list <OUTPUT_PEAK_BED_FILE> \
    --peakcaller <PEAK_TOOL>
  ```

  **Example**:
  ```bash
  /data/OpenOmics/dev/datasets/outputs/test_homer/bin/DiffBind_v2_load.R \
    --csvfile /data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsBroad/IFN0h_vs_IFN24h-macsBroad_Diffbind_prep.csv \
    --counts /data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsBroad/IFN0h_vs_IFN24h-macsBroad_Diffbind_counts.rds \
    --list /data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsBroad/IFN0h_vs_IFN24h-macsBroad_Diffbind_fullList.bed \
    --peakcaller macsBroad
  ```

2. Execute differential comparisons using `DiffBind v2.15.2`.

   1. Determine which differential application suits your needs: `DeSeq2 or edgeR`
     - large variance in library size
     - library size is confounding
     - prior expectations that there should be an equivalent number of peaks on both sides of the contrast

   2. Verify if your experimential design requires the uses of `blocking` or `no blocking`.

   3. Find relevant script (`<ANALYSIS_SCRIPT>`): 
     - <WORKING_DIR>/bin/DiffBind_v2_Deseq2.Rmd  
     - <WORKING_DIR>/bin/DiffBind_v2_Deseq2_block.Rmd  
     - <WORKING_DIR>/bin/DiffBind_v2_EdgeR.Rmd  
     - <WORKING_DIR>/bin/DiffBind_v2_EdgeR_block.Rmd

   4. Establish group contrast from experimental setup: "{group1}_vs_{group2}" [string] = `<INPUT_CONTRASTS>`

   5. Execute script:
      > The `<PEAK_TOOL>` token should be one of macsNarrow, Genrich, macsBroad, SEACR

      **Usage**:
      ```bash
      Rscript -e 'rmarkdown::render("<ANALYSIS_SCRIPT>", 
                  output_file="<OUTPUT_DIFFBIND_REPORT_FILE>",
                  params=list(csvfile="<INPUT_CSV_FILE>", 
                  peakcaller="<PEAK_TOOL>", 
                  list_file="<OUTPUT_PEAK_BED_LIST>",
                  up_file="<OUTPUT_UP_REGULATED_FILE>", 
                  down_file="<OUTPUT_DOWN_REGULATED_FILE>", 
                  contrasts="<INPUT_CONTRASTS>", 
                  counts="<INPUT_PEAK_COUNTS_CSV_FILE>"))'
      ```

      **Example**:
      ```bash
      Rscript -e 'rmarkdown::render("/data/OpenOmics/dev/datasets/outputs/test_homer/bin/DiffBind_v2_Deseq2.Rmd", 
                  output_file="/data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_DeSeq2.html", 
                  params=list(csvfile="/data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_prep.csv", 
                  peakcaller="macsNarrow", 
                  list_file="/data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_DeSeq2_peak_list.tab", 
                  up_file="/data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_DeSeq2_up.bed", 
                  down_file="/data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_DeSeq2_down.bed", 
                  contrasts="IFN0h_vs_IFN24h", counts="/data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_counts.rds"))'
      ```

## Stage 3: Peak Annotation with UROPA

### Purpose
Annotate genomic peaks with gene information, regulatory elements, and genomic context using UROPA.

### Input Files Required
- Peak files (BED or narrowPeak format) `<PEAK_FILE>`
- GTF annotation file `<GTF_FILE>`
- UROPA configuration file `<UROPA_CONFIG>`

The GTF file used within the pipeline can be found in the config.json under [references][<SPECIES>][GTFFILE].

### Steps

1. Create UROPA configuration file. The example below is what is currently used for chrom-seek
    when assay is not cfchip.

  ```json
  {
    "queries": [
        {
            "feature": ["gene"],
            "filter_attribute": "gene_type",
            "attribute_values": ["protein_coding"],
            "feature_anchor": ["start"],
            "distance": [3000, 1000],
            "name": "query_1"
        },
        {
            "feature": ["gene"],
            "filter_attribute": "gene_type",
            "attribute_values": ["protein_coding"],
            "feature_anchor": ["start"],
            "distance": [10000],
            "name": "query_2"
        },
        {
            "feature": ["gene"],
            "filter_attribute": "gene_type",
            "attribute_values": ["protein_coding"],
            "feature_anchor": ["start"],
            "distance": [100000],
            "name": "query_3"
        }
    ],
    "show_attributes": [
        "gene_id",
        "gene_name",
        "gene_type"
    ],
    "priority": "True",
    "gtf": "<GTF_FILE>",
    "bed": "<PEAK_FILE>"
  }
  ```

2. Run UROPA on Biowulf. This example uses the uropa module, but a path to a singularity 
object is also available within the config.json. We suggest 4 threads and 10G memory.

  ```bash
  module load uropa/4.0.2
  uropa -i sample_config.json -t 4
  ```

## Stage 4: Merging Annotations with Differential Results

### Purpose
Combine differential binding results from DiffBind with gene annotations from UROPA to 
create comprehensive results tables.

### Input Files Required
- Differential peaks results from Stage 2
- UROPA annotation results from Stage 3

### Steps

1. **Merge differential results with annotations:**
  1. Collect DiffBind TAB output from **Differential Peak Calling Analysis**: `<OUTPUT_PEAK_BED_LIST>` 
  2. Collect UROPA *_finalhits.txt output for input to merge script: `<UROPA_FINAL_HITS>`
  3. Collect merge file output file location: `<MERGE_OUTPUT>`
  2. Decide if you want to filter at a higher FDR or fold change than default UROPA settings [Optional]: `<FDR_THRESHOLD>`, `<FC_THRESHOLD>`

  **Usage**:
  ```bash
  python $working_dir/bin/merge_diffbind_uropa.py \
    --diffbind <DIFFBIND_PEAK_LIST_FILE> \
    --uropa <UROPA_FINAL_HITS> \
    --fdr <FDR_THRESHOLD> \
    --fold <FC_THRESHOLD> \
    --output <MERGE_OUTPUT>
  ```

  **Example**:
  ```bash
  python /data/OpenOmics/dev/datasets/outputs/test_homer/bin/merge_uropa_diffbind.py \
    --uropa /data/OpenOmics/dev/datasets/outputs/test_homer/UROPA_annotations/DiffBind/IFN0h_vs_IFN24h-macsNarrow-EdgeR/IFN0h_vs_IFN24h_macsNarrow_EdgeR_protTSS_uropa_finalhits.txt \
    --diffbind /data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_EdgeR_peak_list.tab \
    --output /data/OpenOmics/dev/datasets/outputs/test_homer/UROPA_DIFFBIND_TBLS/IFN0h_vs_IFN24h-macsNarrow-EdgeR_protTSS_UROPA_DIFFBIND_JOIN.txt \
    --fdr 0.05 \
    --fold 0
  ```

## Getting Help

- See [biowulf user guide](https://hpc.nih.gov/docs/userguide.html) for more help with `sinteractive`.
- See [Diffbind v2.0.2](https://www.rdocumentation.org/packages/DiffBind/versions/2.0.2) 
documentation to learn about the functions used.  You may also find the 
[Diffbind v1.2.4](https://www.bioconductor.org/packages//2.10/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf)
vignette informative.
- For more information on UROPA run parameters, see: https://uropa-manual.readthedocs.io/

