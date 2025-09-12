# Chrom-seek Manual: DiffBind Analyses and Downstream Steps

This guide provides step-by-step instructions for performing post-pipeline analysis 
workflows outside of Snakemake, using the same working directory and scripts available 
in the `/bin` directory of chrom-seek.

## Overview

After running the chrom-seek pipeline, you may want to perform additional analysis steps 
manually, especially if you do not wish to retain the raw FASTQs and sorted BAMs. This 
guide covers four main stages:

1. [Creating new DiffBind CSV files](#stage-1-creating-new-diffbind-csv-files)  
2. [Running DiffBind in QC or differential peak calling mode](#stage-2-running-diffbind-analysis)  
3. [Annotating peaks with UROPA](#stage-3-peak-annotation-with-uropa)  
4. [Merging annotations with differential peak calls](#stage-4-merging-annotations-with-differential-results)  

## Prerequisites

- Completed chrom-seek pipeline run  
- Access to the pipeline working directory  
- Access to the Docker/Singularity image containing DiffBind v2  
- Python environment (for merging results)  

---

## Stage 1: Creating New DiffBind CSV Files

### Purpose
DiffBind requires a CSV file with metadata (in a specific format), including the locations 
of Q5DD BAMs and peak files. Manually consolidating this information can be time-consuming 
and error-prone, but the chrom-seek DiffBindQC rules generate CSVs with all sample 
information in one place.

### Steps
1. Retrieve the relevant CSV from `<PIPELINE_OUTPUT_DIR>/PeakQC/DB_QC`.  
2. Remove rows for samples that should not be included in the analysis.  
3. Update **only** the following three columns:  
   - **CONDITION:** Group information used by DiffBind to define contrasts.  
     - *Note:* Do not start values with numbers.  
   - **TREATMENT:** (Optional) Used for blocking/paired analysis when available.  
     - *Note:* Do not start values with numbers.  
   - **REPLICATE:** A simple count (1, 2, 3, …) per sample.  
4. For [differential analysis](#stage-2-running-diffbind-analysis), sort rows so that the first group (`group1`) appears before 
   the second (`group2`) for contrasts (`group1-group2`).  

[⬆ Back to top](#overview)

---

## Stage 2: Running DiffBind Analysis

### Purpose
Use the CSV created in [Stage 1](#stage-1-creating-new-diffbind-csv-files) and DiffBind to identify consensus peaks and normalize 
read counts.  

- **QC mode analysis** produces a TMM-normalized counts table, suitable for use with 
  external differential peak callers or for generating heatmaps. It also calculates PCA 
  and UMAP coordinates.  
- **Differential peak calling analysis** uses DiffBind to identify differential peaks, 
  as done in the chrom-seek pipeline.  

### Input Files Required
- DiffBind CSV file (from [Stage 1](#stage-1-creating-new-diffbind-csv-files))  
- BAM (and associated BAI) files and peak files referenced in the CSV  

### Prepare Run Conditions

1. Request an interactive session (adjust resources if needed):  
   ```bash
   sinteractive -N 1 -n 1 --time=1-12:00:00 --mem=100G --gres=lscratch:200 --cpus-per-task=4
   ```

2. Initialize the Singularity container for DiffBind v2.  

   > **Note**  
   > These conditions ensure lscratch remains accessible and that the Stage 1 CSV file can 
   retain full paths without causing errors.  

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

1. Run `DiffBind_v2_QC.Rmd`.  

   > **Notes**  
   > - Rmarkdown requires full paths; otherwise it looks for `<INPUT_CSV_FILE>` in `<WORKING_DIR>/bin`.  
   > - The `<PEAK_TOOL>` token must be one of: `macsNarrow`, `Genrich`, `macsBroad`, `SEACR`.  

   **Usage**:  
   ```bash
   Rscript -e 'rmarkdown::render("<WORKING_DIR>/bin/DiffBind_v2_QC.Rmd", 
               output_file="<OUTPUT_HTML_FILE>",
               params=list(csvfile="<INPUT_CSV_FILE>", counts_bed="<OUTPUT_BED_FILE>", 
               counts_csv="<OUTPUT_COUNTS_CSV_FILE>", peakcaller="<PEAK_TOOL>"))'
   ```

   **Example**:  
   ```bash
   Rscript -e 'rmarkdown::render("/data/OpenOmics/dev/datasets/outputs/test_homer/bin/DiffBind_v2_QC.Rmd", 
     output_file="/data/OpenOmics/dev/datasets/outputs/test_homer/PeakQC/DB_QC/AllSamples-macsNarrow/AllSamples-macsNarrow_DiffBindQC.html",
     params=list(csvfile="/data/OpenOmics/dev/datasets/outputs/test_homer/PeakQC/DB_QC/AllSamples-macsNarrow/AllSamples-macsNarrow_DiffBind_prep.csv", 
     counts_bed="/data/OpenOmics/dev/datasets/outputs/test_homer/PeakQC/DB_QC/AllSamples-macsNarrow/AllSamples-macsNarrow_DiffBindQC_TMMcounts.bed", 
     counts_csv="/data/OpenOmics/dev/datasets/outputs/test_homer/PeakQC/DB_QC/AllSamples-macsNarrow/AllSamples-macsNarrow_DiffBindQC_TMMcounts.csv", 
     peakcaller="macsNarrow"))'
   ```

### Differential Peak Calling Analysis

1. Load counts using `DiffBind_v2_load.R`.  

   > **Note**  
   > The `<PEAK_TOOL>` token must be one of: `macsNarrow`, `Genrich`, `macsBroad`, `SEACR`.

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

2. Execute differential comparisons using **DiffBind v2.15.2**.

   1) Choose the method: **DESeq2** or **edgeR**. DESeq2 generally suggested. Consider EdgeR 
   when there is large variance in library size, library size is confounding, or there are prior 
   knowledge suggesting that there should be an equivalent number of peaks on both sides of the contrast.  
   2) Decide whether your design requires  **no blocking** or **blocking** (paired or batch effects).  
   3) Select the relevant script (`<ANALYSIS_SCRIPT>`):  
      - `<WORKING_DIR>/bin/DiffBind_v2_Deseq2.Rmd`  
      - `<WORKING_DIR>/bin/DiffBind_v2_Deseq2_block.Rmd`  
      - `<WORKING_DIR>/bin/DiffBind_v2_EdgeR.Rmd`  
      - `<WORKING_DIR>/bin/DiffBind_v2_EdgeR_block.Rmd`  
   4) Define the group contrast from your experimental setup: `"{group1}_vs_{group2}"` → `<INPUT_CONTRASTS>`.  
   5) Render the report:  

      > **Note**  
      > The `<PEAK_TOOL>` token must be one of: `macsNarrow`, `Genrich`, `macsBroad`, `SEACR`.

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
                  counts="<INPUT_PEAK_COUNTS_RDS_FILE>"))'
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
        contrasts="IFN0h_vs_IFN24h", 
        counts="/data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_counts.rds"))'
      ```

[⬆ Back to top](#overview)

---

## Stage 3: Peak Annotation with UROPA

### Purpose
Annotate genomic peaks with gene information, regulatory elements, and genomic context 
using UROPA.

### Input Files Required
- Peak file (BED or narrowPeak format) `<PEAK_FILE>`  
- GTF annotation file `<GTF_FILE>`  
- UROPA configuration file `<UROPA_CONFIG>`  

The GTF file used in the pipeline can be found in `config.json` under `[references][<SPECIES>][GTFFILE]`.

### Steps

1. Create a UROPA configuration file. The example below is what is currently used for chrom-seek when assay is not `cfchip`.
   
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

2. Run UROPA on Biowulf (or in your container). This example uses the `uropa` module. We suggest 4 threads and 10G memory.
   
   **Usage**:
   ```bash
   module load uropa/4.0.2
   uropa -i <UROPA_CONFIG> -t 4
   ```
   
   **Example**:
   ```bash
   module load uropa/4.0.2
   uropa -i sample_config.json -t 4
   ```

[⬆ Back to top](#overview)

---

## Stage 4: Merging Annotations with Differential Results

### Purpose
Combine differential binding results from [Stage 2](#stage-2-running-diffbind-analysis) with gene annotations from [Stage 3](#stage-3-peak-annotation-with-uropa) to 
generate comprehensive results tables.

### Input Files Required
- DiffBind peak list (`<OUTPUT_PEAK_BED_LIST>`) from [Stage 2](#stage-2-running-diffbind-analysis)  
- UROPA results (`<UROPA_FINAL_HITS>`) from [Stage 3](#stage-3-peak-annotation-with-uropa)  

### Steps

1. Merge differential results with annotations. Optionally filter at a stricter FDR or 
fold-change threshold than defaults. Set FDR to 1 for all peaks analyzed.
   
   **Usage**:
   ```bash
   python <WORKING_DIR>/bin/merge_diffbind_uropa.py \
     --diffbind <DIFFBIND_PEAK_LIST_FILE> \
     --uropa <UROPA_FINAL_HITS> \
     --fdr <FDR_THRESHOLD> \
     --fold <FC_THRESHOLD> \
     --output <MERGE_OUTPUT>
   ```
   
   **Example**:
   ```bash
   python /data/OpenOmics/dev/datasets/outputs/test_homer/bin/merge_diffbind_uropa.py \
     --uropa /data/OpenOmics/dev/datasets/outputs/test_homer/UROPA_annotations/DiffBind/IFN0h_vs_IFN24h-macsNarrow-EdgeR/IFN0h_vs_IFN24h_macsNarrow_EdgeR_protTSS_uropa_finalhits.txt \
     --diffbind /data/OpenOmics/dev/datasets/outputs/test_homer/DiffBind/IFN0h_vs_IFN24h-macsNarrow/IFN0h_vs_IFN24h-macsNarrow_Diffbind_EdgeR_peak_list.tab \
     --output /data/OpenOmics/dev/datasets/outputs/test_homer/UROPA_DIFFBIND_TBLS/IFN0h_vs_IFN24h-macsNarrow-EdgeR_protTSS_UROPA_DIFFBIND_JOIN.txt \
     --fdr 0.05 \
     --fold 0
   ```

[⬆ Back to top](#overview)

---

## Getting Help

- [Biowulf user guide](https://hpc.nih.gov/docs/userguide.html) (for `sinteractive`)  
- [DiffBind v2.0.2 documentation](https://www.rdocumentation.org/packages/DiffBind/versions/2.0.2)  
- [DiffBind v1.2.4 vignette (archival)](https://www.bioconductor.org/packages//2.10/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf)  
- [UROPA manual](https://uropa-manual.readthedocs.io/)  
