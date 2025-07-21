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

1. Grab the relevant CSV from within ./PeakQC/DB_QC  
2. Remove rows with samples that should not be included in the analysis.  
3. Update information in ONLY the following three columns:  
      a. CONDITION:  
          This column contains the group information. This is the column
       DiffBind will use to define the contrast when running the differential binding
       analysis scripts. Do not start the values with numbers.  
      b. TREATMENT:  
          Feel free to add this column as needed. The differential binding analysis
       scripts use this column for blocking/paired analysis when available. Do not start the values with numbers.  
      c. REPLICATE:  
          This is a simple counting device for the tool and nothing more.
4. For differential analysis, you will want to sort your rows so that the top rows are group1 for the contrast group1-group2.

## Workflow 2: Running DiffBind Analysis

### Purpose
Perform either quality control analysis or differential binding analysis using DiffBind.

### Input Files Required
- DiffBind CSV file (from Workflow 1)
- BAM files and peak files referenced in the CSV

### QC Mode Analysis

1. **Run DiffBind QC analysis:**

    a. Get an interactive session with 100G mem, 200G lscratch, and 4 threads.
     **[Needs filled in]**  
    b. Initialize the docker containing Diffbind v2.  
   **[Needs filled in; make sure it still can use lscratch]**  
    c. Run ./bin/DiffBind_v2_QC.Rmd using example code.  
    **[Needs filled in; Example code from line 521 of dba.smk]**  
`Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", counts_bed="{output.countsbed}", 
            counts_csv="{output.countscsv}", peakcaller="{params.peak_tool}"))'`

### Differential Peak Calling Mode ###

    a. Get an interactive session with 100G mem, 200G lscratch, and 4 threads.
     **[Needs filled in]**
    b. Initialize the docker containing Diffbind v2.  
   **[Needs filled in; make sure it still can use lscratch]**  
   c. Run ./bin/DiffBind_v2_load.R  
**[Needs filled in; Example code in rule diffbind_count, begins line 50 of dba.smk]**  
   d. Run one of the four Rmds:  
       - ./bin/DiffBind_v2_Deseq2.Rmd  
       - ./bin/DiffBind_v2_Deseq2_block.Rmd  
       - ./bin/DiffBind_v2_EdgeR.Rmd  
       - ./bin/DiffBind_v2_EdgeR_block.Rmd  
**[Needs filled in; Example code in rule diffbind_edger, begins line 218 of dba.smk]**

## Workflow 3: Peak Annotation with UROPA

### Purpose
Annotate genomic peaks with gene information, regulatory elements, and genomic context using UROPA.

### Input Files Required
- Peak files (BED or narrowPeak format)
- GTF annotation file
- UROPA configuration file

### Steps

1. **Create UROPA configuration file:**
   ```bash
   ./bin/create_uropa_config.py \
     --gtf resources/genomes/hg38/genes.gtf \
     --output uropa_config.json \
     --distance 3000
   ```

2. **Run UROPA annotation:**
   ```bash
   # For single peak file
   python ./bin/run_uropa_annotation.py \
     --peaks results/peak_calling/Sample1_peaks.narrowPeak \
     --config uropa_config.json \
     --output uropa_annotations/Sample1_annotated.txt

   # For multiple peak files (batch processing)
   ./bin/batch_uropa_annotation.sh \
     --peak-dir results/peak_calling \
     --config uropa_config.json \
     --output-dir uropa_annotations \
     --cores 4
   ```

3. **Expected outputs:**
   - `uropa_annotations/Sample1_annotated.txt` - Detailed annotations
   - `uropa_annotations/Sample1_summary.txt` - Annotation summary statistics
   - `uropa_annotations/annotation_plots.pdf` - Visualization plots

### Configuration Options
The UROPA config file supports various annotation priorities:
- **Promoter regions**: TSS ± 3kb
- **Gene body**: Exons and introns
- **Intergenic regions**: Regions between genes
- **Custom features**: User-defined genomic regions

### Example Configuration
```json
{
  "queries": [
    {
      "feature": "gene",
      "feature.anchor": "start",
      "distance": 3000,
      "strand": "same"
    },
    {
      "feature": "gene",
      "feature.anchor": "start",
      "distance": 10000,
      "strand": "same"
    }
  ],
  "priority": "True",
  "gtf": "resources/genomes/hg38/genes.gtf",
  "bed": "input.bed"
}
```

## Workflow 4: Merging Annotations with Differential Results

### Purpose
Combine differential binding results from DiffBind with gene annotations from UROPA to create comprehensive results tables.

### Input Files Required
- Differential peaks results from Workflow 2
- UROPA annotation results from Workflow 3

### Steps

1. **Merge differential results with annotations:**
   ```bash
   python ./bin/merge_diffbind_uropa.py \
     --diffbind-results diffbind_differential_results/differential_peaks.csv \
     --uropa-annotations uropa_annotations/ \
     --output merged_results.tsv \
     --include-plots
   ```

2. **Generate summary report:**
   ```bash
   Rscript ./bin/generate_summary_report.R \
     --merged-results merged_results.tsv \
     --output-dir final_results \
     --create-plots TRUE
   ```

3. **Expected outputs:**
   - `merged_results.tsv` - Combined differential and annotation data
   - `final_results/summary_report.html` - Interactive HTML report
   - `final_results/enrichment_plots.pdf` - Functional enrichment plots
   - `final_results/top_genes_table.csv` - Table of top differentially bound genes

### Output Columns in Merged Results
- **Peak information**: chr, start, end, peak_name
- **Differential binding**: log2FoldChange, pvalue, padj, baseMean
- **Annotation**: gene_name, gene_id, feature_type, distance_to_tss
- **Genomic context**: promoter, exon, intron, intergenic
- **Additional**: peak_width, summit_position, fold_enrichment

## Advanced Analysis Options

### Custom Filtering and Subsetting

```bash
# Filter results by fold change and significance
python ./bin/filter_results.py \
  --input merged_results.tsv \
  --output filtered_results.tsv \
  --min-lfc 2 \
  --max-padj 0.01 \
  --feature-types promoter,exon

# Create subsets by genomic features
python ./bin/subset_by_features.py \
  --input merged_results.tsv \
  --output-dir feature_subsets \
  --group-by feature_type
```

### Pathway Enrichment Analysis

```bash
# Run pathway enrichment on differentially bound genes
Rscript ./bin/pathway_enrichment.R \
  --gene-list final_results/top_genes_table.csv \
  --output-dir pathway_analysis \
  --organism human \
  --databases GO,KEGG,Reactome
```

### Visualization and Plotting

```bash
# Create custom visualizations
python ./bin/create_custom_plots.py \
  --results merged_results.tsv \
  --output-dir custom_plots \
  --plot-types volcano,ma,heatmap,genomic_distribution
```

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
