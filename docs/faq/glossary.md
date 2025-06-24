# Glossary

### Important file name extensions  

 The pipeline generates a lot of output files! Many intermediate output files have a unique name or extension to denote a special meaning. In a pipeline output directory, you may see one or many files contain the same base prefix but a different set of extensions, e.g.: `sorted`, `Q5`, or `Q5DD`.
 
_What do each of these file extension(s) mean?_

  - **`sorted`**: Indicates that the data has been sorted, but no further filtering has been applied. Note: blacklisted reads are filtered before alignment.
  - **`Q5`**: Denotes that reads with a mapping quality below 5 have been filtered out.
  - **`Q5DD`**: Indicates both low mapping quality filtering and deduplication (removal of PCR duplicates). For paired-end data, this means that no fragments with the same exact start and end position occur more than once. For single-end data, a negative binomial distribution cutoff filter created by MACS is used.
  - **`RPGC`**: Stands for "reads per genomic content," a normalization method based on library size.
  - **`FRiP`**: Represents the "fraction of reads in peaks," calculated as the proportion of aligned reads falling within the peaks called by a specific tool for a given sample.
  
### Annotation options:

The specific parameters chosen are listed in the associated json and pdf files. All Uropa options listed here follow an iterative query approach.  

Currently only protTSS is active. Run parameters for other choices available upon request.  

_Here are a list of different annotation options:_

  - **`genes`**: Analyze all genes in the GTF file using the most lenient parameters from Uropa.
  - **`prot`**: Analyze only protein-coding genes in the GTF file using the most lenient parameters from Uropa.
  - **`protSEC`**: Focuses on protein-coding genes, utilizing Uropa's multi-step annotation approach, and annotates sequentially based on gene start, end, center, and anywhere within the gene.
  - **`protTSS`**: This is the most popular option and is ideal for most projects. It focuses exclusively on protein-coding genes and centers around transcription start sites (TSS).

_protTSS run conditions:_  
  - **`For most assays`**:  
      - Query1: Peak center must be within 3kb upstream or 1kb downstream of the TSS (based upon gene orientation)   
      - Query2: Peak center must be within 10kb of the TSS in either direction  
      - Query3: Peak center must be within 100kb of the TSS in either direction  
  - **`cfChIP`**:  
      - Query1: Peak center must be within 3kb of the TSS in either direction  
      - Query2: Peak center must be within 10kb of the TSS in either direction  
      - Query3: Peak center must be within 100kb of the TSS in either direction  


### File format types

Details about many of these file formats can be found on this [UCSC page](https://genome.ucsc.edu/FAQ/FAQformat.html).

_Here is a short description of import file types/formats created by the pipeline:_

  - **`bw`**: Short for bigwig. Binary file containing normalized pile-up patterns of data along chromosomes, viewable and adjustable across different window sizes.
  - **`wig`**: Short for wiggle. Non-binary format of BigWig, with fixed-step and variable-step variants, each with specific formatting requirements.
  - **`bed`**: Minimum 3-column file (chromosome, start, end), extendable up to 12 columns, used for various purposes due to standardized column order. Lacks a header.
  - **`tagAlign`**: Simple bed file format, typically with each row representing a read, seldom used for peak information
  - **`bedgraph`**: 3-column bed file with an additional score column.
  - **`narrowPeak`**: Similar to bed file with additional columns for peak quality and summit location.
  - **`broadPeak`**: Similar to narrowPeak but lacks the summit location column.

### Peak callers

Peak callers are use to distinguish biological signal from noise within your dataset.

_Here are a list of peak callers the pipeline uses:_

  - **`macsNarrow`**: The macs2 caller optimized for narrow peaks, widely recognized as the most popular peak calling algorithm. Typically used in large databases, it identifies peaks within the range of 150bp to 10kb. Originally designed to handle peaks with a single maxima/summit, its false discovery rate (FDR) has been greatly improved with the addition of an "input" control. It is generally more accurate than most other peak callers, even without controls. https://github.com/macs3-project/MACS/
  - **`macsBroad`**: The macs2 caller for slightly broader peaks, sharing a similar algorithm with macsNarrow. It is particularly useful when peaks exhibit more than one maxima/summit. https://github.com/macs3-project/MACS/
  - **`SICER`**: SICER is a broad peak caller that can be highly effective for certain histone marks. However, it may not perform well for extra broad domains such as lamins or some repressive marks. It allows for a small amount of gaps between peaks, and users may need to adjust window and gap parameters for optimal results. https://zanglab.github.io/SICER2/
  - **`Genrich`**: Designed with ATAC-seq data in mind, Genrich can yield excellent results. However, it may not be universally favored by all collaborators due to its lack of formal publication or review. https://github.com/jsh58/Genrich

### Other important tools we use

_Here are a list of other important tools the pipeline uses:_

  - **`Deeptools`**: This tool is employed for visualizations and quality control (QC) purposes. You can find the documentation for Deeptools at [link](https://deeptools.readthedocs.io/en/develop/index.html).
  - **`DiffBind version2`**: Used for conducting differential peak calling analyses, this tool integrates with Deseq2 and EdgeR for analysis. Here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) to DiffBind's documentation.
  - **`Uropa`**: Uropa is utilized for peak annotations, providing comprehensive annotation features. Here is a [link](https://uropa-manual.readthedocs.io/introduction.html) to Uropa's documentation.
  - **`MEME suite`**: Employed for motif analysis, the MEME suite includes MEME-ChIP for *de novo* motif discovery and AME for known motif analysis. Note that the Centrimo subcomponent of MEME-ChIP may produce inaccurate results for broad peak calling tools. Here is a [link](https://meme-suite.org/meme/index.html) to MEME suite's documentation.
