# Glossary

### Important file name extensions  

 The pipeline generate a lot of output files! Many intermediate output files have a unique name or extension to denote a special meaning. In a pipeline output directory, you may see one or many files contain the same base prefix but a different set of extensions, e.g.: `sorted`, `Q5`, or `Q5DD`.
 
_What do each of these file extension(s) mean?_

  - **`sorted`**: The data has been sorted, but has not been filtered in any way. Note: blacklisted reads are filtered before alignment.
  - **`Q5`**: Reads with a mapping quality below 5 have been filtered.
  - **`Q5DD`**: Low mapping filter and deduplication (removal of PCR duplicates). For paired-end data, this means no fragments with the same exact start and end position occur more than once. For single-end data, we use a negative binomial distribution cutoff filter created by MACS.
  - **`RPGC`**: Reads per genomic content, a method for normalizing to library size.
  - **`FRiP`**: Fraction of reads in peaks. A calculation of the proportion of aligned reads that fall within the peaks called by a particular tool for a given sample.

### Annotation options:

The specific parameters chosen are listed in the associated json and pdf files.

_Here are a list of different annotation options:_

  - **`genes`**: Analyze all genes in the gtf with the most lax Uropa parameters
  - **`prot`**: Analyze only protein-coding genes in the gtf with the most lax Uropa parameters
  - **`protSEC`**: Only protein-coding genes. Using Uropa's multi-step annotation approach, annotate based upon: 1) gene start, 2) gene end, 3) gene center, and 4) anywhere
  - **`protTSS`**: Most popular option, ideal for most projects! Only protein-coding genes, focus around TSS sites.

### File format types

Details about many of these file format can be found on this [UCSC page](https://genome.ucsc.edu/FAQ/FAQformat.html).

_Here is a short description of import file types/formats created by the pipeline:_

  - **`bw`**: Short for bigwig. This is a binary file containing the pile-up patterns of the data along the chromosomes. The data is typically normalized and when viewed can be averaged across different window sizes depending on the size of the region being mapped to screen.
  - **`wig`**: Short for wiggle. This is a non-binary form of a bigwig. There are two flavors of this file type, fixed-step and variable-step, each wth different formatting requirements.
  - **`bed`**: A minimum 3 column file of chromosome, start, and end. The files can be up to 12 columns. Since there is no header to these files, the exact order of columns is standardized across the field. A very multi-purpose file format.
  - **`tagAlign`**: Another name for a simple bed file, but typically each row is a read. Rarely used for peak information. 
  - **`bedgraph`**: A 3 column bed file with a fourth column of score.
  - **`narrowPeak`**: The first 6 columns are the same as the bed file, but the remaining 4 columns contain information about the quality of the peak and the location of the peak summit.
  - **`broadPeak`**: The same as a narrowPeak file, but with the 10th column (referring to the summit), missing.

### Peak callers

Peak callers are use to distinguish biological signal from noise within your dataset.

_Here are a list of peak callers the pipeline uses:_

  - **`macsNarrow`**: The macs2 caller for narrow peaks. The most popular peak calling algorithm; this is typically used in most of the large databases. Can only call peaks between 150bp-10kb. Originally designed to handle peaks with only a single maxima/summit. FDR greatly improved with the addition of an "input" control, but generally still more accurate than most other peak callers out there. https://github.com/macs3-project/MACS/
  - **`macsBroad`**: The macs2 caller for slightly broader peaks. Very similar algorithm to macsNarrow, but sometimes works better than macsNarrow when peaks have more than one maxima/summit. https://github.com/macs3-project/MACS/
  - **`SICER`**: A broad peak caller. Can be really useful for some histone marks. Doesn't work well for extra broad domains like lamins, DNA damage markers, or some repressive marks. Allows for a small amount of dips/gaps between peaks. Window and gap parameters may need to be adjusted to improve calls. https://zanglab.github.io/SICER2/
  - **`Genrich`**: Designed with ATAC-seq data in mind. Can work really well, but not all  collaborators like it as it hasn't been published or reviewed. https://github.com/jsh58/Genrich

### Other important tools we use

_Here are a list of other important tools the pipeline uses:_

  - **`Deeptools`**: Visualizations and QC. Here is a [link](https://deeptools.readthedocs.io/en/develop/index.html) to Deeptools documentation.
  - **`DiffBind version2`**: Differential peak calling. Analysis run with Deseq2 and EdgeR. Here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) to DiffBind's documentation.
  - **`Uropa`**: Peak annotations. Here is a [link](https://uropa-manual.readthedocs.io/introduction.html) to Uropa's documentation.
  - **`MEME suite`**: Motif analysis. We use MEME-ChIP for *de novo* motif calling and AME for known motif calling. Note: the Centrimo subportion of MEME-ChIP will give false results for broad peak calling tools. Here is a [link](https://meme-suite.org/meme/index.html) to MEME suite's documentation.


