
# Tools and versions

## Quality-control assessment tools

The chrom-seek pipeline runs a series of quality-control (QC) tools to assess the sequencing and library quality of each sample. The tools used in the pipeline are listed below, along with their versions and a brief description of their purpose.

| **Tool**                | **Version**  | **Notes**                                                                      |
| ----------------------- | :----------: | ------------------------------------------------------------------------------ |
| FastQC<sup>1<sup>       | 0.11.9       | Assess sequencing quality, run before and after adapter trimming               |
| Kraken<sup>2<sup>       | 2.1.2        | Assess microbial taxonomic composition                                         |
| KronaTools<sup>3<sup>   | 2.8.1        | Visualize kraken output                                                        |
| FastQ Screen<sup>4<sup> | 0.9.3        | Assess contamination; additional dependencies: `bowtie2/2.5.1`, `perl/5.36`    |
| Preseq<sup>5<sup>       | 3.1.2        | Estimate library complexity                                                    |
| MultiQC<sup>6<sup>      | 1.14         | Aggregate sample statistics and quality-control information across all samples |

**References**

<sup>1. **FastQC:** Andrews, S. (2010). FastQC: a quality control tool for high
throughput sequence data.https://www.bioinformatics.babraham.ac.uk/projects/fastqc</sup>  
<sup>2. **Kraken:** Wood, D. E. and S. L. Salzberg (2014). "Kraken: ultrafast
metagenomic sequence classification using exact alignments." Genome
Biol 15(3): R46. http://ccb.jhu.edu/software/kraken/</sup>  
<sup>3. **Krona:** Ondov, B. D., et al. (2011). "Interactive metagenomic
visualization in a Web browser." BMC Bioinformatics 12(1): 385.
https://github.com/marbl/Krona/wiki</sup>  
<sup>4. **FastQ Screen:** Wingett, S. and S. Andrews (2018). "FastQ Screen: A
tool for multi-genome mapping and quality control." F1000Research 7(2):
1338. https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/</sup>  
<sup>5. **Preseq:** Daley, T. and A.D. Smith (2013). Predicting the molecular
complexity of sequencing libraries. Nat Methods 10(4): 325-7.
http://smithlabresearch.org/software/preseq/</sup>  
<sup>6. **MultiQC:** Ewels, P., et al. (2016). "MultiQC: summarize analysis
results for multiple tools and samples in a single report."
Bioinformatics 32(19): 3047-3048. https://multiqc.info/docs/</sup>  

## Data processing tools

The pipeline is composed of a series of data processing steps that include adapter trimming, read alignment, duplicate removal, and bigwig creation. The tools used in the pipeline are listed below, along with their versions and a brief description of their purpose.


| **Tool**                | **Version** | **Notes**                                              |
| ----------------------- | :---------: | :----------------------------------------------------- |
| Cutadapt<sup>1</sup>    | 4.4         | Remove adapter sequences and perform quality trimming  |
| BWA mem<sup>2</sup>     | 0.7.17      | Read alignment, first to identify reads aligning to blacklisted regions and later for the remainder of the genome      |
| Picard<sup>3</sup>      | 2.27.3      | Run SamToFastq (for blacklist read removal) and MarkDuplicates (to remove PCR duplicates in PE data)                   |
| SAMtools<sup>4</sup>    | 1.17        | Remove reads with mapQ less than 6. Also run flagstat and idxstats to calculate alignment statistics.                  |
| MACS<sup>5</sup>        | 2.2.7.1     | Run filterdup on SE data (`--keep-dup="auto"`) to remove PCR duplicates                                                |
| Bedtools<sup>6</sup>    | 2.27.1      | Run intersect and bedtobam to convert `.tag.Align.gz` to `.bam` for use with deeptools (specific to SE data) and MEME  |
| ppqt<sup>7,8</sup>      | 1.1.2       | Also known as phantompeakqualtools, used to calculate estimated fragment length (used for bigwig and peak calling for SE data). Also produces QC metrics: NSC and RSC. |
| deepTools<sup>9</sup>   | 3.5.1       | Used for bigwig creation and multiple QC metrics. Use bamcoverage to create RPGC-normalized data: `--binSize 25 --smoothLength 75 --normalizeUsing RPGC`. For PE data, add `--centerReads`. For control SE, add `-e 200`. For ChIP SE, add `-e [estimated fragmentlength]`. For control subtraction (inputnorm), use bigwigCompare: `--binSize 25 --operation 'subtract'`. Run multiBigWigSummary, plotCorrelation, plotPCA, plotFingerprint, computeMatrix, plotHeatmap, and plotProfile for QC plots. Note: not all these have been reincorporated into the pipeline. |

**References**

<sup>1. **Cutadapt:** Martin, M. (2011). "Cutadapt removes adapter sequences
from high-throughput sequencing reads." EMBnet 17(1): 10-12.
https://cutadapt.readthedocs.io/en/stable/ </sup>  
<sup>2. **BWA:** Li H. and Durbin R. (2009) Fast and accurate short read
alignment with Burrows-Wheeler Transform. Bioinformatics 25: 1754-60.
http://bio-bwa.sourceforge.net/bwa.shtml </sup>  
<sup>3. **Picard:** The Picard toolkit.
https://broadinstitute.github.io/picard/</sup>  
<sup>4. **SAMtools:** Danecek, Petr et al. “Twelve years of SAMtools and BCFtools.” 
GigaScience vol. 10,2 (2021): giab008. </sup>  
<sup>5. **MACS:** Zhang, Y., et al. (2008). Model-based Analysis of ChIP-Seq
(MACS). Genome Biol 9: R137. https://github.com/macs3-project/MACS</sup>  
<sup>6. **Bedtools:** Quinlan, A.R. (2014). BEDTools: The Swiss‐Army Tool for
Genome Feature Analysis. Current Protocols in Bioinformatics, 47:
11.12.1-11.12.34. https://bedtools.readthedocs.io/en/latest/index.html</sup>  
<sup>7. **Ppqt:** Landt S.G., et al. (2012). ChIP-seq guidelines and practices
of the ENCODE and modENCODE consortia. Genome Res 22(9): 1813-31.
https://github.com/kundajelab/phantompeakqualtools</sup>  
<sup>8. **Ppqt:** Kharchenko P.K., et al. (2008). Design and analysis of
ChIP-seq experiments for DNA-binding proteins Nat Biotechnol 26(12):
1351-9.</sup>  
<sup>9. **deepTools:** Ramírez, F., et al. (2016). deepTools2: A next Generation
Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Research
44(W1): W160--W165. https://deeptools.readthedocs.io/en/develop/</sup>  

## Peak calling and differential peak calling tools

The pipeline includes peak calling and differential peak calling tools to identify enriched regions of interest in the genome. The tools used in the pipeline are listed below, along with their versions and a brief description of their purpose.

| **Tool**                | **Version** | **Notes**                                              |
| ----------------------- | :---------: | :----------------------------------------------------- |
| MACS<sup>1</sup>        | 2.2.7.1     | **`macsNarrow`**: The macs2 caller optimized for narrow peaks, widely recognized as the most popular peak calling algorithm. Typically used in large databases, it identifies peaks within the range of 150bp to 10kb. Originally designed to handle peaks with a single maxima/summit, its false discovery rate (FDR) has been greatly improved with the addition of an "input" control. It is generally more accurate than most other peak callers, even without controls. **`macsBroad`**: The macs2 caller for slightly broader peaks, sharing a similar algorithm with macsNarrow. It is particularly useful when peaks exhibit more than one maxima/summit. |
| Sicer<sup>2</sup>       | 2-1.0.3     | Sicer is a broad peak caller that can be highly effective for certain histone marks. However, it may not perform well for extra broad domains such as lamins or some repressive marks. It allows for a small amount of gaps between peaks, and users may need to adjust window and gap parameters for optimal results. |
| Genrich<sup>3</sup>     | 0.6         | Designed with ATAC-seq data in mind, Genrich can yield excellent results. However, it may not be universally favored by all collaborators due to its lack of formal publication or review. |
| SEACR<sup>4</sup>       | 1.3         | Peak caller for CUT & RUN or CUT & TAG only; can utilize IgG as if it was an input control, but assumes IgG is sparse. |
| MANorm<sup>5</sup>      | 1.1.4       | Differential peak calling when no replicates. This tool has not been incorporated into the pipeline. |
| DiffBind<sup>6,7</sup>  | 2.15.2      | Used for conducting differential peak calling analyses, this tool integrates with Deseq2 and EdgeR for analysis. Here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) to DiffBind's documentation. |

**References**

<sup>1. **MACS:** Zhang, Y., et al. (2008). Model-based Analysis of ChIP-Seq
(MACS). Genome Biol 9: R137. https://github.com/macs3-project/MACS</sup>  
<sup>2. **Sicer:** Xu S., et al. (2014). Spatial Clustering for Identification
of ChIP-Enriched Regions (SICER) to Map Regions of Histone Methylation
Patterns in Embryonic Stem Cells. Methods Mol Biol 1150: 97--111.
https://zanglab.github.io/SICER2/</sup>  
<sup>3. **Genrich:** Gaspar,J.M. (2018) Genrich: Detecting sites of genomic
enrichment. https://github.com/jsh58/Genrich.</sup>  
<sup>4. **SEACR:** Meers MP, Tenenbaum D, Henikoff S. (2019). Peak calling by Sparse 
Enrichment Analysis for CUT&RUN chromatin profiling. Epigenetics and Chromatin 12(1):42. 
https://github.com/FredHutch/SEACR</sup>  
<sup>5. **MANorm:** Shao, Z., et al. (2012). MAnorm: a robust model for
quantitative comparison of ChIP-Seq data sets. Genome Biology 13: R16.
https://manorm.readthedocs.io/en/latest/index.html</sup>  
<sup>6. **DiffBind:** Ross-Innes C.S., et al. (2012). Differential oestrogen
receptor binding is associated with clinical outcome in breast cancer.
Nature 481: 389--393.</sup>  
<sup>7. **DiffBind:** Stark R. and G. Brown. (2011). DiffBind: differential
binding analysis of ChIP-Seq peak data.
http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf</sup>  
<sup>8. **Deseq2:** Love, Michael I et al. “Moderated estimation of fold change and dispersion 
for RNA-seq data with DESeq2.” Genome biology vol. 15,12 (2014): 550.</sup>  


## Annotations, motifs, and QC metrics

The pipeline includes tools for peak annotation, motif calling, and quality-control metrics. The tools used in the pipeline are listed below, along with their versions and a brief description of their purpose.

| **Tool**                | **Version** | **Notes**                                              |
| ----------------------- | :---------: | :----------------------------------------------------- |
| Uropa<sup>1</sup>       | 4.0.3       | Uropa is utilized for peak annotations, providing comprehensive annotation features. Here is a [link](https://uropa-manual.readthedocs.io/introduction.html) to Uropa's documentation. See the glossary for options in this pipeline: [link](https://openomics.github.io/chrom-seek/faq/glossary/#annotation-options). |
| Homer<sup>2</sup>       | 4.11.1      | Homer is being used for motif calling                  |
| MEME<sup>3</sup>        | 5.5.5       | Employed for motif analysis, the MEME suite includes MEME-ChIP for *de novo* motif discovery and AME for known motif analysis. Note that the Centrimo subcomponent of MEME-ChIP may produce inaccurate results for broad peak calling tools. Here is a [link](https://meme-suite.org/meme/index.html) to MEME suite's documentation. |
| IDR<sup>4</sup>         | 2.0.3       |  One method for identifying consensus peaks. Only works for 2 replicates. This tool is not currently in the pipeline. |
| Jaccard                 | NA          | Calculation of peak call consistency between two conditions. This is currently is not included in the pipeline. Requires: `pybedtools` `pysam`. |
| FRiP                    | NA          | Represents the _"fraction of reads in peaks"_ calculated as the proportion of aligned reads falling within the peaks called by specific tool for a given sample. Requires: `pybedtools` `pysam` |

**References**

<sup>1. **Uropa:** Kondili M., et al. (2017). UROPA: a tool for Universal RObust
Peak Annotation. Scientific Reports 7: 2593.</sup>  
<sup>2. **Homer:** Heinz S., et al. (2010). Simple combinations of
lineage-determining transcription factors prime cis-regulatory elements
required for macrophage and B cell identities. Mol Cell 38(4): 576--589.</sup>  
<sup>3. **MEME:** Bailey T.L., et al. (2015). The MEME Suite. Nucleic Acids Res
43(W1):W39-W49.</sup>  
<sup>4. **IDR:** Li Q., et al. (2011). Measuring reproducibility of
high-throughput experiments. Ann Appl Stat 5(3): 1752-1779.</sup>  


## Reference genome resources
**Blacklist citation**: The ENCODE Project Consortium., Moore, J.E., Purcaro, M.J. et al. 
Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 583, 699–710 (2020).

| **Reference**  | **Version**       | **Website**                                            |
| -------------- | :---------------: | :----------------------------------------------------- |
| hg19           | GrCh37.p13     | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz |
|                | Gencode v19    | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz |
|                | Blacklist v1   | https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz |
|                | Enhancer reference   | https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz |
| hg38           | GRCh38.p12     | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz |
|                | Gencode v28    | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz |
|                | Blacklist v2   | http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz |
|                | Enhancer reference   | https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/extra/CAGE_peaks/hg38_fair_CAGE_peaks_combined_phase1+2.bed.gz |
| mm10           | GRCm38.p6      | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.p6.genome.fa.gz |
|                | Gencode vM18   | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.annotation.gtf.gz |
|                | Blacklist v2   | http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz |
|                | Enhancer reference   | https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_v9/extra/CAGE_peaks/mm10_fair+new_CAGE_peaks_phase1and2.bed.gz |
| rheMac10       | Ensembl 110    | Chr1 example: https://ftp.ensembl.org/pub/release-110/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.1.fa.gz |
|                | Ensembl 110    | https://ftp.ensembl.org/pub/release-110/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.110.chr.gtf.gz  |

Note: All fastas and gtfs in built into this pipeline only include reference chromosomes, ie autosomes and sex chromosomes.  
See: https://github.com/Boyle-Lab/Blacklist/tree/master/lists for all v1 and v2 blacklists available.

## Main differences between assays
- **ChIP**
   - Only assay that can handle SE data  
   - Peak callers are: macs Narrow and macs Broad  
- **cfChIP**
   - Peak caller is: macs Narrow  
   - Only runs on hg19
   - Filters chromosome X
   - Changes to annotation conditions
   - UMAP plot created in the DiffBindQC report  
   - Some rules run in the other versions do not run here
- **ATAC**
   - Peak callers are: macs Narrow and Genrich  
- **cutnrun**  
   - Peak callers are: macs Narrow and SEACR  

## Acknowledgements

### Biowulf

If you [utilized NIH's Biowulf cluster](https://hpc.nih.gov/Research/) to run `chrom-seek`, *please do not forget to provide an acknowlegement*! The continued growth and support of NIH's Biowulf cluster is dependent upon its demonstrable value to the NIH Intramural Research Program. If you publish research that involved significant use of Biowulf, please cite the cluster.

**Suggested citation text**

```
This work utilized the computational resources of the NIH HPC Biowulf cluster. (http://hpc.nih.gov)
```
