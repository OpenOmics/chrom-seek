# Best practices

## Our take-home message

1. Always remember to consider the specific immunoprecipation or methodology for each analysis.
2. Bigwigs serve as our closest approximation to truth sets. It's essential to cross-check results against them.
3. The boundaries of peaks are primarily influenced by the tool used for analysis rather than being inherent biological features.
4. This type of data often contains significant noise and a high background level.
5. The majority of reads are not detected within peaks. Among the remaining reads, some may represent background artifacts, while others have biological significance. However, distinguishing between the two is challenging. Additionally, the proportion of reads within peaks is influenced by biological factors.

## Key points to consider

**Not all ChIP-seq data can be treated equally!**

_Here are some main classes of proteins to consider when analyzing ChIP-seq data:_  

- Transcription factors
- Promoter-associated histones  
- Enhancer-associated histones  
- Repressive histones  
- RNA Polymerase II  
- Repressive histone marks  
- Chromosome structure proteins  
- DNA damage proteins  

Transcription factors are often regarded as the most straightforward targets in ChIP-seq experiments due to their high abundance, strong binding affinity, and well-characterized binding motifs. They typically exhibit clear and distinct peaks in ChIP-seq data, facilitating their detection and analysis. Moreover, transcription factors are known for their resilience to variations in processing and analysis protocols compared to other genomic features. They also demonstrate reduced dependence on ideal experimental controls. Additionally, transcription factors are frequently associated with established gold standards, which aid in sensitivity and specificity assessments in ChIP-seq assays. Their specific binding to DNA sequences further simplifies validation in laboratory experiments.

Databases like [ENCODE](https://www.encodeproject.org/) and [Cistrome](http://cistrome.org/) are valuable resources for studying specific marks. These databases employ robust data processing methods, although the sensitivity and specificity of processed features may be relatively lower compared to other methods. Some features may be poorly identified, particularly in a subset of cases.

High noise and background levels affect all these methods, influenced by both biological factors and sample processing. Considering this influence is crucial at every stage of experimental design and data processing.

ATAC-seq presents two primary challenges. Firstly, it explores two types of open chromatin regions, each characterized by distinct patterns and functions. These regions encompass histones flanked by spaces, typically found within active genes, and larger open regions devoid of stable histones, such as promoters. Secondly, ATAC-seq lacks a comprehensive global truth set for validation. Consequently, there is a scarcity of literature on best practices for ATAC-seq, with methodologies often relying on knowledge from similar methods.

_When annotating peaks, it's essential to remember:_  

- A single peak may influence multiple genes.
- Unless it occurs within a gene's promoter or within 1kb of a transcriptional start site, determining if the nearest gene is actually affected by the peak is challenging without additional information.
- Activating transcription factors within gene bodies are unlikely to affect that gene directly.
- Peak boundaries are arbitrary.
- While summits are likely binding sites for transcription factors, they may not hold true for other features.

## Recommendations

### General considerations

Using QC metrics as the sole criterion for sample quality assessment is unreliable. While QC values may be informative for certain features or laboratories, they do not universally apply. While they provide a starting point, they should **not** be used as an automatic decision-making tool for sample exclusion!

Normalized counts matrices and derived results should not be solely relied upon, as they depend on the quality of peak calls, consensus peak identification, and normalization methods used. It's crucial to confirm the validity of results using normalized bigwigs or directly generated plots.  

_Bigwigs offer enhanced visibility as they have been processed relative to the bam files in several ways:_  

- Normalization to library size, utilizing only useful reads.
- Position adjustment as needed, particularly significant for single-end ChIP.
- Normalization against input to mitigate one form of biological noise, if available.
- 75 bp smoothing, aligning with the resolution of most methods (50-250bp).

### Peak calling

Here's a structured approach for choosing a peak caller that is most suitable for your data. Remember to consider the specific immunoprecipation or methodology for each analysis. The following guidelines are based on the most common types of data, but they may not be applicable to all scenarios. Always consult the literature, consider the specific characteristics of your data, and visually inspect the peak calls along with their relevent tracks before making a descision.

1. _Transcription Factor (TF) Analysis_
     - Inspect the macsNarrow calls in IGV.
          - If macsNarrow calls appear satisfactory, use macsNarrow.
          - If not, consider using macsBroad.
2. _ATAC-seq Analysis_
     - Determine if publication status matters.
        - If publication status matters, utilize macsNarrow.
        - If not, compare results from macsNarrow and Genrich to determine the better option.
3. _CUN & RUN or CUT & TAG Analysis_
     - Compare results from macsNarrow and SEACR to determine the better option.
4. _Narrow Histone Mark Analysis_
     - Uncertain about whether the data represents a narrow histone mark? Consult the [ENCODE documentation](https://www.encodeproject.org/chip-seq/histone/).
     - Compare results from macsNarrow and macsBroad.
         - If one outperforms the other consistently, choose the better option.
         - If not, either option is viable.
5. _Other Marks Analysis_
     - Check the macsBroad calls. 
         - If they appear suitable, use macsBroad.
     - If macsBroad calls don't meet expectations, assess the pipeline default sicer calls in the chrome-seek pipeline.
         - If sicer calls are acceptable, utilize those.
         - If not, adjust sicer conditions as necessary, noting that peak calling may not be feasible for all marks.

### Differential analyses

Establishing a robust method for identifying consensus peaks is essential. To effectively utilize any differential analysis tool, it's imperative to have a shared set of peaks across all samples. This step facilitates the removal of peaks with low counts, so avoid filtering outside of this process. Keep in mind that peak calling introduces noise, and peak boundaries lack biological significance. Additionally, be aware that peaks in one sample may overlap with multiple peaks in another sample.

Normalization plays a crucial role in differential peak calling. Numerous studies indicate that both the consensus peak list and normalization method significantly impact results. However, determining the optimal approach can be challenging. A consistent finding in the literature is that normalizing solely to reads within peaks, without external information, consistently yields subpar results. This is likely because a significant proportion of reads (60-99%) often reside outside peaks, representing biologically relevant but rarer binding sites. I prefer scaling to the number of reads aligning to the genome following deduplication. In DiffBind v2, the TMM output takes the raw counts, corrects for inputs when applicable, and employs a DESeq2 function to scale the samples to an external value (the number of reads in the bam file) rather than relying solely on the internal values in the counts table.

In ChIP-seq and ATAC-seq analyses, simpler models tend to yield better results. This is likely due to the various sources of noise in the data, including biological, technical, and computational factors. It's essential to scrutinize results from any method used to ensure that background noise isn't being erroneously identified.

In the majority of projects, skews in global fold-changes towards one condition or the other are likely to be biologically meaningful. It's essential not to artificially equalize the number of peaks on both sides of the contrast unless there's prior knowledge or a substantial difference in the size of the original libraries (exceeding 10X).

It's not uncommon to observe small differential fold-changes, even among the most significant peaks. Many papers assessing differential peak calls, as well as the default conditions in tools like DiffBind, often do not include fold-change filtering. However, in my experience, it can be advantageous to incorporate a low fold-change cutoff. A reasonable starting point is an absolute fold-change cutoff of 1.5 (equivalent to a log2 value of 0.585). Alternatively, using a cutoff of 4 (log2 value of 2) will typically remove almost all peaks unless the analysis is focused on transcription factors.

**Please note:** The pipeline opts to use a depreciated version of DiffBind for several reasons. This is by design!

_Here is why:_

1. DiffBind v2 has consistently demonstrated superior performance compared to other differential peak callers across a wide range of features documented in the literature.
2. DiffBind v3 lacks rigorous testing, and the author does not provide rationale for changes in default counting or normalization parameters. It is feasible to achieve results similar to DiffBind v2, especially for ATAC data or in cases lacking input samples, by configuring the pipeline correctly. However, it consistently treats input samples differently regardless of configuration.
3. The tool automatically generates a set of consensus peaks by combining peaks identified in at least two samples, regardless of condition, allowing for partial overlap. This results in each contrast having a unique set of comparison peaks, but it enhances the accuracy of differential peak calling. However, adjustments may be necessary for large cohorts with numerous samples.
4. In DiffBind v2, the default normalization method scales the entire library to all reads in the alignment files. This approach is generally considered one of the most effective normalization options available.

## Useful links

Here are a set of useful links for learning more about ChIP-seq and ATAC-seq.

### ENCODE standards and recommendations

- [Transcription factors](https://www.encodeproject.org/chip-seq/transcription_factor/)
- [Histone marks](https://www.encodeproject.org/chip-seq/histone/)
- [ATAC-seq](https://www.encodeproject.org/atac-seq/)
