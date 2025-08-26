# Quality control metrics and reports
# ~~~~
# Generally applicable quality control rules
from scripts.common import allocated
from textwrap import dedent
from os.path import join


# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
paired_end                      = False if config['project']['nends'] == 1 else True
samples                         = config['samples']
ends                            = [1] if not paired_end else [1, 2]
assay                           = config['options']['assay']
PeakTools                       = get_peaktools(assay)
chips                           = config['project']['peaks']['chips']

# ~~ directories
qc_dir                          = join(workpath, "QC")
kraken_dir                      = join(workpath, 'kraken')
deeptools_dir                   = join(workpath, 'deeptools')
peakqc_dir                      = join(workpath, "PeakQC")
extra_fingerprint_dir           = join(deeptools_dir, 'sorted_fingerprint')
tmpdir                          = config['options']['tmp_dir']


rule preseq:
    """
    Quality step to estimate library complexity. Low library complexity may indicate
    an issue with library preparation or sample storage (FFPE samples) where very
    little input RNA was over-amplified or the sample may be highly degraded.
    @Input:
        Sorted, duplicate marked genomic BAM file (scatter)
    @Output:
        Logfile containing library complexity information
    """
    input:
        bam                     = join(bam_dir, "{name}.sorted.bam"),
    output:
        ccurve                  = join(qc_dir, "{name}.ccurve"),
    params:
        rname                   = "preseq",
        preseqver               = config['tools']['PRESEQVER'],
    shell: 
        """
        module load {params.preseqver};
        preseq c_curve \\
            -B \\
            -o {output.ccurve} \\
            {input.bam}            
        """


rule NRF:
    """
    Quality step computes the expected yield for theoretical 
    larger experiments and the associated confidence intervals
    @Input:
        Sorted BAM file (scatter)
    @Output:
        Output is a text file with four columns. The total number 
        of reads, average expected number of distinct reads, and
        the lower and upper limits of the confidence interval.
        And pyhton code produces NRF = distinct_reads/tot_reads,
        PBC1 = one_pair/distinct_reads, and PBC2 = one_pair/two_pair.
    """
    input:
        bam                     = join(bam_dir, "{name}.sorted.bam"),
    output:
        preseq                  = join(qc_dir, "{name}.preseq.dat"),
        preseqlog               = join(qc_dir, "{name}.preseq.log"),
        nrf                     = join(qc_dir, "{name}.nrf"),
    params:
        rname                   = 'NRF',
        nrfscript               = join(bin_path, "nrf.py"),
    threads: 16
    container: config['images']['preseq']
    shell: 
        """
        preseq lc_extrap \\
            -P \\
            -B \\
            -D \\
            -o {output.preseq} \\
            {input.bam} \\
            -seed 12345 \\
            -v \\
            -l 100000000000 \\
        2> {output.preseqlog}
        python {params.nrfscript} \\
            {output.preseqlog} \\
        > {output.nrf}
        """


rule rawfastqc:
    """
    Quality-control step to assess sequencing quality of the raw data prior removing
    adapter sequences. FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Raw FastQ files (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """
    input:
        expand(join(workpath, "{name}.R{rn}.fastq.gz"), name=samples, rn=list(map(str, ends)))
    output:
        expand(join(qc_dir, "rawfastQC", "{name}.R{rn}_fastqc.html"), name=samples, rn=ends),
    params:
        rname                   = 'rawfastqc',
        outdir                  = join(qc_dir, "rawfastQC"),
        tmpdir                  = tmpdir,
    envmodules: 
        config['tools']['FASTQCVER']
    threads:
        int(allocated("threads", "rawfastqc", cluster))
    shell: 
        """
        # fastqc storage on lscratch b/c nfs bug
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        fastqc \\
            {input} \\
            -t {threads} \\
            -o "${{tmp}}"
            
        find "${{tmp}}" \\
            -type f \\
            \\( -name '*.html' -o -name '*.zip' \\) \\
            -exec cp {{}} {params.outdir} \\; 
        """


rule fastqc:
    """
    Quality-control step to assess sequencing quality of the raw data after removing
    adapter sequences. This step is run after trim_pe rule. FastQC is run after adapter
    trimming to evalute if the adapter sequences were properly removed.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Trimmed FastQC reports and zip file containing data quality information
    """
    input:
        expand(join(trim_dir, "{name}.R{rn}.trim.fastq.gz"), name=samples, rn=ends)
    output:
        expand(join(qc_dir, 'fastQC', "{name}.R{rn}.trim_fastqc.html"), name=samples, rn=ends),
    params:
        rname                   = 'fastqc',
        outdir                  = join(qc_dir, "fastQC"),
        tmpdir                  = tmpdir,
    envmodules: 
        config['tools']['FASTQCVER']
    threads:
        int(allocated("threads", "fastqc", cluster))
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        # Running fastqc with local
        # disk or a tmpdir, fastqc
        # has been observed to lock
        # up gpfs filesystems, adding
        # this on request by HPC staff
        fastqc \\
            {input} \\
            -t {threads} \\
            -o "${{tmp}}"
        
        # Copy output files from tmpdir
        # to output directory
        find "${{tmp}}" \\
            -type f \\
            \\( -name '*.html' -o -name '*.zip' \\) \\
            -exec cp {{}} {params.outdir} \\;
        """


rule multiqc:
    """
    Reporting step to aggregate sample statistics and quality-control information
    across all samples. This will be one of the last steps of the pipeline. The inputs
    listed here are to ensure that this step runs last. During runtime, MultiQC will
    recurively crawl through the working directory and parse files that it supports.
    @Input:
        List of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report and a QC metadata table
    """
    input: 
        expand(join(qc_dir, "FQscreen", "{name}.R1.trim_screen.txt"), name=samples),
        expand(join(qc_dir, "FQscreen2", "{name}.R1.trim_screen.txt"), name=samples),
        expand(join(qc_dir, "{name}.ccurve"), name=samples),
        expand(join(qc_dir, "rawfastQC", "{name}.R1_fastqc.html"), name=samples),
        expand(join(qc_dir, "fastQC", "{name}.R1.trim_fastqc.html"), name=samples),
        expand(join(kraken_dir, "{name}.trim.kraken_bacteria.taxa.txt"), name=samples),
        expand(join(bam_dir, "{name}.Q5DD.bam.flagstat"), name=samples),
        expand(join(bam_dir, "{name}.Q5.bam.flagstat"), name=samples),
        join(deeptools_dir, "spearman_readcounts.Q5DD.tab"),
        join(deeptools_dir, "fingerprint.raw.Q5DD.tab"),
	join(deeptools_dir,"TSS_profile.Q5DD.tab")
    output:
        join(workpath, "multiqc_report.html")
    params:
        rname                   = "multiqc",
        multiqc                 = config['tools']['MULTIQCVER'],
	    qcconfig                = join(workpath, config['shared_resources']['MULTIQC_CONFIG']),
	    excludedir              = join(workpath, extra_fingerprint_dir),
    shell: 
        """
        module load {params.multiqc}
        multiqc \\
            -f \\
            -c {params.qcconfig} \\
            --interactive \\
            -e cutadapt \\
            --ignore {params.excludedir} \\
            -d """ + workpath + """
        """


rule deeptools_QC:
    input:
        [ join(bw_dir, name + ".Q5DD.RPGC.bw") for name in samples ] 
    output:
        heatmap                 = join(deeptools_dir, "spearman_heatmap.Q5DD.pdf"),
        pca                     = join(deeptools_dir, "pca.Q5DD.pdf"),
        npz                     = temp(join(deeptools_dir, "Q5DD.npz")),
	    mqc                     = join(deeptools_dir, "spearman_readcounts.Q5DD.tab"),
	    png                     = join(deeptools_dir, "spearman_heatmap.Q5DD_mqc.png")
    params:
        rname                   = "deeptools_QC",
        parent_dir              = deeptools_dir,
        deeptoolsver            = config['tools']['DEEPTOOLSVER'],
        # this should be the sample names to match the bigwigs in the same order
        labels                  = samples
    threads: 24
    shell: 
        """    
        module load {params.deeptoolsver}
        if [ ! -d "{params.parent_dir}" ]; then mkdir "{params.parent_dir}"; fi
        multiBigwigSummary bins -b {input} -p {threads} -l {params.labels} -out {output.npz}
        plotCorrelation -in {output.npz} -o {output.heatmap} -c 'spearman' -p 'heatmap' \\
               --skipZeros --removeOutliers --outFileCorMatrix {output.mqc}
        plotCorrelation -in {output.npz} -o {output.png} -c 'spearman' -p 'heatmap' --skipZeros --removeOutliers
        plotPCA -in {output.npz} -o {output.pca}
        """


rule deeptools_fingerprint:
    input:
        [ join(bam_dir, name + ".Q5DD.bam") for name in samples ] 
    output:
        image=join(deeptools_dir, "fingerprint.Q5DD.pdf"),
        raw=temp(join(deeptools_dir, "fingerprint.raw.Q5DD.tab")),
        metrics=join(deeptools_dir, "fingerprint.metrics.Q5DD.tsv"),
    params:
        rname                   = "deeptools_fingerprint",
        parent_dir              = deeptools_dir,
        deeptoolsver            = config['tools']['DEEPTOOLSVER'],
        # this should be the sample names to match the bigwigs in the same order
        labels                  = samples,
        ext                     = "" if paired_end else "-e 200",
    threads: int(allocated("threads", "deeptools_fingerprint", cluster)),
    shell: 
        dedent("""    
        module load {params.deeptoolsver}
        if [ ! -d "{params.parent_dir}" ]; then mkdir "{params.parent_dir}"; fi
        plotFingerprint \\
         -b {input} \\
         --labels {params.labels} \\
         -p {threads} \\
         --skipZeros \\
         --outQualityMetrics {output.metrics} \\
         --plotFile {output.image} \\
         --outRawCounts {output.raw} {params.ext}
        """)


rule deeptools_gene_all:
    input:
        [ join(bw_dir, name + ".Q5DD.RPGC.bw") for name in samples ] 
    output:
        bed                     = temp(join(deeptools_dir, "geneinfo.Q5DD.bed")),
        TSSline                 = join(deeptools_dir, "TSS_profile.Q5DD.pdf"),
        TSSmat                  = temp(join(deeptools_dir, "TSS.Q5DD.mat.gz")),
        TSSheat                 = join(deeptools_dir, "TSS_heatmap.Q5DD.pdf"),
	    mqc                     = join(deeptools_dir, "TSS_profile.Q5DD.tab"),
        metamat                 = temp(join(deeptools_dir, "metagene.Q5DD.mat.gz")),
        metaline                = join(deeptools_dir, "meta_profile.Q5DD.pdf"),
        metaheat                = join(deeptools_dir, "metagene_heatmap.Q5DD.pdf"),
    params:
        rname                   = "deeptools_gene_all",
        parent_dir              = deeptools_dir,
        deeptoolsver            = config['tools']['DEEPTOOLSVER'],
        labels                  = samples,
        prebed                  = config['references'][genome]['GENEINFO'],
    threads: 4
    # eventually threads should be 16
    shell: 
        dedent("""    
        module load {params.deeptoolsver}
        if [ ! -d "{params.parent_dir}" ]; then mkdir "{params.parent_dir}"; fi
        grep --line-buffered 'protein_coding' {params.prebed} | awk -v OFS='\t' -F'\t' '{{print $1, $2, $3, $5, ".", $4}}' > {output.bed}
        # TSS
        computeMatrix reference-point \\
            -S {input} \\
            -R {output.bed} \\
            -p {threads} \\
            --referencePoint TSS \\
            --upstream 3000 \\
            --downstream 3000 \\
            --skipZeros \\
            -o {output.TSSmat} \\
            --samplesLabel {params.labels}
        plotProfile \\
            -m {output.TSSmat} \\
            -out {output.TSSline} \\
            --yAxisLabel 'average RPGC' \\
            --plotType 'se' \\
            --numPlotsPerRow 5 \\
            --outFileNameData {output.mqc}
        plotHeatmap -m {output.TSSmat} \\
            -out {output.TSSheat} \\
            --colorMap 'BuPu' \\
            --yAxisLabel 'average RPGC' \\
            --regionsLabel 'genes'
        # metagene
        computeMatrix scale-regions \\
            -S {input} \\
            -R {output.bed} \\
            -p {threads} \\
            --upstream 1000 \\
            --regionBodyLength 2000 \\
            --downstream 1000 \\
            --skipZeros \\
            -o {output.metamat} \\
            --samplesLabel {params.labels}
        plotHeatmap \\
            -m {output.metamat} \\
            -out {output.metaheat} \\
            --colorMap 'BuGn' \\
            --yAxisLabel 'average RPGC' \\
            --regionsLabel 'genes'
        plotProfile \\
            -m {output.metamat} \\
            -out {output.metaline} \\
            --yAxisLabel 'average RPGC' \\
            --plotType 'se' \\
            --numPlotsPerRow 5
        """)


rule enhancer_plot:
    input:
        [ join(bw_dir, name + ".Q5DD.RPGC.bw") for name in samples ]
    output:
        heatmap                 = join(deeptools_dir, "enhancer_heatmap.Q5DD.pdf"),
        matrix                  = join(deeptools_dir, "enhancer_matrix.Q5DD.tsv"),
        line                    = join(deeptools_dir, "enhancer_profile.Q5DD.pdf")
    params:
        rname                   = "enhancer_plot",
        enhancer_ref            = enhancer_ref,
        parent_dir              = deeptools_dir,
        tmpdir                  = tmpdir,
        labels                  = samples,
        deeptoolsver            = config['tools']['DEEPTOOLSVER']
    threads: int(allocated("threads", "enhancer_plot", cluster))
    shell: 
        dedent("""
        module load {params.deeptoolsver}
        if [ ! -d "{params.parent_dir}" ]; then mkdir "{params.parent_dir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        
        computeMatrix reference-point \\
            -S {input} \\
            -R {params.enhancer_ref} \\
            -p {threads} \\
            --referencePoint 'center' \\
            --upstream 3000 \\
            --downstream 3000 \\
            --skipZeros \\
            -o {output.matrix} \\
            --samplesLabel {params.labels}
        plotProfile \\
            -m {output.matrix} \\
            -out {output.line} \\
            --yAxisLabel 'average RPGC' \\
            --plotType 'se' \\
            --numPlotsPerRow 5
        plotHeatmap -m {output.matrix} \\
            -out {output.heatmap} \\
            --yAxisLabel 'average RPGC' \\
            --regionsLabel 'enhancers' \\
            --colorList '#FFFFE5,#004529'
        """)


rule FRiP_macsN:
    input:
        peaks                   = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        tbl                     = join(peakqc_dir, "FRiP", "macsNarrow", "macsNarrow.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_macsN",
        tblscript               = join(bin_path, "frip.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir,
        this_config             = join(workpath, 'config.json')
    container: 
        config['images']['python']
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.tblscript} \\
            -p {input.peaks} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {output.tbl} \\
            -x 16
        """


rule FRiP_Genrich:
    input:
        peaks                   = join(genrich_dir, "{name}", "{name}.narrowPeak"),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        tbl                     = join(peakqc_dir, "FRiP", "Genrich", "Genrich.{name}.Q5DD.FRiP_table.txt")
    params:
        rname                   = "FRiP_Genrich",
        tblscript               = join(bin_path, "frip.py"),
        plotscript              = join(bin_path, "FRiP_plot.R"),
        this_config             = join(workpath, 'config.json'),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir,
    container: 
        config['images']['python']
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.tblscript} \\
            -p {input.peaks} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {output.tbl} \\
            -x 16
        """


rule FRiP_macsB:
    input:
        peaks                   = expand(join(macsB_dir, "{name}", "{name}_peaks.broadPeak"), name=chips),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        tbl                     = join(peakqc_dir, "FRiP", "macsBroad", "macsBroad.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_macsB",
        tblscript               = join(bin_path, "frip.py"),
        plotscript              = join(bin_path, "FRiP_plot.R"),
        this_config             = join(workpath, 'config.json'),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir,
    container: 
        config['images']['python']
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.tblscript} \\
            -p {input.peaks} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {output.tbl} \\
            -x 16
        """


rule FRiP_SEACR:
    input:
        peaks                   = join(seacr_dir, "{name}", "{name}.stringent.bed"),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        tbl                     = join(peakqc_dir, "FRiP", "SEACR", "SEACR.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_SEACR",
        tblscript               = join(bin_path, "frip.py"),
        plotscript              = join(bin_path, "FRiP_plot.R"),
        this_config             = join(workpath, 'config.json'),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir,
    container: 
        config['images']['python']
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.script} \\
            -p {input.peaks} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {output.tbl} \\
            -x 16
        """


rule frip_summary:
    input:
        tbl                     = expand(join(peakqc_dir, "FRiP", "{pktool}", "{pktool}.{name}.Q5DD.FRiP_table.txt"), pktool=PeakTools, name=samples)
    output:
        heatmap                 = join(peakqc_dir, "FRiP", "summary.Q5DD.FRiP_heatmap.pdf"),
        bar_plot                = join(peakqc_dir, "FRiP", "summary.Q5DD.FRiP_barplot.pdf"),
        # scatter_plot            = join(peakqc_dir, "FRiP", "summary.Q5DD.FRiP_scatter.pdf")
    params:
        rname                   = "frip_summary",
        plotscript              = join(bin_path, "frip_summary_plot.py"),
        this_config             = join(workpath, 'config.json')
    shell:
        dedent("""
        {params.plotscript} \\
            --hm {output.heatmap} \\
            -t {input.tbl} \\
            -b {output.bar_plot} \\
            -c {params.this_config} \\
            -v
        """)


rule jaccard_genrich:
    input:
        expand(join(genrich_dir, "{SID}", "{SID}.narrowPeak"), SID=chips),
    output:
        table                   = join(peakqc_dir, "jaccard", 'Genrich_jaccard.txt'),
        pcaplot                 = join(peakqc_dir, "jaccard", 'Genrich_jaccard_pca.pdf'),
        pcatab                  = join(peakqc_dir, "jaccard", 'Genrich_jaccard_pca.tsv'),
        heatmap                 = join(peakqc_dir, "jaccard", 'Genrich_jaccard_heatmap.pdf'),
        heatmaptab              = join(peakqc_dir, "jaccard", 'Genrich_jaccard_heatmap.tsv'),
    params:
        rname                   = "jaccard_genrich",
        outroot                 = join(peakqc_dir, "jaccard"),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \\
            -i "{input}" \\
            --caller Genrich \\
            --pcatab {output.pcatab} \\
            --outtable {output.table} \\
            --pcaplot {output.pcaplot} \\
            --outheatmap {output.heatmap} \\
            --tabheatmap {output.heatmaptab} \\
            -g {params.genome}
        """)


rule jaccard_macsbroad:
    input:
        expand(join(macsB_dir, "{SID}", "{SID}_peaks.broadPeak"), SID=chips),
    output:
        table                   = join(peakqc_dir, "jaccard", 'macsBroad_jaccard.txt'),
        pcaplot                 = join(peakqc_dir, "jaccard", 'macsBroad_jaccard_pca.pdf'),
        pcatab                  = join(peakqc_dir, "jaccard", 'macsBroad_jaccard_pca.tsv'),
        heatmap                 = join(peakqc_dir, "jaccard", 'macsBroad_jaccard_heatmap.pdf'),
        heatmaptab              = join(peakqc_dir, "jaccard", 'macsBroad_jaccard_heatmap.tsv'),
    params:
        rname                   = "jaccard_macsbroad",
        outroot                 = join(peakqc_dir, "jaccard"),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \\
            -i "{input}" \\
            --caller macsBroad \\
            --pcatab {output.pcatab} \\
            --outtable {output.table} \\
            --pcaplot {output.pcaplot} \\
            --outheatmap {output.heatmap} \\
            --tabheatmap {output.heatmaptab} \\
            -g {params.genome}
        """)


rule jaccard_macsnarrow:
    input:
        expand(join(macsN_dir, "{SID}", "{SID}_peaks.narrowPeak"), SID=chips),
    output:
        table                   = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard.txt'),
        pcaplot                 = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard_pca.pdf'),
        pcatab                  = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard_pca.tsv'),
        heatmap                 = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard_heatmap.pdf'),
        heatmaptab              = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard_heatmap.tsv'),
    params:
        rname                   = "jaccard_macsnarrow",
        outroot                 = join(peakqc_dir, "jaccard"),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \\
            -i "{input}" \\
            --caller macsNarrow \\
            --pcatab {output.pcatab} \\
            --outtable {output.table} \\
            --pcaplot {output.pcaplot} \\
            --outheatmap {output.heatmap} \\
            --tabheatmap {output.heatmaptab} \\
            -g {params.genome}
        """)


rule jaccard_seacr:
    input:
        expand(join(seacr_dir, "{SID}", "{SID}.stringent.bed"), SID=chips),
    output:
        table                   = join(peakqc_dir, "jaccard", 'SEACR_jaccard.txt'),
        pcaplot                 = join(peakqc_dir, "jaccard", 'SEACR_jaccard_pca.pdf'),
        pcatab                  = join(peakqc_dir, "jaccard", 'SEACR_jaccard_pca.tsv'),
        heatmap                 = join(peakqc_dir, "jaccard", 'SEACR_jaccard_heatmap.pdf'),
        heatmaptab              = join(peakqc_dir, "jaccard", 'SEACR_jaccard_heatmap.tsv')
    params:
        rname                   = "jaccard_seacr",
        outroot                 = join(peakqc_dir, "jaccard"),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \\
            -i "{input}" \\
            --caller SEARC \\
            --pcatab {output.pcatab} \\
            --outtable {output.table} \\
            --pcaplot {output.pcaplot} \\
            --outheatmap {output.heatmap} \\
            --tabheatmap {output.heatmaptab} \\
            -g {params.genome}
        """)


rule jaccard_summary:
    input:
        pca_cords               = expand(join(peakqc_dir, "jaccard", '{pkcaller}_jaccard_pca.tsv'), pkcaller=PeakTools),
        hm_cords                = expand(join(peakqc_dir, "jaccard", '{pkcaller}_jaccard_heatmap.tsv'), pkcaller=PeakTools),
    output:
        pcaplot                 = join(peakqc_dir, "jaccard", 'jaccard_summary_pca.pdf'),
        hmplot                 = join(peakqc_dir, "jaccard", 'jaccard_summary_heatmap.pdf'),
    params:
        rname                   = "jaccard_summary",
        genome                  = config['references'][genome]['REFLEN'],
        script                  = join(bin_path, "jaccard_summary.py"),
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \
            --pca {input.pca_cords} \
            --hm {input.hm_cords}
        """)


rule summary_tbls:
    input:
        # aligns to glob(f"{workpath}/**/*flagstat", recursive=True) in summary_tables.py
        flagstat_sorted             = expand(join(bam_dir, "{name}.sorted.bam.flagstat"), name=samples),
        flagstat_q5                 = expand(join(bam_dir, "{name}.Q5.bam.flagstat"), name=samples),
        flagstat_q5dd               = expand(join(bam_dir, "{name}.Q5DD.bam.flagstat"), name=samples),
        # aligns to glob(f"{workpath}/**/*Q5DD.insert_size_metrics.txt") in summary_tables.py
        insert_size                 = expand(join(qc_dir, "{name}.Q5DD.insert_size_metrics.txt"), name=samples),
        # glob(f"{workpath}/PeakQC/FRiP/**/*.FRiP_table.txt") in summary_tables.py
        frip                        = expand(join(peakqc_dir, "FRiP", "{pktool}", "{pktool}.{name}.Q5DD.FRiP_table.txt"), pktool=PeakTools, name=samples)
    output:
        mapping                     = join(workpath, "MappingSummary.txt"),
        encode                      = join(workpath, "EncodeQC.txt") if paired_end else [],
        peak                        = join(workpath, "PeakSummary.txt"),
    container: 
        config['images']['python']
    params:
        rname                   = "summary_tbls",
        script                  = join(bin_path, "summary_tables.py"),
        workpath                = workpath,
        ended_flag              = "--paired" if paired_end else "", 
        tmpdir                  = tmpdir
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} {params.ended_flag} {params.workpath}
        """)
