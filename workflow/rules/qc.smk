# Common quality-control rules
# Includes the following:
#   - preseq
#   - NRF
#   - rawfastqc
#   - fastqc
#   - fastq_screen
#   - multiQC
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
        bam = join(workpath,bam_dir,"{name}.sorted.bam"),
    output:
        ccurve = join(workpath,qc_dir,"{name}.ccurve"),
    params:
        rname = "preseq",
        preseqver=config['tools']['PRESEQVER'],
    shell: """
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
        bam=join(workpath,bam_dir,"{name}.sorted.bam"),
    output:
        preseq=join(workpath,qc_dir,"{name}.preseq.dat"),
        preseqlog=join(workpath,qc_dir,"{name}.preseq.log"),
        nrf=temp(join(workpath,qc_dir,"{name}.nrf")),
    params:
        rname='NRF',
        samtoolsver=config['tools']['SAMTOOLSVER'],
        rver=config['tools']['RVER'],
        preseqver=config['tools']['PRESEQVER'],
        nrfscript=join(workpath,"workflow","scripts","atac_nrf.py "),
    threads: 16
    shell: """
    module load {params.preseqver};
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
        expand(join(workpath,"{name}.R1.fastq.gz"), name=samples) if \
            not paired_end else \
            expand(join(workpath,"{name}.R{rn}.fastq.gz"), name=samples,rn=[1,2])
    output:
        expand(join(workpath,'rawfastQC',"{name}.R1_fastqc.html"),name=samples),
    params:
        rname='rawfastqc',
        outdir=join(workpath,"rawfastQC"),
    envmodules: 
        config['tools']['FASTQCVER']
    threads:
        int(allocated("threads", "rawfastqc", cluster))
    shell: """
    fastqc \\
        {input} \\
        -t {threads} \\
        -o {params.outdir}
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
        expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),name=samples) if \
            not paired_end else \
            expand(join(workpath,trim_dir,"{name}.R{rn}.trim.fastq.gz"), name=samples,rn=[1,2])
    output:
        expand(join(workpath,'fastQC',"{name}.R1.trim_fastqc.html"),name=samples),
    params:
        rname='fastqc',
        outdir=join(workpath,"fastQC"),
    envmodules: 
        config['tools']['FASTQCVER']
    threads:
        int(allocated("threads", "fastqc", cluster))
    shell: """
    fastqc \\
        {input} \\
        -t {threads} \\
        -o {params.outdir}
    """

rule fastq_screen:
    """
    Quality-control step to screen for different sources of contamination.
    FastQ Screen compares your sequencing data to a set of different reference
    genomes to determine if there is contamination. It allows a user to see if
    the composition of your library matches what you expect.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        FastQ Screen report and logfiles
    """
    input:
        join(workpath,trim_dir,"{name}.R1.trim.fastq.gz") if not paired_end else \
            expand(join(workpath,trim_dir,"{name}.R{rn}.trim.fastq.gz"),name=samples,rn=[1,2])
    output:
        join(workpath,"FQscreen","{name}.R1.trim_screen.txt") if not paired_end else \
            expand(join(workpath,"FQscreen","{name}.R{rn}.trim_screen.txt"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen","{name}.R1.trim_screen.png") if not paired_end else \
            expand(join(workpath,"FQscreen","{name}.R{rn}.trim_screen.png"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen2","{name}.R1.trim_screen.txt") if not paired_end else \
            expand(join(workpath,"FQscreen2","{name}.R{rn}.trim_screen.txt"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen2","{name}.R1.trim_screen.png") if not paired_end else \
            expand(join(workpath,"FQscreen2","{name}.R{rn}.trim_screen.png"),name=samples,rn=[1,2]),
    params:
        rname   = 'fqscreen',
        outdir  = join(workpath,"FQscreen"),
        outdir2 = join(workpath,"FQscreen2"),
        # Exposed Parameters: modify resources/fastq_screen{_2}.conf 
        # to change defaults locations to bowtie2 indices
        fastq_screen         = config['bin']['FASTQ_SCREEN'],
        fastq_screen_config1 = config['shared_resources']['FASTQ_SCREEN_CONFIG_P1'],
        fastq_screen_config2 = config['shared_resources']['FASTQ_SCREEN_CONFIG_P2'],
    envmodules:
        config['tools']['BOWTIE2VER'],
        config['tools']['PERLVER'],
    threads: 
        int(allocated("threads", "fastq_screen", cluster))
    shell: """
    # First pass of contamination screening
    {params.fastq_screen} \\
        --conf {params.fastq_screen_config1} \\
        --outdir {params.outdir} \\
        --threads {threads} \\
        --subset 1000000 \\
        --aligner bowtie2 \\
        --force \\
        {input}
    # Second pass of contamination screening
    {params.fastq_screen} \\
        --conf {params.fastq_screen_config2} \\
        --outdir {params.outdir2} \\
        --threads {threads} \\
        --subset 1000000 \\
        --aligner bowtie2 \\
        --force \\
        {input}
    """

rule kraken:
    """
    Quality-control step to assess for potential sources of microbial contamination.
    If there are high levels of microbial contamination, Kraken will provide an
    estimation of the taxonomic composition. Kraken is used in conjunction with
    Krona to produce an interactive reports.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Kraken logfile and interative krona report
    """
    input:
        fq1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        fq2=provided(join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"), paired_end)
    output:
        krakenout = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.out.txt"),
        krakentaxa = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.taxa.txt"),
        kronahtml = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.krona.html"),
    params:
        rname='kraken',
        outdir=join(workpath,kraken_dir),
        bacdb=config['shared_resources']['KRAKENBACDB'],
        tmpdir=tmpdir,
        paired_end = paired_end
    threads: int(allocated("threads", "kraken_pe", cluster)),
    envmodules:
        config['tools']['KRAKENVER'],
        config['tools']['KRONATOOLSVER'],
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    
    # Copy kraken2 db to /lscratch or temp 
    # location to reduce filesystem strain
    cp -rv {params.bacdb} ${{tmp}}/;
    kdb_base=$(basename {params.bacdb})
    if [ '{params.paired_end}' == True ]; then
        kraken2 --db ${{tmp}}/${{kdb_base}} \\
            --threads {threads} --report {output.krakentaxa} \\
            --output {output.krakenout} \\
            --gzip-compressed \\
            --paired {input.fq1} {input.fq2}
    else
        kraken2 --db ${{tmp}}/${{kdb_base}} \\
            --threads {threads} --report {output.krakentaxa} \\
            --output {output.krakenout} \\
            --gzip-compressed \\
            {input.fq1}
    fi
    
    # Generate Krona Report
    cut -f2,3 {output.krakenout} | \\
        ktImportTaxonomy - -o {output.kronahtml}
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
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.krona.html"),name=samples),
        expand(join(workpath,qc_dir,"{name}.ccurve"), name=samples),
        expand(join(workpath,bam_dir,"{name}.Q5DD.bam.flagstat"), name=samples),
        expand(join(workpath,bam_dir,"{name}.Q5.bam.flagstat"), name=samples),
        # join(workpath,qc_dir,"QCTable.txt"),
        expand(join(workpath,"rawfastQC","{name}.R1_fastqc.html"),name=samples),
        expand(join(workpath,"fastQC","{name}.R1.trim_fastqc.html"),name=samples),
        # expand(join(workpath,deeptools_dir,"{group}.fingerprint.raw.Q5DD.tab"),group=groups),
        join(workpath,deeptools_dir,"spearman_heatmap.Q5DD_mqc.png")
    output:
        join(workpath,"multiqc_report.html")
    params:
        rname="multiqc",
        multiqc=config['tools']['MULTIQCVER'],
	    qcconfig=join(workpath, config['shared_resources']['MULTIQC_CONFIG']),
	    excludedir=join(workpath,extra_fingerprint_dir),
        dir=workpath
    shell: """
    module load {params.multiqc}
    multiqc \\
        -f \\
        -c {params.qcconfig} \\
        --interactive \\
        -e cutadapt \\
        --ignore {params.excludedir} \\
        -d {params.dir}
    """

rule insert_size:
    """
    Quality step calculates number of reads per insert size.
    @Input:
        Sorted only bam file, and also bam files that were sorted, 
        filtered by mapQ a value, and deduplicated (extensions: sorted and Q5DD),
        for all samples.
    @Output:
        Number of reads per insert size and their histogram
    """
    input:
        bam = lambda w : join(workpath,bam_dir,w.name + "." + w.ext + "." + extensionsDict[w.ext])
    output:
        txt= join(workpath,qc_dir,"{name}.{ext}.insert_size_metrics.txt"),
        pdf= join(workpath,qc_dir,"{name}.{ext}.insert_size_histogram.pdf"),
    params:
        rname="insert_size",
        picardver=config['tools']['PICARDVER'],
        rver=config['tools']['RVER'],
        javaram='16g',
    shell: """
    module load {params.picardver} {params.rver};
    java -Xmx{params.javaram} -jar ${{PICARDJARPATH}}/picard.jar CollectInsertSizeMetrics \\
        -I {input.bam} \\
        -O {output.txt} \\
        -H {output.pdf}
    """

rule deeptools_QC:
    input:
        [ join(workpath, bw_dir, name + ".Q5DD.RPGC.bw") for name in samples ] # this should be all bigwigs
    output:
        heatmap=join(workpath,deeptools_dir,"spearman_heatmap.Q5DD.pdf"),
        pca=join(workpath,deeptools_dir,"pca.Q5DD.pdf"),
	    npz=temp(join(workpath,deeptools_dir,"Q5DD.npz")),
	    png=join(workpath,deeptools_dir,"spearman_heatmap.Q5DD_mqc.png")
    params:
        rname="deeptools_QC",
        deeptoolsver=config['tools']['DEEPTOOLSVER'],
        labels=samples # this should be the sample names to match the bigwigs in the same order
    shell: """    
    module load {params.deeptoolsver}
    multiBigwigSummary bins -b {input} -l {params.labels} -out {output.npz}
    plotCorrelation -in {output.npz} -o {output.heatmap} -c 'spearman' -p 'heatmap' --skipZeros --removeOutliers
    plotCorrelation -in {output.npz} -o {output.png} -c 'spearman' -p 'heatmap' --skipZeros --removeOutliers
    plotPCA -in {output.npz} -o {output.pca}
    """

