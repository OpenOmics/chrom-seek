# Common quality-control rules
# Includes the following:
#   - preseq
#   - NRF
#   - rawfastqc
#   - fastqc
#   - fastq_screen
#   - multiQC
rule preseq:
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
            se == "yes" else \
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
            se == "yes" else \
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
        join(workpath,trim_dir,"{name}.R1.trim.fastq.gz") if se == "yes" else \
            expand(join(workpath,trim_dir,"{name}.R{rn}.trim.fastq.gz"),name=samples,rn=[1,2])
    output:
        join(workpath,"FQscreen","{name}.R1.trim_screen.txt") if se == "yes" else \
            expand(join(workpath,"FQscreen","{name}.R{rn}.trim_screen.txt"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen","{name}.R1.trim_screen.png") if se == "yes" else \
            expand(join(workpath,"FQscreen","{name}.R{rn}.trim_screen.png"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen2","{name}.R1.trim_screen.txt") if se == "yes" else \
            expand(join(workpath,"FQscreen2","{name}.R{rn}.trim_screen.txt"),name=samples,rn=[1,2]),
        join(workpath,"FQscreen2","{name}.R1.trim_screen.png") if se == "yes" else \
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

rule multiqc:
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
        # join(workpath,deeptools_dir,"spearman_heatmap.Q5DD.RPGC.pdf"),
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