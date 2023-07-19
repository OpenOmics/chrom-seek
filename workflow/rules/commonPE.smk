# Common paired-end rules:
#   - trim_pe
#   - kraken_pe
#   - BWA_PE
#   - insert_size

# trim, remove PolyX and remove blacklist reads
rule trim_pe:
    input:
        file1=join(workpath,"{name}.R1.fastq.gz"),
        file2=join(workpath,"{name}.R2.fastq.gz"),
    output:
        outfq1=temp(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz")),
        outfq2=temp(join(workpath,trim_dir,"{name}.R2.trim.fastq.gz")),
    params:
        rname="trim",
        cutadaptver=config['tools']['CUTADAPTVER'],
        workpath=config['project']['workpath'],
        fastawithadaptersetd=join(workpath, config['shared_resources']['ADAPTERS_FASTA']),
        blacklistbwaindex=config['references'][genome]['BLACKLISTBWAINDEX'],
        picardver=config['tools']['PICARDVER'],
        bwaver=config['tools']['BWAVER'],
        samtoolsver=config['tools']['SAMTOOLSVER'],
        minlen=35,
        leadingquality=10,
        trailingquality=10,
        javaram="64g",
        sample="{name}"
    threads: 16
    shell: """
    module load {params.cutadaptver};
    if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
    cd /lscratch/$SLURM_JOBID
    cutadapt \\
        --pair-filter=any \\
        --nextseq-trim=2 \\
        --trim-n \\
        -n 5 \\
        -O 5 \\
        -q {params.leadingquality},{params.trailingquality} \\
        -m {params.minlen}:{params.minlen} \\
        -b file:{params.fastawithadaptersetd} \\
        -B file:{params.fastawithadaptersetd} \\
        -j {threads} \\
        -o {params.sample}.R1.cutadapt.fastq \\
        -p {params.sample}.R2.cutadapt.fastq \\
        {input.file1} {input.file2}
    
    module load {params.bwaver};
    module load {params.samtoolsver};
    module load {params.picardver};
    bwa mem -t {threads} \\
        {params.blacklistbwaindex} \\
        {params.sample}.R1.cutadapt.fastq \\
        {params.sample}.R2.cutadapt.fastq \\
        | samtools view -@{threads} \\
            -f4 \\
            -b \\
            -o {params.sample}.bam
    
    java -Xmx{params.javaram} -jar $PICARDJARPATH/picard.jar SamToFastq \\
        VALIDATION_STRINGENCY=SILENT \\
        INPUT={params.sample}.bam \\
        FASTQ={params.sample}.R1.cutadapt.noBL.fastq \\
        SECOND_END_FASTQ={params.sample}.R2.cutadapt.noBL.fastq \\
        UNPAIRED_FASTQ={params.sample}.unpaired.noBL.fastq
    
    pigz -p {threads} {params.sample}.R1.cutadapt.noBL.fastq;
    pigz -p {threads} {params.sample}.R2.cutadapt.noBL.fastq;
    
    mv /lscratch/$SLURM_JOBID/{params.sample}.R1.cutadapt.noBL.fastq.gz {output.outfq1};
    mv /lscratch/$SLURM_JOBID/{params.sample}.R2.cutadapt.noBL.fastq.gz {output.outfq2};
    """


rule kraken_pe:
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
        fq2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
    output:
        krakenout = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.out.txt"),
        krakentaxa = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.taxa.txt"),
        kronahtml = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.krona.html"),
    params:
        rname='kraken',
        outdir=join(workpath,kraken_dir),
        bacdb=config['shared_resources']['KRAKENBACDB'],
        tmpdir='/lscratch/$SLURM_JOBID',
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
    kraken2 --db ${{tmp}}/${{kdb_base}} \\
        --threads {threads} --report {output.krakentaxa} \\
        --output {output.krakenout} \\
        --gzip-compressed \\
        --paired {input.fq1} {input.fq2}
    
    # Generate Krona Report
    cut -f2,3 {output.krakenout} | \\
        ktImportTaxonomy - -o {output.kronahtml}
    """


rule BWA_PE:
    input:
        infq1 = join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        infq2 = join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
    params:
        d=join(workpath,bam_dir),
        rname='bwa',
        reference=config['references'][genome]['BWA'],
        bwaver=config['tools']['BWAVER'],
        samtoolsver=config['tools']['SAMTOOLSVER'],
        script=join(workpath,"workflow","scripts","bam_filter_by_mapq.py"),
        pythonver="python/3.5"
    output:
        outbam1=join(workpath,bam_dir,"{name}.sorted.bam"), 
        outbam2=temp(join(join(workpath,bam_dir,"{name}.Q5.bam"))),
        flagstat1=join(workpath,bam_dir,"{name}.sorted.bam.flagstat"),
        idxstat1=join(workpath,bam_dir,"{name}.sorted.bam.idxstat"),
        flagstat2=join(workpath,bam_dir,"{name}.Q5.bam.flagstat"),
        idxstat2=join(workpath,bam_dir,"{name}.Q5.bam.idxstat"),
    threads: 32
    shell: """
    module load {params.bwaver};
    module load {params.samtoolsver};
    module load {params.pythonver};
    cd /lscratch/$SLURM_JOBID;
    bwa mem -t {threads} {params.reference} {input.infq1} {input.infq2} \\
        | samtools sort -@{threads} -o {output.outbam1}
    
    samtools index {output.outbam1}
    samtools flagstat {output.outbam1} > {output.flagstat1}
    samtools idxstats {output.outbam1} > {output.idxstat1}
    #samtools view -b -q 6 {output.outbam1} -o {output.outbam2}
    
    python {params.script} -i {output.outbam1} -o {output.outbam2} -q 6
    samtools index {output.outbam2}
    samtools flagstat {output.outbam2} > {output.flagstat2}
    samtools idxstats {output.outbam2} > {output.idxstat2}
    """

rule picard_dedup:
    input: 
        bam2=join(workpath,bam_dir,"{name}.Q5.bam")
    output:
        out5=join(workpath,bam_dir,"{name}.Q5DD.bam"),
        out5f=join(workpath,bam_dir,"{name}.Q5DD.bam.flagstat"),
        out5i=join(workpath,bam_dir,"{name}.Q5DD.bam.idxstat"),
        out6=join(workpath,bam_dir,"{name}.bwa.Q5.duplic") 
    params:
        rname='dedup',
        picardver=config['tools']['PICARDVER'],
        samtoolsver=config['tools']['SAMTOOLSVER'],
        javaram='16g'
    shell: """
    module load {params.samtoolsver};
    module load {params.picardver}; 
    if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
    cd /lscratch/$SLURM_JOBID
    
    java -Xmx{params.javaram} \\
      -jar $PICARDJARPATH/picard.jar MarkDuplicates \\
      INPUT={input.bam2} \\
      OUTPUT={output.out5} \\
      TMP_DIR=/lscratch/$SLURM_JOBID \\
      VALIDATION_STRINGENCY=SILENT \\
      REMOVE_DUPLICATES=true \\
      METRICS_FILE={output.out6}
    
    
    samtools index {output.out5}
    samtools flagstat {output.out5} > {output.out5f}
    samtools idxstats {output.out5} > {output.out5i}
    """

rule insert_size:
    input:
        bam = lambda w : join(workpath,bam_dir,w.name + "." + w.ext + "." + extensions3[w.ext + "."])
    output:
        txt= join(workpath,qc_dir,"{name}.{ext}.insert_size_metrics.txt"),
        pdf= temp(join(workpath,qc_dir,"{name}.{ext}.insert_size_histogram.pdf")),
    params:
        rname="insert_size",
        picardver=config['tools']['PICARDVER'],
        javaram='16g',
    shell: """
    module load {params.picardver};
    java -Xmx{params.javaram} -jar ${{PICARDJARPATH}}/picard.jar CollectInsertSizeMetrics \\
        INPUT={input.bam} \\
        OUTPUT={output.txt} \\
        H={output.pdf}
    """
