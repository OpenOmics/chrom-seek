# Common paired-end rules:
#   - trim_pe
#   - kraken_pe
#   - BWA_PE
#   - insert_size

# trim, remove PolyX and remove blacklist reads
rule trim:
    """
    Data-processing step to remove adapter sequences and perform quality trimming
    prior to alignment the reference genome.  Adapters are composed of synthetic
    sequences and should be removed prior to alignment. Bwa mem aligns adapter-free 
    fastqs to ba lacklist that has anomalous, unstructured, or high signal in 
    next-generation sequencing experiments independent of cell line or experiment.
    Samtools view -f4 selects for reads unmapped and outpute a blacklist-sequences-free bam file.
    SamToFastq convers BAM file to FASTQs. All file processing is done in
    "tmp_dir": "/lscratch/$SLURM_JOBID/", meaning all files are lost at job completion
    except for final R1.trim.fastq.gz and R2.trim.fastq.gz

    @Input:
        Raw FastQ files
    @Output:
        Trimmed and blacklist-sequences-free FastQ files
    """
    input:
        file1=join(workpath,"{name}.R1.fastq.gz"),
        file2=provided(join(workpath,"{name}.R2.fastq.gz"),paired_end)
    output:
        outfq1=temp(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz")),
        outfq2=provided(temp(join(workpath,trim_dir,"{name}.R2.trim.fastq.gz")),paired_end)
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
        sample="{name}",
        tmpdir=tmpdir,
        paired_end = paired_end
    threads: 16
    shell: """
    module load {params.cutadaptver};
    module load {params.bwaver};
    module load {params.samtoolsver};
    module load {params.picardver};
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    if [ '{params.paired_end}' == True ];then
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
            -o ${{tmp}}/{params.sample}.R1.trim.fastq.gz \\
            -p ${{tmp}}/{params.sample}.R2.trim.fastq.gz \\
            {input.file1} {input.file2}
        
        if [ "{params.blacklistbwaindex}" != ""];
        then bwa mem -t {threads} \\
            {params.blacklistbwaindex} \\
            ${{tmp}}/{params.sample}.R1.trim.fastq.gz \\
            ${{tmp}}/{params.sample}.R2.trim.fastq.gz \\
            | samtools view -@{threads} \\
                -f4 \\
                -b \\
                -o ${{tmp}}/{params.sample}.bam;
        rm ${{tmp}}/{params.sample}.R1.trim.fastq.gz;
        rm ${{tmp}}/{params.sample}.R2.trim.fastq.gz;
        
        java -Xmx{params.javaram} -jar $PICARDJARPATH/picard.jar SamToFastq \\
            -VALIDATION_STRINGENCY SILENT \\
            -INPUT ${{tmp}}/{params.sample}.bam \\
            -FASTQ ${{tmp}}/{params.sample}.R1.trim.fastq \\
            -SECOND_END_FASTQ ${{tmp}}/{params.sample}.R2.trim.fastq \\
            -UNPAIRED_FASTQ ${{tmp}}/{params.sample}.unpaired.noBL.fastq
            
        rm ${{tmp}}/{params.sample}.bam;
        
        pigz -p {threads} ${{tmp}}/{params.sample}.R1.trim.fastq;
        pigz -p {threads} ${{tmp}}/{params.sample}.R2.trim.fastq;
        fi
        mv ${{tmp}}/{params.sample}.R1.trim.fastq.gz {output.outfq1};
        mv ${{tmp}}/{params.sample}.R2.trim.fastq.gz {output.outfq2};
    else
        cutadapt \\
            --nextseq-trim=2 \\
            --trim-n \\
            -n 5 \\
            -O 5 \\
            -q {params.leadingquality},{params.trailingquality} \\
            -m {params.minlen} \\
            -b file:{params.fastawithadaptersetd} \\
            -j {threads} \\
            -o ${{tmp}}/{params.sample}.R1.trim.fastq.gz \\
            {input.file1}
        
        if [ "{params.blacklistbwaindex}" != ""];
        then bwa mem -t {threads} \\
            {params.blacklistbwaindex} \\
            ${{tmp}}/{params.sample}.R1.trim.fastq.gz \\
            | samtools view -@{threads} \\
                -f4 \\
                -b \\
                -o ${{tmp}}/{params.sample}.bam;
        rm ${{tmp}}/{params.sample}.R1.trim.fastq.gz;
        
        java -Xmx{params.javaram} -jar $PICARDJARPATH/picard.jar SamToFastq \\
            -VALIDATION_STRINGENCY SILENT \\
            -INPUT ${{tmp}}/{params.sample}.bam \\
            -FASTQ ${{tmp}}/{params.sample}.R1.trim.fastq
            
        rm ${{tmp}}/{params.sample}.bam;
        
        pigz -p {threads} ${{tmp}}/{params.sample}.R1.trim.fastq;

        fi
        mv ${{tmp}}/{params.sample}.R1.trim.fastq.gz {output.outfq1};
    fi
        """

rule BWA:
    """
    Data processing rule to align trimmed and blacklist-sequences-free reads 
    to the reference genome using bwa mem aligner. Samtools sort sort alignments 
    by coordinates to produce bam file. Samtools flagstats summarizes reads by QC
    and number of alignments for each FLAG type. Samtools idxstats outputs
    stats by chromosomesequence,length, # mapped read-segments and 
    # unmapped read-segments. Resulting bam file is filtered by a mapQ value.

    @Input:
        Trimmed and blacklist-sequences-free FastQ files
    @Output:
        Bam file that have read aligned: sorted.bam
        Bam file that has reads aligned and filted by mapQ a value: Q5.bam
    """
    input:
        infq1 = join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        infq2 = provided(join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"), paired_end)
    params:
        d=join(workpath,bam_dir),
        rname='bwa',
        reference=config['references'][genome]['BWA'],
        bwaver=config['tools']['BWAVER'],
        samtoolsver=config['tools']['SAMTOOLSVER'],
        script=join(workpath,"workflow","scripts","bam_filter_by_mapq.py"),
        pythonver=config['tools']['PYTHONVER'],
        paired_end = paired_end
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
    if [ '{params.paired_end}' == True ];then
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
    else
        bwa mem -t {threads} {params.reference} {input.infq1} \\
            | samtools sort -@{threads} -o {output.outbam1}
        
        samtools index {output.outbam1}
        samtools flagstat {output.outbam1} > {output.flagstat1}
        samtools idxstats {output.outbam1} > {output.idxstat1}
        samtools view -b -q 6 {output.outbam1} -o {output.outbam2}
        
        samtools index {output.outbam2}
        samtools flagstat {output.outbam2} > {output.flagstat2}
        samtools idxstats {output.outbam2} > {output.idxstat2}
    fi
    """

rule picard_dedup:
    """
    Picard MarkDuplicates removes duplicates from bam file.
    Samtools flagstats summarizes reads by QC
    and number of alignments for each FLAG type. Samtools idxstats outputs
    stats by chromosomesequence,length, # mapped read-segments and 
    # unmapped read-segments. If assay is cfchip, deduplicated bam file is
    filtered for chromosomes 1 through 22 to produce: Q5DD.bam; 
    the original deduplicated file also convereted to bed file of fragments 
    via granges and rtracklayer: Q5DD.tagAlign.

    @Input:
        Bam file that has reads aligned and filted by mapQ a value: Q5.bam
    @Output:
        Deduplicated Q5DD.bam for all assays, plus Q5DD.tagAlign if cfchip assay

    """
    input: 
        bam2=join(workpath,bam_dir,"{name}.Q5.bam")
    output:
        out5=join(workpath,bam_dir,"{name}.Q5DD.bam"),
        out5f=join(workpath,bam_dir,"{name}.Q5DD.bam.flagstat"),
        out5i=join(workpath,bam_dir,"{name}.Q5DD.bam.idxstat"),
        out6=join(workpath,bam_dir,"{name}.bwa.Q5.duplic"),
    params:
        rname='dedup',
        picardver=config['tools']['PICARDVER'],
        samtoolsver=config['tools']['SAMTOOLSVER'],
        rver=config['tools']['RVER'],
        javaram='16g',
        tmpdir=tmpdir,
        tmpBam="{name}.Q5DD.withXY.bam",
        rscript=join(config['references'][genome]['cfChIP_TOOLS_SRC'], "bam2fragment.R"),
        out7 = lambda w: join(workpath,bam_dir, w.name+".Q5DD.tagAlign") \
            if assay=="cfchip" else ""
    shell: """
    module load {params.samtoolsver};
    module load {params.picardver};
    module load {params.rver}; 
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    
    if [ "{assay}" == "cfchip" ];then
      java -Xmx{params.javaram} \\
        -jar $PICARDJARPATH/picard.jar MarkDuplicates \\
        -I {input.bam2} \\
        -O {params.tmpBam} \\
        -TMP_DIR ${{tmp}} \\
        -VALIDATION_STRINGENCY SILENT \\
        -REMOVE_DUPLICATES true \\
        -METRICS_FILE {output.out6};
      samtools index {params.tmpBam};
      samtools view -b {params.tmpBam} chr{{1..22}} > {output.out5};
      Rscript {params.rscript} {params.tmpBam} {params.out7};
      rm {params.tmpBam} {params.tmpBam}.bai;
      samtools index {output.out5};
      samtools flagstat {output.out5} > {output.out5f};
      samtools idxstats {output.out5} > {output.out5i}; 
    else
      java -Xmx{params.javaram} \\
        -jar $PICARDJARPATH/picard.jar MarkDuplicates \\
        -I {input.bam2} \\
        -O {output.out5} \\
        -TMP_DIR ${{tmp}} \\
        -VALIDATION_STRINGENCY SILENT \\
        -REMOVE_DUPLICATES true \\
        -METRICS_FILE {output.out6};
      samtools index {output.out5};
      samtools flagstat {output.out5} > {output.out5f};
      samtools idxstats {output.out5} > {output.out5i}; 
    fi
    """

rule bam2bw:
    """
    bamCoverage converts bams to bigwig files for read visialization
    in UCSC Genome Browser or IGV.
    @Input:
        Files with extensions: sorted.bam and Q5DD.bam,
        for all samples (treatment and input controls)
    @Output:
       bigWig file with contains coordinates for an interval and 
       an associated score, RPGC
    """
    input:
        bam=join(workpath,bam_dir,"{name}.{ext}.bam"),
        # ppqt=join(workpath,bam_dir,"{ext}.ppqt.txt"),
    output:
        outbw=join(workpath,bw_dir,"{name}.{ext}.RPGC.bw"),
    params:
        rname="bam2bw",
        pe=pe,
        effectivegenomesize=config['references'][genome]['EFFECTIVEGENOMESIZE'],
    threads: int(allocated("threads", "bam2bw", cluster)),
    envmodules: config['tools']['DEEPTOOLSVER'],
    shell: """
    if [ "{params.pe}" == "yes" ]; then
        bamCoverage \\
            --bam {input.bam} \\
            -o {output.outbw} \\
            --binSize 25 \\
            --smoothLength 75 \\
            --numberOfProcessors {threads} \\
            --normalizeUsing RPGC \\
            --effectiveGenomeSize {params.effectivegenomesize} \\
            --centerReads;   
    fi
    """    

rule inputnorm:
    """
    bigwigCompare (deepTools) subracts input control from treatment bigWig file,
    normalizing treatment bigWig.
    @Input:
        Treatment sample, which has input control: extension Q5DD.RPGC.bw
        and input its control: extension Q5DD.RPGC.bw
    @Output:
       bigWig file of treatmment sample normalizes with its input control
    """
    input:
        chip = join(workpath,bw_dir,"{name}.Q5DD.RPGC.bw"),
        ctrl = lambda w : join(workpath,bw_dir,chip2input[w.name] + ".Q5DD.RPGC.bw")
    output:
        join(workpath,bw_dir,"{name}.Q5DD.RPGC.inputnorm.bw")
    params:
        rname="inputnorm",
    threads: int(allocated("threads", "inputnorm", cluster)),
    envmodules: config['tools']['DEEPTOOLSVER'],
    shell: """
    bigwigCompare \\
        --binSize 25 \\
        --outFileName {output} \\
        --outFileFormat 'bigwig' \\
        --bigwig1 {input.chip} \\
        --bigwig2 {input.ctrl} \\
        --operation 'subtract' \\
        --skipNonCoveredRegions \\
        -p {threads}
    """
