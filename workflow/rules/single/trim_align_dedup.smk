# Trimming, alignment, and redundancy reduction rules
# ~~~~
# Rules for single-ended reads

# ~~ workflow endpoints
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']

# ~~ directories
trim_dir                        = join(workpath, 'trim')
tmpdir                          = config['options']['tmp_dir']
bam_dir                         = join(workpath, "bam")
bw_dir                          = join(workpath, "bigwig")
ppqt_dir                        = join(bam_dir, "ppqt")


wildcard_constraints:
    # - regex to avoid backslashes in name wild cards
    #   this corrects routing for `ppqt` and 
    #   `ppqt_tagalign` rules
    # - if this isn't included, ppqt file paths (because
    #   they're nested in the bam directory) get mixed up
    #   with regular bam paths
    # - also possible to use negative look around regex:
    #   "^((?!\/).)*$"
    name                                = "[A-Za-z0-9_-]+"


rule trim:
    """
    Data-processing step to remove adapter sequences and perform quality trimming
    prior to alignment the reference genome.  Adapters are composed of synthetic
    sequences and should be removed prior to alignment. Bwa mem aligns adapter-free 
    fastqs to ba lacklist that has anomalous, unstructured, or high signal in 
    next-generation sequencing experiments independent of cell line or experiment.
    Samtools view -f4 selects for reads unmapped and outpute a blacklist-sequences-free bam file.
    SamToFastq convers BAM file to FASTQs. All file processing is done in
    "tmp_dir": "/lscratch/$SLURM_JOB_ID/", meaning all files are lost at job completion
    except for final R1.trim.fastq.gz and R2.trim.fastq.gz

    @Input:
        Raw FastQ files
    @Output:
        Trimmed and blacklist-sequences-free FastQ files
    """
    input:
        join(workpath, "{name}.R1.fastq.gz"),
    output:
        temp(join(trim_dir, "{name}.R1.trim.fastq.gz")),
    params:
        rname                               = "trim",
        cutadaptver                         = config['tools']['CUTADAPTVER'],
        workpath                            = config['project']['workpath'],
        fastawithadaptersetd                = join(workpath, config['shared_resources']['ADAPTERS_FASTA']),
        blacklistbwaindex                   = config['references'][genome]['BLACKLISTBWAINDEX'],
        picardver                           = config['tools']['PICARDVER'],
        bwaver                              = config['tools']['BWAVER'],
        samtoolsver                         = config['tools']['SAMTOOLSVER'],
        minlen                              = 35,
        leadingquality                      = 10,
        trailingquality                     = 10,
        javaram                             = "64g",
        sample                              = "{name}",
        tmpdir                              = tmpdir,
    threads: 
        32
    shell: 
        """
        module load {params.cutadaptver};
        module load {params.bwaver};
        module load {params.samtoolsver};
        module load {params.picardver};
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

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
            {input}
        
        if [ "{params.blacklistbwaindex}" != "" ]; then 
            thr=$(echo {threads}/2 | bc)
            bwa mem -t ${{thr}} \\
                {params.blacklistbwaindex} \\
                ${{tmp}}/{params.sample}.R1.trim.fastq.gz \\
                | samtools view -@${{thr}} \\
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
        mv ${{tmp}}/{params.sample}.R1.trim.fastq.gz {output};
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
        join(trim_dir, "{name}.R1.trim.fastq.gz"),
    params:
        d                                   = join(bam_dir),
        rname                               = 'bwa',
        reference                           = config['references'][genome]['BWA'],
        bwaver                              = config['tools']['BWAVER'],
        samtoolsver                         = config['tools']['SAMTOOLSVER'],
        script                              = join(workpath, "bin", "bam_filter_by_mapq.py"),
        pythonver                           = config['tools']['PYTHONVER'],
    output:
        outbam1                             = join(bam_dir, "{name}.sorted.bam"), 
        outbam2                             = temp(join(bam_dir, "{name}.Q5.bam")),
        flagstat1                           = join(bam_dir, "{name}.sorted.bam.flagstat"),
        idxstat1                            = join(bam_dir, "{name}.sorted.bam.idxstat"),
        flagstat2                           = join(bam_dir, "{name}.Q5.bam.flagstat"),
        idxstat2                            = join(bam_dir, "{name}.Q5.bam.idxstat"),
    threads: 
        32
    shell: 
        """
        module load {params.bwaver};
        module load {params.samtoolsver};
        module load {params.pythonver};
        
        bwa mem -t {threads} {params.reference} {input} \\
            | samtools sort -@{threads} -o {output.outbam1}
        samtools index {output.outbam1}
        samtools flagstat {output.outbam1} > {output.flagstat1}
        samtools idxstats {output.outbam1} > {output.idxstat1}
        samtools view -b -q 6 {output.outbam1} -o {output.outbam2}
        
        samtools index {output.outbam2}
        samtools flagstat {output.outbam2} > {output.flagstat2}
        samtools idxstats {output.outbam2} > {output.idxstat2}
        """


rule dedup:
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
        join(bam_dir, "{name}.Q5.bam")
    output:
        out_bam                             = join(bam_dir, "{name}.Q5DD.bam"),
        flagstat                            = join(bam_dir, "{name}.Q5DD.bam.flagstat"),
        idxstat                             = join(bam_dir, "{name}.Q5DD.bam.idxstat"),
        tagalign                            = join(bam_dir, "{name}.Q5DD_tagAlign.gz")
    params:
        rname                               = 'dedup',
        samtoolsver                         = config['tools']['SAMTOOLSVER'],
        bedtoolsver                         = config['tools']['BEDTOOLSVER'],
        macsver                             = config['tools']['MACSVER'],
        gsize                               = config['references'][genome]['EFFECTIVEGENOMESIZE'],
        genomefile                          = config['references'][genome]['REFLEN'],
        tmpdir                              = tmpdir,
    shell: 
        """
        module load {params.samtoolsver};
        module load {params.bedtoolsver};
        module load {params.macsver};
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT
        
        macs2 filterdup -i {input} -g {params.gsize} --keep-dup="auto" -o ${{tmp}}/TmpTagAlign;
        awk -F"\\t" -v OFS="\\t" '{{if ($2>0 && $3>0) {{print}}}}' ${{tmp}}/TmpTagAlign > ${{tmp}}/TmpTagAlign2;
        awk -F"\\t" -v OFS="\\t" '{{print $1,1,$2}}' {params.genomefile} | sort -k1,1 -k2,2n > ${{tmp}}/GenomeFileBed;
        bedtools intersect -wa -f 1.0 -a ${{tmp}}/TmpTagAlign2 -b ${{tmp}}/GenomeFileBed > ${{tmp}}/TmpTagAlign3;
        bedtools bedtobam -i ${{tmp}}/TmpTagAlign3 -g {params.genomefile} | samtools sort -@4 -o {output.out_bam};
        gzip ${{tmp}}/TmpTagAlign3;
        mv ${{tmp}}/TmpTagAlign3.gz {output.tagalign};
        samtools index {output.out_bam};
        samtools flagstat {output.out_bam} > {output.flagstat}
        samtools idxstats {output.out_bam} > {output.idxstat}
        """


rule ppqt:
    input:
        bam                                 = join(bam_dir, "{name}.{ext}.bam")
    output:                                          
        ppqt                                = join(ppqt_dir, "{name}.{ext}.ppqt.txt"),
        frag_length                         = join(ppqt_dir, "{name}.{ext}.fragment.length"),
        pdf                                 = join(ppqt_dir, "{name}.{ext}.pdf"),
    params:
        rname                               = "ppqt",
        samtoolsver                         = config['tools']['SAMTOOLSVER'],
        rver                                = config['tools']['RVER'],
        tmpdir                              = tmpdir,
        ppqt_script                         = join(bin_path, 'ppqt_process.py')
    container: 
        config['images']['ppqt']
    threads:
        int(cluster['ppqt'].get('threads', cluster['__default__']['threads']))
    shell: 
        """
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        run_spp.R \\
            -c={input.bam} \\
            -savp={output.pdf} \\
            -out={output.ppqt} \\
            -tmpdir=${{tmp}} \\
            -p={threads}
        python {params.ppqt_script} {output.ppqt} > {output.frag_length}
        """


rule ppqt_tagalign:
    input:
        tagalign                            = join(bam_dir, "{name}.Q5DD_tagAlign.gz")
    output:                                          
        ppqt                                = join(ppqt_dir, "{name}.Q5DD_tagAlign.ppqt.txt"),
        frag_length                         = join(ppqt_dir, "{name}.Q5DD_tagAlign.fragment.length"),
        pdf                                 = join(ppqt_dir, "{name}.Q5DD_tagAlign.pdf"),
    params:
        rname                               = "ppqt_tagalign",
        samtoolsver                         = config['tools']['SAMTOOLSVER'],
        rver                                = config['tools']['RVER'],
        tmpdir                              = tmpdir,
        ppqt_script                         = join(bin_path, 'ppqt_process.py')
    container: 
        config['images']['ppqt']
    threads:
        int(cluster['ppqt_tagalign'].get('threads', cluster['__default__']['threads']))
    shell: 
        """
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT
        # ppqt will not work unless file name ends with ".tagAlign"
        #   "_tagAlign" will not work
        gunzip -c {input.tagalign} > ${{tmp}}/sample.tagAlign
        run_spp.R \\
            -c=${{tmp}}/sample.tagAlign \\
            -savp={output.pdf} \\
            -out={output.ppqt} \\
            -p={threads} \\
            -tmpdir=${{tmp}} \\
            -rf
        python {params.ppqt_script} {output.ppqt} > {output.frag_length}
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
        bam                                 = join(bam_dir, "{name}.{ext}.bam"),
        ppqt                                = join(ppqt_dir, "{name}.{ext}.ppqt.txt"),
        frag_len                            = join(ppqt_dir, "{name}.{ext}.fragment.length"),
    output:
        outbw                               = join(bw_dir, "{name}.{ext}.RPGC.bw"),
    params:
        rname                               = "bam2bw",
        name                                = "{name}",
        effectivegenomesize                 = config['references'][genome]['EFFECTIVEGENOMESIZE'],
        tmpdir                              = tmpdir,
        frag_len_script                     = join(bin_path, "ppqt_process.py"),
    threads: 
        int(allocated("threads", "bam2bw", cluster)),
    envmodules: 
        config['tools']['DEEPTOOLSVER'],
    shell: 
        """
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        ppqt_len=$(cat {input.frag_len})
        
        bamCoverage \\
            --bam {input.bam} \\
            -o {output.outbw} \\
            --binSize 25 \\
            --smoothLength 75 \\
            --numberOfProcessors {threads} \\
            --normalizeUsing RPGC \\
            --effectiveGenomeSize {params.effectivegenomesize} \\
            -e ${{ppqt_len}}
        """