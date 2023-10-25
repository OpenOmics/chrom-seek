
# Helper functions
def get_input_bam(wildcards):
    """
    Returns a ChIP samples input BAM file,
    see chip2input for ChIP, Input pairs.
    """
    input_sample = chip2input[wildcards.name]
    if input_sample:
        # Runs in a ChIP, input mode
        return join(workpath, bam_dir, "{0}.Q5DD.bam".format(input_sample))
    else:
        # Runs in ChIP-only mode
        return []

rule sortByRead:
    """
    Bam files(extension: sorted.bam) need to be sorted by read name
    for Genrich.
    @Input:
        Bam file (extension: sorted.bam)
    @Output:
        Bam file sorted by read name (extension: sortByRead.bam)
    """
    input:
        join(workpath,bam_dir,"{name}.sorted.bam")
    output:
        temp(join(workpath,bam_dir,"{name}.sortedByRead.bam"))
    params:
        rname="sortByRead",
        samtools=config['tools']['SAMTOOLSVER'],
        mem=allocated("mem", "sortByRead", cluster)
    threads: int(allocated("threads", "sortByRead", cluster))
    shell: """
    module load {params.samtools}
    samtools sort {input} -n \\
        -@ {threads} \\
        -o {output}
    """

rule genrich:
    """
    In ATAC-mode, identifies open chromatic regions at intervals 
    centered on transposase cut sites (ends of the reads/fragments).
    @Input:
        Bam file sorted by read name (extension: sortByRead.bam)
    @Output:
        Output peak file (in ENCODE narrowPeak format):
        chrom name, star pos of peak, end pos of peak, peak name,
        AUC score, strand, area under AUC, summit -log(p-value),
        summit -log(q-value), and summit position.
    """
    input: 
        join(workpath,bam_dir,"{name}.sortedByRead.bam")
    output: 
        join(workpath,"Genrich","{name}","{name}.narrowPeak")
    params:
        rname="genrich",
        genrich_ver=config['tools']['GENRICHVER']
    shell: """
    module load {params.genrich_ver}
    Genrich \\
        -t {input} \\
        -o {output} \\
        -j \\
        -y \\
        -r \\
        -v \\
        -d 150 \\
        -m 5 \\
        -e chrM,chrY
    """

# INDIVIDUAL RULES
rule MACS2_narrow:
    input:
        chip = lambda w: join(workpath,bam_dir, w.name+".Q5DD.bam") \
        if paired_end else join(workpath,bam_dir, w.name+".Q5DD_tagAlign.gz"),
        txt = lambda w: join(workpath,bam_dir, ppqt_dir, w.name+".Q5DD.ppqt.txt") \
        if paired_end else join(workpath,bam_dir, ppqt_dir, w.name+".Q5DD_tagAlign.ppqt.txt"),
        c_option_pe = provided( lambda w: "{0}.Q5DD.bam".format(
            join(workpath, bam_dir, chip2input[w.name])
        ) if chip2input[w.name] else "", paired_end == True),
        c_option_se = provided(lambda w: "{0}.Q5DD_tagAlign.gz".format(
            join(workpath, bam_dir, chip2input[w.name])
        ) if chip2input[w.name] else "", paired_end == False) # Building optional argument for paired input option, input: '-c /path/to/input.Q5DD.bam', No input: ''
    output:
        join(workpath,macsN_dir,"{name}","{name}_peaks.narrowPeak"),
    params:
        rname='MACS2_narrow',
        gsize=config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver=config['tools']['MACSVER'],
        paired_end = paired_end,
        flag= lambda w: "-c" if chip2input[w.name] else ""
    shell: """
    module load {params.macsver};
    if [ '{params.paired_end}' == True ]; then
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option_pe} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {workpath}/{macsN_dir}/{wildcards.name} \\
            -q 0.01 \\
            --keep-dup="all" \\
            -f "BAMPE"
    else
        ppqt_len=$(awk '{{print $1}}' {input.txt})
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option_se} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {workpath}/{macsN_dir}/{wildcards.name} \\
            -q 0.01 \\
            --keep-dup="all" \\
            --nomodel \\
            --extsize $ppqt_len
    fi
    """

rule MACS2_broad:
    input:
        chip = join(workpath,bam_dir,"{name}.Q5DD.bam"),
        ctrl = get_input_bam,
    output:
        join(workpath,macsB_dir,"{name}","{name}_peaks.broadPeak"),
    params:
        rname='MACS2_broad',
        gsize=config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver=config['tools']['MACSVER'],
        # Building optional argument for paired input option,
        # input: '-c /path/to/input.Q5DD.bam', No input: ''
        c_option = lambda w: "-c {0}.Q5DD.bam".format(
            join(workpath, bam_dir, chip2input[w.name])
        ) if chip2input[w.name] else "",
    shell: """
    module load {params.macsver};
    macs2 callpeak \\
        -t {input.chip} {params.c_option} \\
        -g {params.gsize} \\
        -n {wildcards.name} \\
        --outdir {workpath}/{macsB_dir}/{wildcards.name} \\
        --broad \\
        --broad-cutoff 0.01 \\
        --keep-dup="all" \\
        -f "BAMPE"
    """

rule SICER:
    input: 
        chip = join(workpath,bam_dir,"{name}.Q5DD.bam"),
        ctrl = lambda w : join(workpath,bam_dir,chip2input[w.name] + ".Q5DD.bam"),
        fragLen= join(workpath,"QC","{name}.Q5DD.insert_size_metrics.txt"),
    output:
        bed = join(workpath,sicer_dir,"{name}","{name}_broadpeaks.bed"),
    params:
        rname='SICER',
        sicerver=config['tools']['SICERVER'],
        bedtoolsver=config['tools']['BEDTOOLSVER'],
        genomever = config['options']['genome'],
        name="{name}",
        sicer_dir=join(workpath,sicer_dir,"{name}"),
        tmpdir=tmpdir
    shell: """
    module load {params.sicerver}
    module load {params.bedtoolsver}
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    MEAN_INSERT_SIZE=$(cat {input.fragLen} | awk '/MEDIAN_INSERT_SIZE/{{f=1;next}} /## HISTOGRAM/{{f=0}} f' | cut -f 6)
    mean_insert_size=$(printf "%.0f" $MEAN_INSERT_SIZE)

    if [ -f "{input.ctrl}" ]; 
    then
        sicer \\
        -t {input.chip} \\
        -c {input.ctrl} \\
        -s {params.genomever} \\
        -rt 100 \\
        -w 300 \\
        -f ${{mean_insert_size}} \\
        -egf 0.75 \\
        -g 600 \\
        -fdr 1E-2 \\
        -cpu 30 \\
        -o .
        mv {params.name}.Q5DD-W300-G600-FDR0.01-island.bed {output.bed};
        mv {params.name}.Q5DD-W300-G600-islands-summary {params.sicer_dir}
    else
        sicer \\
        -t {input.chip} \\
        -s {params.genomever} \\
        -rt 100 \\
        -w 300 \\
        -f ${{mean_insert_size}} \\
        -egf 0.75 \\
        -g 600 \\
        -e 100 \\
        -cpu 30 \\
        -o .
        mv {params.name}.Q5DD-W300-G600.scoreisland {output.bed}
    fi
    """