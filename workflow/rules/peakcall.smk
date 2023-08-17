
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

if assay == "atac": 
    rule genrich_all:
        input:
            expand(join(workpath,"Genrich","{name}","{name}.narrowPeak"),name=chips)

rule sortByRead:
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
        chip = join(workpath,bam_dir,"{name}.Q5DD.bam"),
        ctrl = get_input_bam,
    output:
        join(workpath,macsN_dir,"{name}","{name}_peaks.narrowPeak"),
    params:
        rname='MACS2_narrow',
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
        --outdir {workpath}/{macsN_dir}/{wildcards.name} \\
        -q 0.01 \\
        --keep-dup="all" \\
        -f "BAMPE"
    """
