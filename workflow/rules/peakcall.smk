# INDIVIDUAL RULES
rule MACS2_narrow:
    input:
        chip = join(workpath,bam_dir,"{name}.Q5DD.bam"),
        ctrl = get_input_bam,
    output:
        join(workpath,macsN_dir,"{name}","{name}_peaks.narrowPeak"),
    params:
        rname='MACS2_narrow',
        gsize=config['references']['EFFECTIVEGENOMESIZE'],
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