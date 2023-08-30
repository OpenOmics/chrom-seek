rule FRiP:
    input:
        bed = lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ],
        bam = join(workpath,bam_dir,"{name}.Q5DD.bam"),
    output:
        join(workpath,"PeakQC","{PeakTool}.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname="frip",
        outroot = lambda w: join(workpath,"PeakQC",w.PeakTool),
        script=join(workpath,"workflow","scripts","frip.py"),
        genome = config['references'][genome]['REFLEN'],
        tmpdir = tmpdir,
    container: config['images']['python']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    python {params.script} \\
        -p {input.bed} \\
        -b {input.bam} \\
        -g {params.genome} \\
        -o {params.outroot}
    """

rule jaccard:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ],
    output:
        join(workpath,qc_dir,'{PeakTool}_jaccard.txt'),
    params:
        rname="jaccard",
        outroot = lambda w: join(workpath,qc_dir,w.PeakTool),
        script=join(workpath,"workflow","scripts","jaccard_score.py"),
        genome = config['references'][genome]['REFLEN']
    envmodules:
        config['tools']['BEDTOOLSVER']
    shell: """
    python {params.script} \\
        -i "{input}" \\
        -o "{params.outroot}" \\
        -g {params.genome}
    """
