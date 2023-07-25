rule FRiP:
    input:
        bed = lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ],
        bam = join(workpath,bam_dir,"{Sample}.Q5DD.bam"),
    output:
        join(workpath,qc_dir,"{PeakTool}.{Sample}.Q5DD.FRiP_table.txt"),
    params:
        rname="frip",
        pythonver="python/3.5",
        outroot = lambda w: join(workpath,qc_dir,w.PeakTool),
        script=join(workpath,"workflow","scripts","frip.py"),
        genome = config['references']['REFLEN'],
        tmpdir = tmpdir,
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    module load {params.pythonver}
    export TMPDIR="${{tmp}}"  # pysam writes to $TMPDIR
    python {params.script} \\
        -p "{input.bed}" \\
        -b "{input.bam}" \\
        -g {params.genome} \\
        -o "{params.outroot}"
    """

rule FRiP_plot:
    input:
        expand(join(workpath,qc_dir,"{PeakTool}.{Sample}.Q5DD.FRiP_table.txt"), PeakTool=PeakTools, Sample=samples),
    output:
        expand(join(workpath, qc_dir, "{Group}.FRiP_barplot.pdf"),Group=groups),
    params:
        rname="frip_plot",
        Rver="R/4.2",
        script=join(workpath,"workflow","scripts","FRiP_plot.R"),
    shell: """
    module load {params.Rver}
    Rscript {params.script} \\
        {workpath}
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
        genome = config['references']['REFLEN']
    envmodules:
        config['tools']['BEDTOOLSVER']
    shell: """
    python {params.script} \\
        -i "{input}" \\
        -o "{params.outroot}" \\
        -g {params.genome}
    """