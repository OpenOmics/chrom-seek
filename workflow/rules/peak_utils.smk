
# Calling read coverage as peaks 
# ~~~~
# Genrally applicable rules
from os.path import join


# ~~ workflow configuration
workpath                        = config['project']['workpath']
genome                          = config['options']['genome']
paired_end                      = False if config['project']['nends'] == 1 else True
chip2input                      = config['project']['peaks']['inputs']
tmpdir                          = config['options']['tmp_dir']
assay                           = config['options']['assay']

# Directory end points
bam_dir                         = join(workpath, "bam")
ppqt_dir                        = join(bam_dir, "ppqt")
genrich_dir                     = join(workpath, "Genrich")
macsN_dir                       = join(workpath, "macsNarrow")
macsB_dir                       = join(workpath, "macsBroad")
sicer_dir                       = join(workpath, "sicer")


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
        chip                                = lambda wc: join(bw_dir, f"{wc.name}.Q5DD.RPGC.bw"),
        ctrl                                = lambda wc: join(bw_dir, f"{chip2input[wc.name]}.Q5DD.RPGC.bw") 
                                                if wc.name in chip2input and chip2input[wc.name] else []
    output:
        join(bw_dir, "{name}.Q5DD.RPGC.inputnorm.bw")
    params:
        rname                               = "inputnorm",
        bigwig_declare                      = lambda wc, input: f"--bigwig1 {input.chip} --bigwig2 {input.ctrl}",
    threads: 
        int(allocated("threads", "inputnorm", cluster)),
    envmodules: 
        config['tools']['DEEPTOOLSVER'],
    shell: 
        """
        bigwigCompare \\
            --binSize 25 \\
            --outFileName {output} \\
            --outFileFormat 'bigwig' \\
            {params.bigwig_declare} \\
            --operation 'subtract' \\
            --skipNonCoveredRegions \\
            -p {threads}
        """


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
        join(bam_dir, "{name}.sorted.bam")
    output:
        temp(join(bam_dir, "{name}.sortedByRead.bam"))
    params:
        rname                           = "sortByRead",
        samtools                        = config['tools']['SAMTOOLSVER'],
        mem                             = allocated("mem", "sortByRead", cluster)
    threads: 
        int(allocated("threads", "sortByRead", cluster))
    shell: 
        """
        module load {params.samtools}
        samtools sort {input} -n \\
            -@ {threads} \\
            -o {output}
        """


rule MEME:
    input:
        bed                             = lambda w: join(workpath, w.PeakTool, w.name, w.name + PeakExtensions[w.PeakTool])
    output:
        meme_out                        = join(MEME_dir, "{PeakTool}", "{name}_meme", "meme-chip.html"),
        ame_out                         = join(MEME_dir, "{PeakTool}", "{name}_ame", "ame.html")
    params:
        rname                           = 'MEME',
        ref_fa                          = config['references'][genome]['GENOME'],
        meme_vertebrates_db             = config['references'][genome]['MEME_VERTEBRATES_DB'],
        meme_euk_db                     = config['references'][genome]['MEME_EUKARYOTE_DB'],
        meme_genome_db                  = config['references'][genome]['MEME_GENOME_DB'],
        oc                              = join(MEME_dir, "{PeakTool}", "{name}"),
        outfa                           = "{name}.fa",
        ntasks                          = int(28)
    shell: 
        """
        module load meme
        module load bedtools
        if [ ! -d \"""" + str(tmpdir) + """\" ]; then mkdir -p \"""" + str(tmpdir) + """\"; fi
        tmp=$(mktemp -d -p \"""" + str(tmpdir) + """\")
        trap 'rm -rf "${{tmp}}"' EXIT

        bedtools getfasta -fi {params.ref_fa} -bed {input.bed} -fo ${{tmp}}/{params.outfa}
        meme-chip \\
        --oc {params.oc}_meme \\
        -db {params.meme_vertebrates_db} \\
        -meme-searchsize 34000000 \\
        -meme-p {params.ntasks} \\
        ${{tmp}}/{params.outfa}

        ame \\
        --oc {params.oc}_ame ${{tmp}}/{params.outfa} \\
        {params.meme_euk_db} {params.meme_vertebrates_db} {params.meme_genome_db}
        """