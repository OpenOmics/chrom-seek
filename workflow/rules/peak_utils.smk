
# Calling read coverage as peaks 
# ~~~~
# Genrally applicable rules
from os.path import join
from collections import defaultdict

# ~~ workflow configuration
workpath                        = config['project']['workpath']
genome                          = config['options']['genome']
paired_end                      = False if config['project']['nends'] == 1 else True
chip2input                      = config['project']['peaks']['inputs']
tmpdir                          = config['options']['tmp_dir']
assay                           = config['options']['assay']
blocking                        = False if set(blocks.values()) in ({None}, {''}) else True
block_add                       = "_block" if blocking else ""
homer_output_targets            = ['homerMotifs.all.motifs', 'motifFindingParameters.txt', 
                                'knownResults.txt', 'seq.autonorm.tsv', 'homerResults.html']

# Directory end points
bam_dir                         = join(workpath, "bam")
ppqt_dir                        = join(bam_dir, "ppqt")
genrich_dir                     = join(workpath, "Genrich")
macsN_dir                       = join(workpath, "macsNarrow")
macsB_dir                       = join(workpath, "macsBroad")
sicer_dir                       = join(workpath, "sicer")
homer_dir                       = join(workpath, "HOMER")


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


pkcaller2homer_size = defaultdict(lambda: "200")
pkcaller2homer_size.update({
    'macsNarow': "200", 
    'macsBroad': "given", 
    'genrich': "given", 
    'SEACR': "given"
})


rule HOMER:
     input:
         up_file                         = join(
                                             diffbind_dir,
                                             "{contrast}-{PeakTool}",
                                             "{contrast}-{PeakTool}_Diffbind" + block_add + "_{differential_app}_up.bed",
                                           ),
         down_file                       = join(
                                             diffbind_dir,
                                             "{contrast}-{PeakTool}",
                                             "{contrast}-{PeakTool}_Diffbind" + block_add + "_{differential_app}_down.bed",
                                           ),
     output:
         down_motifs                     = [join(homer_dir, 'DOWN', fn) for fn in homer_output_targets],
         up_motifs                       = [join(homer_dir, 'UP', fn) for fn in homer_output_targets],
     params:
         rname                           = 'HOMER',
         genome_name                     = genome,
         out_dir_up                      = join(homer_dir, 'UP'),
         out_dir_down                    = join(homer_dir, 'DOWN'),
         # -len <#>[,<#>,<#>...] (motif length, default=8,10,12) [NOTE: values greater 12 may cause the program
         #   to run out of memmory - in these cases decrease the number of sequences analyzed]
         seq_length                      = "8,10",
         # Selecting the size of the region for motif finding (-size # or -size given, default: 200)
         #   This is one of the most important parameters and also a source of confusion for many.  
         #   If you wish to find motifs using your peaks using their exact sizes, use the option "-size given").  
         # However, for Transcription Factor peaks, most of the motifs are found +/- 50-75 bp from the peak center, 
         #   making it better to use a fixed size rather than depend on your peak size.
         motif_finding_region            = pkcaller2homer_size["{PeakTool}"],
         tmpdir                          = tmpdir
     threads:
         int(cluster['HOMER'].get('threads', cluster['__default__']['threads']))
     shell:
         """
         if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
         tmp=$(mktemp -d -p "{params.tmpdir}")
         trap 'rm -rf "${{tmp}}"' EXIT
         export TMPDIR="${{tmp}}" # used by sort
         module load homer
         cd ${{tmp}}
         findMotifsGenome.pl {input.up} {params.genome_name} {params.out_dir_up} -p {threads} -size {params.motif_finding_region} -len {params.seq_len}
         findMotifsGenome.pl {input.down} {params.genome_name} {params.out_dir_down} -p {threads} -size {params.motif_finding_region} -len {params.seq_len}
         """
