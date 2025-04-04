
# Calling read coverage as peaks 
# ~~~~
# Genrally applicable rules
from os.path import join
from collections import defaultdict
from textwrap import dedent

# ~~ workflow configuration
workpath                        = config['project']['workpath']
genome                          = config['options']['genome']
homer_genome                    = config['references'][genome]['HOMER_REF']
genomefa                        = config['references'][genome]['GENOME']
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


pkcaller2homer_size = defaultdict(lambda: "given")
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
         down_motifs                     = temp([join(homer_dir, "DOWN_{contrast}_{PeakTool}_{differential_app}", fn) for fn in homer_output_targets]),
         down_motifs_gz                  = join(homer_dir, "DOWN_{contrast}_{PeakTool}_{differential_app}", "DOWN_{contrast}_{PeakTool}_{differential_app}.tar.gz"),
         up_motifs                       = temp([join(homer_dir, "UP_{contrast}_{PeakTool}_{differential_app}", fn) for fn in homer_output_targets]),
         up_motifs_gz                    = join(homer_dir, "UP_{contrast}_{PeakTool}_{differential_app}", "UP_{contrast}_{PeakTool}_{differential_app}.tar.gz"),
     params:
         rname                           = 'HOMER',
         homer_genome                    = homer_genome,
         genomealias                     = genome,
         genomefa                        = genomefa,
         out_dir_up                      = join(homer_dir, "UP_{contrast}_{PeakTool}_{differential_app}"),
         out_dir_down                    = join(homer_dir, "DOWN_{contrast}_{PeakTool}_{differential_app}"),
         seq_length                      = "8,10",
         motif_finding_region            = pkcaller2homer_size["{PeakTool}"],
         tmpdir                          = tmpdir,
         homer_peak_threshhold           = 20
     threads:
         int(cluster['HOMER'].get('threads', cluster['__default__']['threads']))
     shell:
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT
        export TMPDIR="${{tmp}}" # used by sort
        module load homer/4.11.1
        cd ${{TMPDIR}}
        [ -d "{params.homer_genome}" ] || {{ echo "Homer does not support this genome!" >&2; exit 1; }}
        for each in {params.homer_genome}/preparsed/*; do ln -s $each .; done
        ln -s {params.genomefa} ${{TMPDIR}}/{params.genomealias}
        uppeaks=$(wc -l {input.up_file} | cut -f1 -d$' ')
        downpeaks=$(wc -l {input.down_file} | cut -f1 -d$' ')
        thres=$(("{params.homer_peak_threshhold}"))
        thres=$((${{thres}} + 1))
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{print $4, $1, $2, $3, $6}}' {input.up_file} | sed -e 's/\./+/g' > ${{tmp}}/homer_up_input.bed
        if [ "${{uppeaks}}" -ge ${{thres}} ]; then
            echo "\\n\\n-------- HOMER UP_GENES_{wildcards.contrast}_{wildcards.PeakTool}_{wildcards.differential_app} sample sheet --------"
            cat ${{tmp}}/homer_up_input.bed
            findMotifsGenome.pl ${{tmp}}/homer_up_input.bed \\
                ${{TMPDIR}}/{params.genomealias} \\
                {params.out_dir_up} \\
                -preparsedDir ${{TMPDIR}} \\
                -p {threads} \\
                -mask \\
                -size {params.motif_finding_region} \\
                -len {params.seq_length}
            tar -czf {params.out_dir_up}/UP_{wildcards.contrast}_{wildcards.PeakTool}_{wildcards.differential_app}.tar.gz {params.out_dir_up}
            echo "-------- HOMER UP_GENES_{wildcards.contrast}_{wildcards.PeakTool}_{wildcards.differential_app} sample sheet --------\\n\\n"
        else
            touch {output.up_motifs}
            echo "{input.up_file} has less than 20 peaks; Not running homer!"
        fi
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{print $4, $1, $2, $3, $6}}' {input.down_file} | sed -e 's/\./+/g' > ${{tmp}}/homer_down_input.bed
        if [ "${{downpeaks}}" -ge ${{thres}} ]; then
            echo "\\n\\n-------- HOMER DOWN_GENES_{wildcards.contrast}_{wildcards.PeakTool}_{wildcards.differential_app} sample sheet --------"
            cat ${{tmp}}/homer_down_input.bed
            findMotifsGenome.pl ${{tmp}}/homer_down_input.bed \\
                ${{TMPDIR}}/{params.genomealias} \\
                {params.out_dir_down} \\
                -preparsedDir ${{TMPDIR}} \\
                -mask \\
                -p {threads} \\
                -size {params.motif_finding_region} \\
                -len {params.seq_length}
            tar -czf {params.out_dir_down}/UP_{wildcards.contrast}_{wildcards.PeakTool}_{wildcards.differential_app}.tar.gz {params.out_dir_down}
            echo "-------- HOMER DOWN_GENES_{wildcards.contrast}_{wildcards.PeakTool}_{wildcards.differential_app} sample sheet --------\\n\\n"
        else
            touch {output.down_motifs}
            echo "{input.down_file} has less than 20 peaks; Not running homer!"
        fi
        """)
