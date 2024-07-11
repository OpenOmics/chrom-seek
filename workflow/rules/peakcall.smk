
# Quality control rules
# ~~~~
# Common quality-control rules: preseq, NRF, rawfastqc,
#   fastqc, fastq_screen, multiQC
from os.path import join
from scripts.peakcall import get_control_input, getMacTXT, getMacChip, \
    getSicerChips, getSicerFragLen, get_control_input


# ~~ workflow configuration
workpath                        = config['project']['workpath']
genome                          = config['options']['genome']
paired_end                      = False if config['project']['nends'] == 1 else True
chip2input                      = config['project']['peaks']['inputs']

# Directory end points
bam_dir                         = join(workpath, "bam")
ppqt_dir                        = join(bam_dir, "ppqt")
genrich_dir                     = join(workpath, "Genrich")
macsN_dir                       = join(workpath, "macsNarrow")
macsB_dir                       = join(workpath, "macsBroad")
sicer_dir                       = join(workpath, "sicer")


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
        join(bam_dir, "{name}.sortedByRead.bam")
    output: 
        join(genrich_dir, "{name}", "{name}.narrowPeak")
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
        chip = lambda w: getMacChip(bam_dir, w.name, paired_end),
        txt = lambda w: getMacTXT(ppqt_dir, w.name, paired_end),
        c_option = lambda w: get_control_input(chip2input[w.name], paired_end, bam_dir),
    output:
        join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"),
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
            -t {input.chip} {params.flag} {input.c_option} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {macsN_dir}/{wildcards.name} \\
            -q 0.01 \\
            --keep-dup="all" \\
            -f "BAMPE"
    else
        ppqt_len=$(awk '{{print $1}}' {input.txt})
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {macsN_dir}/{wildcards.name} \\
            -q 0.01 \\
            --keep-dup="all" \\
            --nomodel \\
            --extsize $ppqt_len
    fi
    """

rule MACS2_broad:
    input:
        chip = lambda w: getMacChip(bam_dir, w.name, paired_end),
        txt = lambda w: getMacTXT(ppqt_dir, w.name, paired_end),
        c_option = lambda w: get_control_input(chip2input[w.name], paired_end, bam_dir),
    output:
        join(macsB_dir, "{name}", "{name}_peaks.broadPeak"),
    params:
        rname='MACS2_broad',
        gsize=config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver=config['tools']['MACSVER'],
        paired_end = paired_end,
        flag= lambda w: "-c" if chip2input[w.name] else ""
    shell: """
    module load {params.macsver};
    if [ '{params.paired_end}' == True ]; then
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {macsB_dir}/{wildcards.name} \\
            --broad \\
            --broad-cutoff 0.01 \\
            --keep-dup="all" \\
            -f "BAMPE"
    else 
        ppqt_len=$(awk '{{print $1}}' {input.txt})
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {macsB_dir}/{wildcards.name} \\
            --broad \\
            --broad-cutoff 0.01 \\
            --keep-dup="all" \\
            --nomodel \\
            --extsize $ppqt_len
    fi
    """

rule SICER:
    input: 
        chip = lambda w: getSicerChips(bam_dir, w.name, paired_end),
        fragLen = lambda w: getSicerFragLen(ppqt_dir, qc_dir, w.name, paired_end),
        c_option = lambda w: get_control_input(chip2input[w.name], paired_end, bam_dir),
    output:
        bed = join(sicer_dir, "{name}", "{name}_broadpeaks.bed"),
    params:
        rname='SICER',
        sicerver=config['tools']['SICERVER'],
        bedtoolsver=config['tools']['BEDTOOLSVER'],
        genomever = config['options']['genome'],
        name="{name}",
        sicer_dir=join(sicer_dir,"{name}"),
        tmpdir=tmpdir,
        paired_end = paired_end,
        frac=config['references'][genome]['FRAC'],
        flag= lambda w: "-c" if chip2input[w.name] else "",
    shell: """
    module load {params.sicerver}
    module load {params.bedtoolsver}
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    if [ '{params.paired_end}' == True ]; then
        MEAN_INSERT_SIZE=$(cat {input.fragLen} | awk '/MEDIAN_INSERT_SIZE/{{f=1;next}} /## HISTOGRAM/{{f=0}} f' | cut -f 6)
        mean_insert_size=$(printf "%.0f" $MEAN_INSERT_SIZE)
    else
        mean_insert_size=$(awk '{{print $1}}' {input.fragLen})
    fi
    echo "printing out value of mean-insert-size ${{mean_insert_size}}"
    a={input.c_option}
    echo "Printing input.c_option ${{a}}"
    if [ '{params.paired_end}' == True ]; then
        if [ -f "{input.c_option}" ]; then
            # Copying input to tmpdir due to SICER2
            # bam2bed file conversion, if more than
            # one sample shares the same IP sample
            # than a race condition can occur where 
            # two jobs can concurrent try to write 
            # to the same BED file (during bedtools
            # bam2bed that sicer calls).
            input_bam="$(basename "{input.c_option}")"
            cp {input.c_option} ${{tmp}}
            echo "paired-end with input... ${{tmp}}/${{input_bam}}"
            sicer \\
            -t {input.chip} \\
            -c "${{tmp}}/${{input_bam}}" \\
            -s {params.genomever} \\
            -rt 100 \\
            -w 300 \\
            -f ${{mean_insert_size}} \\
            -egf {params.frac} \\
            -g 600 \\
            -fdr 1E-2 \\
            -cpu 30 \\
            -o ${{tmp}}
            
            mv ${{tmp}}/{params.name}.Q5DD-W300-G600-FDR0.01-island.bed {output.bed};
            mv ${{tmp}}/{params.name}.Q5DD-W300-G600-islands-summary {params.sicer_dir}
        else
            echo "paired-end without input"
            sicer \\
            -t {input.chip} \\
            -s {params.genomever} \\
            -rt 100 \\
            -w 300 \\
            -f ${{mean_insert_size}} \\
            -egf {params.frac} \\
            -g 600 \\
            -e 100 \\
            -cpu 30 \\
            -o ${{tmp}}

            mv ${{tmp}}/{params.name}.Q5DD-W300-G600.scoreisland {params.sicer_dir}
        fi
    else
        if [ -f "{input.c_option}" ]; then
            echo "single-end with input"
            cp {input.chip} ${{tmp}}/chip.bed.gz; gzip -d ${{tmp}}/chip.bed.gz;
            awk 'BEGIN{{FS=OFS="\\t"}} {{gsub(/\./, 0, $5)}} 1' ${{tmp}}/chip.bed > ${{tmp}}/{params.name}.bed;

            cp {input.c_option} ${{tmp}}/input.bed.gz; gzip -d ${{tmp}}/input.bed.gz;
            awk 'BEGIN{{FS=OFS="\\t"}} {{gsub(/\./, 0, $5)}} 1' ${{tmp}}/input.bed > ${{tmp}}/inputV2.bed;
            
            sicer \\
            -t ${{tmp}}/{params.name}.bed \\
            -c ${{tmp}}/inputV2.bed \\
            -s {params.genomever} \\
            -rt 100 \\
            -w 300 \\
            -f ${{mean_insert_size}} \\
            -egf {params.frac} \\
            -g 600 \\
            -fdr 1E-2 \\
            -cpu 30 \\
            -o ${{tmp}}
            mv ${{tmp}}/{params.name}-W300-G600-FDR0.01-island.bed {output.bed};
            mv ${{tmp}}/{params.name}-W300-G600-islands-summary {params.sicer_dir}
        else
            echo "single-end without input"
            cp {input.chip} ${{tmp}}/chip.bed.gz; gzip -d ${{tmp}}/chip.bed.gz;
            awk 'BEGIN{{FS=OFS="\\t"}} {{gsub(/\./, 0, $5)}} 1' ${{tmp}}/chip.bed > ${{tmp}}/{params.name}.bed;
            sicer \\
            -t ${{tmp}}/{params.name}.bed \\
            -s {params.genomever} \\
            -rt 100 \\
            -w 300 \\
            -f ${{mean_insert_size}} \\
            -egf {params.frac} \\
            -g 600 \\
            -e 100 \\
            -cpu 30 \\
            -o ${{tmp}}
            mv ${{tmp}}/{params.name}-W300-G600.scoreisland {output.bed}
        fi
    fi
    """

rule MEME:
    input:
        bed = lambda w: join(workpath, w.PeakTool, w.name, w.name + PeakExtensions[w.PeakTool])
    output:
        meme_out = join(MEME_dir, "{PeakTool}", "{name}_meme", "meme-chip.html"),
        ame_out = join(MEME_dir, "{PeakTool}", "{name}_ame", "ame.html")
    params:
        rname='MEME',
        ref_fa=config['references'][genome]['GENOME'],
        meme_vertebrates_db=config['references'][genome]['MEME_VERTEBRATES_DB'],
        meme_euk_db=config['references'][genome]['MEME_EUKARYOTE_DB'],
        meme_genome_db=config['references'][genome]['MEME_GENOME_DB'],
        oc=join(MEME_dir, "{PeakTool}", "{name}"),
        tmpdir=tmpdir,
        outfa="{name}.fa",
        ntasks=int(28)
    shell: """
    module load meme
    module load bedtools
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
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