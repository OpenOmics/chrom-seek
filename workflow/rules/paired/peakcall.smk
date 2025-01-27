# Calling read coverage as peaks 
# ~~~~
# Paired-end peak calling rules

# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
tmpdir                          = config['options']['tmp_dir']


rule MACS2_broad:
    input:
        chip                            = join(bam_dir, "{name}.Q5DD.bam"),
        c_option                        = lambda w: join(bam_dir, f"{chip2input[w.name]}.Q5DD.bam") 
                                            if chip2input[w.name] else [],
    output:
        join(macsB_dir, "{name}", "{name}_peaks.broadPeak"),
    params:
        rname                           = 'MACS2_broad',
        gsize                           = config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver                         = config['tools']['MACSVER'],
        flag                            = lambda w: "-c" if chip2input[w.name] else "",
    shell: 
        """
        module load {params.macsver};
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {macsB_dir}/{wildcards.name} \\
            --broad \\
            --broad-cutoff 0.01 \\
            --keep-dup="all" \\
            -f "BAMPE"
        """


rule MACS2_narrow:
    input:
        chip                            = join(bam_dir, "{name}.Q5DD.bam"),
        c_option                        = lambda w: join(bam_dir, f"{chip2input[w.name]}.Q5DD.bam")
                                            if chip2input[w.name] else [],
    output:
        join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"),
    params:
        rname                           = 'MACS2_narrow',
        gsize                           = config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver                         = config['tools']['MACSVER'],
        flag                            = lambda w: "-c" if chip2input[w.name] else "",
    shell: 
        """
        module load {params.macsver};
        macs2 callpeak \\
                -t {input.chip} {params.flag} {input.c_option} \\
                -g {params.gsize} \\
                -n {wildcards.name} \\
                --outdir {macsN_dir}/{wildcards.name} \\
                -q 0.01 \\
                --keep-dup="all" \\
                -f "BAMPE"
        """


rule SICER:
    input:
        chip                            = lambda w: join(bam_dir, w.name + ".Q5DD.bam"),
        fragLen                         = lambda w: join(qc_dir, name + ".Q5DD.insert_size_metrics.txt"),
        c_option                        = lambda w: join(bam_dir, f"{chip2input[w.name]}.Q5DD.bam")
                                                        if chip2input[w.name] else [],
    output:
        bed                             = join(sicer_dir, "{name}", "{name}_broadpeaks.bed") if has_inputs else [],
    params:
        rname                           = 'SICER',
        name                            = "{name}",
        sicerver                        = config['tools']['SICERVER'],
        bedtoolsver                     = config['tools']['BEDTOOLSVER'],
        genomever                       = config['options']['genome'],
        this_sicer_dir                  = join(sicer_dir,"{name}"),
        frac                            = config['references'][genome]['FRAC'],
        flag                            = lambda w: "-c" if chip2input[w.name] else "",
    shell: 
        """
        module load {params.sicerver}
        module load {params.bedtoolsver}
        if [ ! -d \"""" + str(tmpdir) + """\" ]; then mkdir -p \"""" + str(tmpdir) + """\"; fi
        tmp=$(mktemp -d -p \"""" + str(tmpdir) + """\")
        trap 'rm -rf "${{tmp}}"' EXIT

        MEAN_INSERT_SIZE=$(cat {input.fragLen} | awk '/MEDIAN_INSERT_SIZE/{{f=1;next}} /## HISTOGRAM/{{f=0}} f' | cut -f 6)
        mean_insert_size=$(printf "%.0f" $MEAN_INSERT_SIZE)

        echo "printing out value of mean-insert-size ${{mean_insert_size}}"
        a={input.c_option}
        echo "Printing input.c_option ${{a}}"
        
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
            mv ${{tmp}}/{params.name}.Q5DD-W300-G600-islands-summary {params.this_sicer_dir}
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

            mv ${{tmp}}/{params.name}.Q5DD-W300-G600.scoreisland {params.this_sicer_dir}
        fi
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
        rname                           = "genrich",
        genrich_ver                     = config['tools']['GENRICHVER']
    shell: 
        """
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