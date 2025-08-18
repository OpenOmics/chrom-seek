# Calling read coverage as peaks 
# ~~~~
# Single-end peak calling rules

# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
tmpdir                          = config['options']['tmp_dir']
seacr_dir                       = join(workpath, "SEACR")
bg_dir                          = join(workpath, "bedgraph")


wildcard_constraints:
    # - regex to avoid backslashes in name wild cards
    #   this corrects routing for `ppqt` and 
    #   `ppqt_tagalign` rules
    # - if this isn't included, ppqt file paths (because
    #   they're nested in the bam directory) get mixed up
    #   with regular bam paths
    # - also possible to use negative look around regex:
    #   "^((?!\/).)*$"
    name                                = "[A-Za-z0-9_-]+"


rule MACS2_broad:
    input:
        chip                                = join(bam_dir, "{name}.Q5DD_tagAlign.gz"),
        txt                                 = join(ppqt_dir, "{name}.Q5DD_tagAlign.ppqt.txt"),
        c_option                            = lambda w: join(bam_dir, f"{chip2input[w.name]}.Q5DD_tagAlign.gz")
                                                if w.name in chip2input and chip2input[w.name] != "" else [],
        frag_length                         = join(ppqt_dir, "{name}.Q5DD.fragment.length"),
    output:
        join(macsB_dir, "{name}", "{name}_peaks.broadPeak"),
    params:
        rname                               = 'MACS2_broad',
        gsize                               = config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver                             = config['tools']['MACSVER'],
        flag                                = lambda w: "-c" if w.name in chip2input else "",
        frag_len_script                     = join(bin_path, "ppqt_process.py"),
    shell: 
        """
        module load {params.macsver};
        ppqt_len=$(cat {input.frag_length})
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {macsB_dir}/{wildcards.name} \\
            --broad \\
            --broad-cutoff 0.01 \\
            --keep-dup="all" \\
            --nomodel \\
            --extsize ${{ppqt_len}}
        """


rule MACS2_narrow:
    input:
        chip                                = join(bam_dir, "{name}.Q5DD_tagAlign.gz"),
        txt                                 = join(ppqt_dir, "{name}.Q5DD.ppqt.txt"),
        c_option                            = lambda w: join(bam_dir, f"{chip2input[w.name]}.Q5DD_tagAlign.gz")
                                                if w.name in chip2input and chip2input[w.name] != "" else [],
        frag_length                         = join(ppqt_dir, "{name}.Q5DD.fragment.length"),
    output:
        join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"),
    params:
        rname                               = 'MACS2_narrow',
        gsize                               = config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver                             = config['tools']['MACSVER'],
        flag                                = lambda w: "-c" if chip2input[w.name] else "",
    shell: 
        """
        module load {params.macsver};
        ppqt_len=$(cat {input.frag_length})
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {macsN_dir}/{wildcards.name} \\
            -q 0.01 \\
            --keep-dup="all" \\
            --nomodel \\
            --extsize ${{ppqt_len}}
        """


rule SICER:
    input:
        chip                            = lambda w: join(bam_dir, w.name + ".Q5DD_tagAlign.gz"),
        fragLen                         = lambda w: join(ppqt_dir, w.name + ".Q5DD_tagAlign.ppqt.txt"),
        c_option                        = lambda w: join(bam_dir, f"{chip2input[w.name]}.Q5DD_tagAlign.gz"),
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

        mean_insert_size=$(awk '{{print $1}}' {input.fragLen})

        echo "printing out value of mean-insert-size ${{mean_insert_size}}"
        a={input.c_option}
        echo "Printing input.c_option ${{a}}"
        
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
            mv ${{tmp}}/{params.name}-W300-G600-islands-summary {params.this_sicer_dir}
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


rule bam2bedgraph:
    input: join(bam_dir, "{name}.Q5DD.bam")
    output: temp(join(bg_dir, "{name}.bedgraph"))
    container: config['images']['seacr']
    params:
        rname                           = 'bam2bedgraph',
        tmpdir                          = tmpdir,
        reference                       = config['references'][genome]['REFLEN'],
    threads:
        int(cluster['bam2bedgraph'].get('threads', cluster['__default__']['threads']))
    shell:
        """
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT
        export TMPDIR="${{tmp}}" # used by sort
        # bam to bedgraph
        samtools sort -@{threads} -n {input} > ${{tmp}}/name_sorted_sample.bam
        bedtools bamtobed -bedpe -i ${{tmp}}/name_sorted_sample.bam > ${{tmp}}/sample.bed
        awk '$1==$4 && $6-$2 < 1000 {{print $0}}' ${{tmp}}/sample.bed > ${{tmp}}/sample.clean.bed
        cut -f 1,2,6 ${{tmp}}/sample.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${{tmp}}/sample.fragments.bed
        bedtools genomecov -bg -i ${{tmp}}/sample.fragments.bed -g {params.reference} > {output}
        """


rule SEACR:
    input:
        exp                             = join(bg_dir, "{name}.bedgraph"),
        control                         = lambda w: join(bg_dir, f"{chip2input[w.name]}.bedgraph")
                                            if w.name in chip2input and chip2input[w.name] != "" else [],
    output:
        peaks                           = join(seacr_dir, "{name}", "{name}.stringent.bed")
    params:
        rname                           = 'SEACR',
        out_dir                         = join(seacr_dir, "{name}"),
        control_flag                    = lambda w, input: input.control if input.control else "0.01",
    container:
        config['images']['seacr']
    shell:
        """
        cd {params.out_dir}
        SEACR_1.3.sh {input.exp} {params.control_flag} non stringent {wildcards.name}
        """