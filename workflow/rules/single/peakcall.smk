# Calling read coverage as peaks 
# ~~~~
# Single-end peak calling rules

# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
tmpdir                          = config['options']['tmp_dir']


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
                                                if w.name and chip2input[w.name] else [],
    output:
        join(macsB_dir, "{name}", "{name}_peaks.broadPeak"),
    params:
        rname                               = 'MACS2_broad',
        gsize                               = config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver                             = config['tools']['MACSVER'],
        flag                                = lambda w: "-c" if chip2input[w.name] else "",
        frag_len_script                     = join(bin_path, "ppqt_process.py"),
    shell: 
        """
        module load {params.macsver};
        ppqt_len=$({params.frag_len_script} {input.txt})
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
        """


rule MACS2_narrow:
    input:
        chip                                = join(bam_dir, "{name}.Q5DD_tagAlign.gz"),
        txt                                 = join(ppqt_dir, "{name}.Q5DD.ppqt.txt"),
        c_option                            = lambda w: join(bam_dir, f"{chip2input[w.name]}.Q5DD_tagAlign.gz")
                                                if w.name and chip2input[w.name] else [],
    output:
        join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"),
    params:
        rname                               = 'MACS2_narrow',
        gsize                               = config['references'][genome]['EFFECTIVEGENOMESIZE'],
        macsver                             = config['tools']['MACSVER'],
        flag                                = lambda w: "-c" if chip2input[w.name] else "",
        frag_len_script                     = join(bin_path, "ppqt_process.py"),
    shell: 
        """
        module load {params.macsver};
        ppqt_len=$({params.frag_len_script} {input.txt})
        macs2 callpeak \\
            -t {input.chip} {params.flag} {input.c_option} \\
            -g {params.gsize} \\
            -n {wildcards.name} \\
            --outdir {macsN_dir}/{wildcards.name} \\
            -q 0.01 \\
            --keep-dup="all" \\
            --nomodel \\
            --extsize $ppqt_len
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


