import json
from os.path import join
from textwrap import dedent
from itertools import combinations
from shlex import quote
from scripts.common import allocated, mk_dir_if_not_exist
from scripts.peakcall import outputIDR, zip_peak_files, \
    calc_effective_genome_fraction, get_manorm_sizes
from scripts.grouping import test_for_block


######### workflow configuration value mapping #########
workpath                        = config["project"]["workpath"]
bin_path                        = config["project"]["binpath"]
genome                          = config["options"]["genome"]
blocks                          = config["project"]["blocks"]
groupdata                       = config["project"]["groups"]
contrast                        = config["project"]["contrast"]
uropaver                        = config["tools"]["UROPAVER"]
gtf                             = config["references"][genome]["GTFFILE"]
chips                           = config['project']['peaks']['chips']


######### Output Directory preparation #########
diffbind_dir2                   = join(workpath, "DiffBind_block")
log_dir                         = join(workpath, "logfiles")
local_log_dir                   = join(log_dir, "local")
diffbind_dir                    = join(workpath, "DiffBind")
peakqc_dir                      = join(workpath, "PeakQC")
diffbind_qc_dir                 = join(peakqc_dir, "DB_QC")
uropa_dir                       = join(workpath, "UROPA_annotations")
uropa_diffbind_dir              = join(uropa_dir, "DiffBind")
uropa_diffbind_join_dir         = join(workpath, "UROPA_DIFFBIND_TBLS")
bam_dir                         = join(workpath, "bam")
seacr_dir                       = join(workpath, "SEACR")
config_path                     = join(workpath, "config.json")


######### workflow configuration flags #########
group_combos                    = []
blocking                        = False if set(blocks.values()) in ({None}, {''}) else True
chip2input                      = config['project']['peaks']['inputs']
has_inputs                      = False if set(chip2input.values()) in ({''}, {None}) else True
blocking                        = False if set(blocks.values()) in ({None}, {""}) else True


######## CSV generation locally executed ########
localrules: diffbind_csv_macsN, diffbind_csv_macsB, diffbind_csv_genrich, diffbind_csv_SEACR


######### Diffbind binding analysis #########
rule diffbind_count:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
    output:
        peak_counts                     = join(
                                            diffbind_dir, 
                                            "{group1}_vs_{group2}-{PeakTool}", 
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"
                                          ),
        peak_list                       = join(
                                            diffbind_dir, 
                                            "{group1}_vs_{group2}-{PeakTool}", 
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_fullList.bed"
                                          ),
    params:
        rname                           = "diffbind_count",
        this_contrast                   = "{group1}_vs_{group2}",
        this_peaktool                   = "{PeakTool}",
        this_script                     = join(bin_path, "DiffBind_v2_load.R"),
    container: config["images"]["cfchip"]
    shell:
        dedent("""
        {params.this_script} \\
            --csvfile {input.csvfile} \\
            --counts {output.peak_counts} \\
            --list {output.peak_list} \\
            --peakcaller {params.this_peaktool}
        """)

        
######### Diffbind input preparation #########
rule diffbind_csv_macsN:
    input:
        bed                            = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips)
    output:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-macsNarrow",
                                            "{group1}_vs_{group2}-macsNarrow_Diffbind_prep.csv",
                                          ),
    params:
        rname                           = "diffbind_csv_macsN",
        this_peaktool                   = "macsNarrow",
        peakcaller                      = "narrowPeak",
        this_peakextension              = "_peaks.narrowPeak",
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        bam_dir                         = bam_dir,
        workpath                        = workpath,
        contrast                        = "{group1}_vs_{group2}",
    log: join(local_log_dir, "diffbind_csv_macsN", "{group1}_vs_{group2}_diffbind_csv.log")
    run:
            shell(dedent(
                """
                python {params.pythonscript} \\
                    --con {params.contrast} \\
                    --wp {params.workpath} \\
                    --pt {params.this_peaktool} \\
                    --pe {params.this_peakextension} \\
                    --bd {params.bam_dir} \\
                    --pc {params.peakcaller} \\
                    --csv {output.csvfile}
                """
            ))


rule diffbind_csv_genrich:
    input:
        bed                            = expand(join(genrich_dir, "{name}", "{name}.narrowPeak"), name=chips)
    output:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-Genrich",
                                            "{group1}_vs_{group2}-Genrich_Diffbind_prep.csv",
                                          ),
    params:
        rname                           = "diffbind_csv_genrich",
        this_peaktool                   = "Genrich",
        peakcaller                      = "narrowPeak",
        this_peakextension              = ".narrowPeak",
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        bam_dir                         = bam_dir,
        workpath                        = workpath,
        contrast                        = "{group1}_vs_{group2}",
    log: join(local_log_dir, "diffbind_csv_genrich", "{group1}_vs_{group2}_diffbind_csv.log")
    shell:
        dedent("""
        python {params.pythonscript} \\
            --con {params.contrast} \\
            --wp {params.workpath} \\
            --pt {params.this_peaktool} \\
            --pe {params.this_peakextension} \\
            --bd {params.bam_dir} \\
            --pc {params.peakcaller} \\
            --csv {output.csvfile} &> {log}
        """)


rule diffbind_csv_macsB:
    input:
        beds                            = expand(join(macsB_dir, "{name}", "{name}_peaks.broadPeak"), name=chips),
    output:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-macsBroad",
                                            "{group1}_vs_{group2}-macsBroad_Diffbind_prep.csv",
                                          ),
    params:
        rname                           = "diffbind_csv_macsB",
        this_peaktool                   = "macsBroad",
        peakcaller                      = "narrowPeak",
        this_peakextension              = "_peaks.broadPeak",
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        bam_dir                         = bam_dir,
        workpath                        = workpath,
        contrast                        = "{group1}_vs_{group2}",
    log: join(local_log_dir, "diffbind_csv_macsB", "{group1}_vs_{group2}_diffbind_csv.log")
    run:
            shell(dedent(
                """
                python {params.pythonscript} \\
                    --con {params.contrast} \\
                    --wp {params.workpath} \\
                    --pt {params.this_peaktool} \\
                    --pe {params.this_peakextension} \\
                    --bd {params.bam_dir} \\
                    --pc {params.peakcaller} \\
                    --csv {output.csvfile}
                """
            ))


rule diffbind_csv_SEACR:
    input:
        beds                            = expand(join(seacr_dir, "{name}", "{name}.stringent.bed"), name=chips),
    output:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-SEACR",
                                            "{group1}_vs_{group2}-SEACR_Diffbind_prep.csv",
                                          ),
    params:
        rname                           = "diffbind_csv_SEACR",
        this_peaktool                   = "SEACR",
        peakcaller                      = "raw",
        this_peakextension              = ".stringent.bed",
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        bam_dir                         = bam_dir,
        workpath                        = workpath,
        contrast                        = "{group1}_vs_{group2}",
    log: join(local_log_dir, "diffbind_csv_SEACR", "{group1}_vs_{group2}_diffbind_csv.log")
    shell: 
        dedent("""
            python {params.pythonscript} \\
                --con {params.contrast} \\
                --wp {params.workpath} \\
                --pt {params.this_peaktool} \\
                --pe {params.this_peakextension} \\
                --bd {params.bam_dir} \\
                --pc {params.peakcaller} \\
                --csv {output.csvfile}
        """)


######### diffbind EdgeR DE analysis #########
rule diffbind_edger:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        peak_counts                     = join(
                                            diffbind_dir, 
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_counts.rds"
                                          ),
    output:
        diffbind_report                 = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_EdgeR.html",
                                          ),
        up_file                         = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_EdgeR_up.bed",
                                          ),
        down_file                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_EdgeR_down.bed",
                                          ),
        peak_list                       = join(
                                           diffbind_dir,
                                           "{contrast}-{PeakTool}",
                                           "{contrast}-{PeakTool}_Diffbind_EdgeR_peak_list.tab",
                                        )
    params:
        rname                           = "diffbind_edger",
        rscript                         = join(bin_path, "DiffBind_v2_EdgeR.Rmd"),
        outdir                          = join(diffbind_dir, "{contrast}-{PeakTool}"),
    container:
        config["images"]["diffbind"]
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\" ]; then mkdir -p \"{tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\")
        trap 'rm -rf "${{tmp}}"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}

        cat <<'EOF' > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.diffbind_report}", 
            params=list(csvfile="{input.csvfile}", peakcaller="{wildcards.PeakTool}", list_file="{output.peak_list}", 
            up_file="{output.up_file}", down_file="{output.down_file}", contrasts="{wildcards.contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        ls -al ${{tmp}}
        sh ${{tmp}}/rscript.sh
        """)


rule diffbind_edger_blocking:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        peak_counts                     = join(
                                            diffbind_dir, 
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_counts.rds"),
    output:
        diffbind_block_report           = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_EdgeR.html",
                                          ),
        up_file                         = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_EdgeR_up.bed",
                                          ),
        down_file                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_EdgeR_down.bed",
                                          ),

        peak_list                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_EdgeR_peak_list.tab",
                                          ),
    params:
        rname                           = "diffbind_edger_block",
        blocking_rscript                = join(bin_path, "DiffBind_v2_EdgeR_block.Rmd"),
        outdir                          = join(diffbind_dir, "{contrast}-{PeakTool}"),
    container:
        config["images"]["diffbind"]
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\" ]; then mkdir -p \"{tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\")
        trap 'rm -rf "${{tmp}}"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}

        cat << EOF > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.blocking_rscript}", output_file="{output.diffbind_block_report}",
            params=list(csvfile="{input.csvfile}", peakcaller="{wildcards.PeakTool}", list_file="{output.peak_list}",
            up_file="{output.up_file}", down_file="{output.down_file}", contrasts="{wildcards.contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        ls -al ${{tmp}}
        sh ${{tmp}}/rscript.sh
        """)


######### diffbind DeSeq2 DE analysis #########
rule diffbind_deseq:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        peak_counts                     = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_counts.rds"
                                          ),
    output:
        diffbind_report                 = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_DeSeq2.html",
                                          ),
        up_file                         = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_DeSeq2_up.bed",
                                          ),
        down_file                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_DeSeq2_down.bed",
                                          ),
        peak_list                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_DeSeq2_peak_list.tab",
                                          ),
    params:
        rname                           = "diffbind_deseq2",
        rscript                         = join(bin_path, "DiffBind_v2_Deseq2.Rmd"),
        outdir                          = join(diffbind_dir, "{contrast}-{PeakTool}"),
    container:
        config["images"]["diffbind"]
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\\diffbind_deseq\" ]; then mkdir -p \"{tmpdir}\\diffbind_deseq\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\\diffbind_deseq")
        trap 'rm -rf "{tmpdir}\\diffbind_deseq"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}

        cat <<'EOF' > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.diffbind_report}", 
            params=list(csvfile="{input.csvfile}", peakcaller="{wildcards.PeakTool}", list_file="{output.peak_list}", 
            up_file="{output.up_file}", down_file="{output.down_file}", contrasts="{wildcards.contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        ls -al ${{tmp}}
        sh ${{tmp}}/rscript.sh
        """)


rule diffbind_deseq_blocking:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        peak_counts                     = join(
                                            diffbind_dir, 
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_counts.rds"
                                          ),
    output:
        diffbind_block_report           = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_DeSeq2.html",
                                          ),
        up_file                         = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_DeSeq2_up.bed",
                                          ),
        down_file                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_DeSeq2_down.bed",
                                          ),
        peak_list                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_DeSeq2_peak_list.tab",
                                          ),
    params:
        rname                           = "diffbind_deseq_block",
        blocking_rscript                = join(bin_path, "DiffBind_v2_Deseq2_block.Rmd"),
        outdir                          = join(diffbind_dir, "{contrast}-{PeakTool}"),
    container:
        config["images"]["diffbind"]
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\" ]; then mkdir -p \"{tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\")
        trap 'rm -rf "${{tmp}}"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}
        cat <<'EOF' > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.blocking_rscript}", output_file="{output.diffbind_block_report}", 
            params=list(csvfile="{input.csvfile}", peakcaller="{wildcards.PeakTool}", list_file="{output.peak_list}", 
            up_file="{output.up_file}", down_file="{output.down_file}", contrasts="{wildcards.contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        sh ${{tmp}}/rscript.sh
        """)


######### diffbindQC rules #########
rule diffbindQC_macsN:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips) if chips else [],
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips if chip2input[sample_name] != ''] if chips and has_inputs else [],
        samples_peaks               = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips) if chips else []
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_TMMcounts.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_macsN",
        peak_tool                   = "macsNarrow",
        peak_type                   = "narrow",
        all_bams                    = lambda w, input: list(set(input.sample_bams + input.control_bams)),
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-macsNarrow"),
        configfile                  = config_path,    
        pythonscript                = join(bin_path, "prep_diffbindQC.py"),
    container:
       config['images']['diffbind']
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\" ]; then mkdir -p \"{tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\")
        trap 'rm -rf "${{tmp}}"' EXIT

        {params.pythonscript} \\
            -p {input.samples_peaks} \\
            -t {params.peak_type} \\
            -o {output.csvfile} \\
            -c {params.configfile} \\
            -b {params.all_bams}

        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", counts_bed="{output.countsbed}", 
            counts_csv="{output.countscsv}", peakcaller="{params.peak_tool}"))'
        """)


rule diffbindQC_macsB:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips if chip2input[sample_name] != ''] if chips and has_inputs else [],
        samples_peaks               = expand(join(macsB_dir, "{name}", "{name}_peaks.broadPeak"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_TMMcounts.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_macsB",
        peak_tool                   = "macsBroad",
        peak_type                   = "narrow",
        all_bams                    = lambda w, input: list(set(input.sample_bams + input.control_bams)),
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-macsBroad"),
        configfile                  = config_path,
        pythonscript                = join(bin_path, "prep_diffbindQC.py"),
    container:
       config['images']['diffbind']
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\" ]; then mkdir -p \"{tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\")
        trap 'rm -rf "${{tmp}}"' EXIT

        {params.pythonscript} \\
            -p {input.samples_peaks} \\
            -t {params.peak_type} \\
            -o {output.csvfile} \\
            -c {params.configfile} \\
            -b {params.all_bams}

        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", counts_bed="{output.countsbed}", 
            counts_csv="{output.countscsv}", peakcaller="{params.peak_tool}"))'
        """)


rule diffbindQC_genrich:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips if chip2input[sample_name] != ''] if chips and has_inputs else [],
        samples_peaks               = expand(join(genrich_dir, "{name}", "{name}.narrowPeak"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_TMMcounts.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_genrich",
        peak_tool                   = "Genrich",
        peak_type                   = "narrow",
        all_bams                    = lambda w, input: list(set(input.sample_bams + input.control_bams)),
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-Genrich"),
        configfile                  = config_path,
        pythonscript                = join(bin_path, "prep_diffbindQC.py")
    container:
       config['images']['diffbind']
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\" ]; then mkdir -p \"{tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\")
        trap 'rm -rf "${{tmp}}"' EXIT
        
        {params.pythonscript} \\
            -p {input.samples_peaks} \\
            -t {params.peak_type} \\
            -o {output.csvfile} \\
            -c {params.configfile} \\
            -b {params.all_bams}

        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", counts_bed="{output.countsbed}", 
            counts_csv="{output.countscsv}", peakcaller="{params.peak_tool}"))'
        """)


rule diffbindQC_SEACR:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips] if chips and has_inputs else [],
        samples_peaks               = expand(join(seacr_dir, "{name}", "{name}.stringent.bed"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_TMMcounts.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_SEACR",
        peak_tool                   = "SEACR",
        peak_type                   = "raw",
        all_bams                    = lambda w, input: list(set(input.sample_bams + input.control_bams)),
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-SEACR"),
        pythonscript                = join(bin_path, "prep_diffbindQC.py"),
        configfile                  = config_path,
        groups                      = json.dumps(config['project']['groups']).replace('"', '\\"')
    container:
       config['images']['diffbind']
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\" ]; then mkdir -p \"{tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\")
        trap 'rm -rf "${{tmp}}"' EXIT

        {params.pythonscript} \\
            -p {input.samples_peaks} \\
            -t {params.peak_type} \\
            -o {output.csvfile} \\
            -c {params.configfile} \\
            -b {params.all_bams}

        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", counts_bed="{output.countsbed}", 
            counts_csv="{output.countscsv}", peakcaller="{params.peak_tool}"))'
        """)


######### UROPA + Diffbind #########
rule join_diffbind_uropa:
    input:
        finalhits_txt                   = join(
                                            uropa_diffbind_dir, 
                                            "{contrast}-{PeakTool}-{differential_app}", 
                                            "{contrast}_{PeakTool}_{differential_app}_{_type}_uropa_finalhits.txt"
                                          ),
        peak_list                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_{differential_app}_peak_list.tab",
                                          ),
    output:
        joined_tbl                      = join(
                                            uropa_diffbind_join_dir, 
                                            "{contrast}-{PeakTool}-{differential_app}_{_type}_UROPA_DIFFBIND_JOIN.txt"
                                          )
    params:
        rname                           = "join_diffbind_uropa",
        this_script                     = join(bin_path, "merge_uropa_diffbind.py"),
        output_dir                      = uropa_diffbind_join_dir,
        fdr_cutoff                      = 0.05,
        foldchange_cutoff               = 0
    container: config["images"]["cfchip"]
    shell:
        dedent("""
        mkdir -p {params.output_dir}
        {params.this_script} \\
            --uropa {input.finalhits_txt} \\
            --diffbind {input.peak_list} \\
            --output {output.joined_tbl} \\
            --fdr {params.fdr_cutoff} \\
            --fold {params.foldchange_cutoff}
        """)
