# Differential binding analysis rules
# ~~~~
from os.path import join
from textwrap import dedent
from itertools import combinations
from scripts.common import allocated, mk_dir_if_not_exist
from scripts.peakcall import outputIDR, zip_peak_files, \
    calc_effective_genome_fraction, get_manorm_sizes
from scripts.grouping import test_for_block


# ~~ workflow configuration & directories ~~
workpath                        = config["project"]["workpath"]
bin_path                        = config["project"]["binpath"]
genome                          = config["options"]["genome"]
blocks                          = config["project"]["blocks"]
groupdata                       = config["project"]["groups"]
contrast                        = config["project"]["contrast"]
uropaver                        = config["tools"]["UROPAVER"]
gtf                             = config["references"][genome]["GTFFILE"]
chips                           = config['project']['peaks']['chips']
diffbind_dir2                   = join(workpath, "DiffBind_block")
diffbind_dir                    = join(workpath, "DiffBind")
uropa_dir                       = join(workpath, "UROPA_annotations")
uropa_diffbind_dir              = join(uropa_dir, "DiffBind")
bam_dir                         = join(workpath, "bam")
ppqt_dir                        = join(bam_dir, "ppqt")
qc_dir                          = join(workpath, "PeakQC")
idr_dir                         = join(workpath, "IDR")
memechip_dir                    = join(workpath, "MEME")
homer_dir                       = join(workpath, "HOMER_motifs")
manorm_dir                      = join(workpath, "MANorm")
downstream_dir                  = join(workpath, "Downstream")
otherDirs                       = [qc_dir, homer_dir, uropa_dir]
cfTool_dir                      = join(workpath, "cfChIPtool")
cfTool_subdir2                  = join(cfTool_dir, "BED", "H3K4me3")
group_combos                    = []

# ~~ workflow config ~~
for (g1, g2) in list(combinations(config['project']['groups'].keys(), 2)):
    group_combos.append(f"{g1}_vs_{g2}")
blocking = False if set(blocks.values()) in ({None}, {""}) else True
if reps == "yes": otherDirs.append(diffbind_dir)
mk_dir_if_not_exist(PeakTools + otherDirs)


localrules: UROPA_prep_in_macsB, UROPA_prep_in_macsN, \
                diffbind_csv_macsN, diffbind_csv_macsB, \
                UROPA_prep_in_diffbind


# ~~ differential binding analysis ~~ #
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
    threads: 6
    shell:
        dedent("""
        {params.this_script} \\
            --csvfile {input.csvfile} \\
            --counts {output.peak_counts} \\
            --list {output.peak_list} \\
            --peakcaller {params.this_peaktool}
        """)


# ~~  macs Narrow peak annotation  ~~ #
rule UROPA_prep_in_macsN:
    input: 
        join(macsN_dir, "{name}", "{name}_peaks.narrowPeak")
    params:
        rname                           = "UROPA_prep_in_macsN",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [join(
                                            uropa_dir, 
                                            "macsNarrow", 
                                            "{name}.macsNarrow." + pktype + ".json"
                                          ) for pktype in peak_types],
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -t {params.peak_types}
        """)


rule UROPA_macsN:
    input: join(uropa_dir, "macsNarrow", "{name}.macsNarrow.{_type}.json")
    params:
        rname                           = "UROPA_macsN",
        outroot                         = join(uropa_dir, "macsNarrow"),
    output: join(uropa_dir, "macsNarrow", "{name}_uropa_{_type}_allhits.txt"),
    threads: int(allocated("threads", "UROPA_macsNarrow", cluster)),
    log: join(uropa_dir, "macsNarrow", "{name}.macsNarrow.{_type}.log"),
    container: config["images"]["uropa"]
    shell: "uropa -i {input} -p {params.outroot} -l {log} -t {threads} -s"
        


# ~~  macs Broad peak annotation  ~~ #
rule UROPA_prep_in_macsB:
    input: 
        join(macsB_dir, "{name}", "{name}_peaks.broadPeak"),
    params:
        rname                           = "UROPA_prep_in_macsB",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [join(
                                            uropa_dir, 
                                            "macsBroad", 
                                            "{name}.macsBroad." + pktype + ".json"
                                          ) for pktype in peak_types],
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -t {params.peak_types}
        """)


rule UROPA_macsB:
    input: join(uropa_dir, "macsBroad", "{name}.macsBroad.{_type}.json")
    params:
        rname                           = "UROPA_macsB",
        outroot                         = join(uropa_dir, "macsBroad"),
    output: join(uropa_dir, "macsBroad", "{name}_uropa_{_type}_allhits.txt"),
    threads: int(allocated("threads", "UROPA_macsB", cluster)),
    log: join(uropa_dir, "macsBroad", "{name}.macsBroad.{_type}.log"),
    container: config["images"]["uropa"]
    shell: "uropa -i {input} -p {params.outroot} -l {log} -t {threads} -s"


# ~~ diffbind EdgeR DE analysis ~~ #
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
        full_list                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_EdgeR_full_list.bed",
                                          ),
    params:
        rname                           = "diffbind_edger",
        this_peakextension              = lambda w: {
                                                        "macsNarrow": "_peaks.narrowPeak",
                                                        "macsBroad": "_peaks.broadPeak",
                                                        "sicer": "_broadpeaks.bed",
                                                        "gem": ".GEM_events.narrowPeak",
                                                        "MANorm": "_all_MA.bed",
                                                        "DiffbindEdgeR": "_DiffbindEdgeR_fullList.bed",
                                                        "DiffbindDeseq2": "_DiffbindDeseq2_fullList.bed",
                                                        "DiffbindEdgeRBlock": "_DiffbindEdgeRBlock_fullList.bed",
                                                        "DiffbindDeseq2Block": "_DiffbindDeseq2Block_fullList.bed",
                                                        "Genrich": ".narrowPeak",
                                                        "DiffBindQC": "_DiffBindQC_TMMcounts.bed",
                                                    }[w.PeakTool],
        peakcaller                      = lambda w: {
                                                        "macsNarrow": "narrowPeak",
                                                        "macsBroad": "narrowPeak",
                                                        "sicer": "bed",
                                                        "gem": "narrowPeak",
                                                        "Genrich": "narrowPeak",
                                                    }[w.PeakTool],
        rscript                         = join(bin_path, "DiffBind_v2_EdgeR.Rmd"),
        outdir                          = join(diffbind_dir, "{contrast}-{PeakTool}"),
    container:
        config["images"]["cfchip"]
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
            params=list(csvfile="{input.csvfile}", peakcaller="{wildcards.PeakTool}", list_file="{output.full_list}", 
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
                                            "{contrast}-{PeakTool}_Diffbind_blocking_EdgeR.html",
                                          ),
        full_list                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_EdgeR_full_list.bed",
                                          ),
    params:
        rname                           = "diffbind_edger_block",
        blocking_rscript                = join(bin_path, "DiffBind_v2_EdgeR_block.Rmd"),
        outdir                          = join(diffbind_dir, "{contrast}-{PeakTool}"),
    container:
        config["images"]["cfchip"]
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
            params=list(csvfile="{input.csvfile}", peakcaller="{wildcards.PeakTool}", list_file="{output.full_list}",
            contrasts="{wildcards.contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        ls -al ${{tmp}}
        sh ${{tmp}}/rscript.sh
        """)


# ~~ diffbind DeSeq2 DE analysis ~~ # 
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
                                            "{contrast}-{PeakTool}_Diffbind_Deseq2_up.bed",
                                          ),
        down_file                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_Deseq2_down.bed",
                                          ),
        full_list                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_Deseq2_full_list.txt",
                                          ),
    params:
        rname                           = "diffbind_deseq2",
        this_peakextension              = lambda w: {
                                                        "macsNarrow": "_peaks.narrowPeak",
                                                        "macsBroad": "_peaks.broadPeak",
                                                        "sicer": "_broadpeaks.bed",
                                                        "gem": ".GEM_events.narrowPeak",
                                                        "MANorm": "_all_MA.bed",
                                                        "DiffbindEdgeR": "_DiffbindEdgeR_fullList.bed",
                                                        "DiffbindDeseq2": "_DiffbindDeseq2_fullList.bed",
                                                        "DiffbindEdgeRBlock": "_DiffbindEdgeRBlock_fullList.bed",
                                                        "DiffbindDeseq2Block": "_DiffbindDeseq2Block_fullList.bed",
                                                        "Genrich": ".narrowPeak",
                                                        "DiffBindQC": "_DiffBindQC_TMMcounts.bed",
                                                    }[w.PeakTool],
        peakcaller                      = lambda w: {
                                                        "macsNarrow": "narrowPeak",
                                                        "macsBroad": "narrowPeak",
                                                        "sicer": "bed",
                                                        "gem": "narrowPeak",
                                                        "Genrich": "narrowPeak",
                                                    }[w.PeakTool],
        rscript                         = join(bin_path, "DiffBind_v2_Deseq2.Rmd"),
        outdir                          = join(diffbind_dir, "{contrast}-{PeakTool}"),
    container:
        config["images"]["cfchip"]
    shell:
        dedent("""
        if [ ! -d \"{tmpdir}\" ]; then mkdir -p \"{tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{tmpdir}\")
        trap 'rm -rf "{tmpdir}"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}

        cat <<'EOF' > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.diffbind_report}", 
            params=list(csvfile="{input.csvfile}", peakcaller="{wildcards.PeakTool}", list_file="{output.full_list}", 
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
                                            "{contrast}-{PeakTool}_Diffbind_blocking_DeSeq2.html",
                                          ),
        full_list                       = join(
                                            diffbind_dir,
                                            "{contrast}-{PeakTool}",
                                            "{contrast}-{PeakTool}_Diffbind_block_Deseq2_full_list.txt",
                                          ),
    params:
        rname                           = "diffbind_deseq_block",
        blocking_rscript                = join(bin_path, "DiffBind_v2_Deseq2_block.Rmd"),
        outdir                          = join(diffbind_dir, "{contrast}-{PeakTool}"),
    container:
        config["images"]["cfchip"]
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
            params=list(csvfile="{input.csvfile}", peakcaller="{wildcards.PeakTool}", list_file="{output.full_list}", 
            contrasts="{wildcards.contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        sh ${{tmp}}/rscript.sh
        """)


# ~~  diffbind peak annotation  ~~ #
rule diffbind_csv_macsN:
    input:
        beds                            = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips),
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
    threads: 32
    run:
        for i, contrast in enumerate(group_combos):
            shell(dedent(
                """
                python {params.pythonscript} \\
                    --con "{wildcards.group1}_vs_{wildcards.group2}" \\
                    --wp {params.workpath} \\
                    --pt {params.this_peaktool} \\
                    --pe {params.this_peakextension} \\
                    --bd {params.bam_dir} \\
                    --pc {params.peakcaller} \\
                    --csv {output.csvfile}
                """
            ))


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
        peakcaller                      = "broadPeak",
        this_peakextension              = "_peaks.broadPeak",
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        bam_dir                         = bam_dir,
        workpath                        = workpath,
    threads: 32
    run:
        for i, contrast in enumerate(group_combos):
            shell(dedent(
                """
                python {params.pythonscript} \\
                    --con "{wildcards.group1}_vs_{wildcards.group2}" \\
                    --wp {params.workpath} \\
                    --pt {params.this_peaktool} \\
                    --pe {params.this_peakextension} \\
                    --bd {params.bam_dir} \\
                    --pc {params.peakcaller} \\
                    --csv {output.csvfile}
                """
            ))


rule UROPA_prep_in_diffbind:
    input: 
        join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind_fullList.bed")
    params:
        rname                           = "UROPA_prep_in_diffbind",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [join(
                                            uropa_dir, 
                                            "DiffBind", 
                                            "{group1}_vs_{group2}.{PeakTool}.DiffBind." + pktype + ".json"
                                          ) for pktype in peak_types],
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -t {params.peak_types}
        """)


rule UROPA_diffbind:
    input: join(uropa_dir, "DiffBind", "{group1}_vs_{group2}.{PeakTool}.DiffBind.{_type}.json")
    params:
        rname                           = "UROPA_diffbind",
        outroot                         = join(uropa_dir, "DiffBind"),
    output: join(uropa_dir, "DiffBind", "{group1}_vs_{group2}_{PeakTool}_uropa_{_type}_allhits.txt"),
    threads: int(allocated("threads", "UROPA_diffbind", cluster)),
    log: join(uropa_dir, "DiffBind", "{group1}_vs_{group2}.{PeakTool}.DiffBind.{_type}.log"),
    container: config["images"]["uropa"]
    shell: "uropa -i {input} -p {params.outroot} -l {log} -t {threads} -s"