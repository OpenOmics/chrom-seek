# Differential binding analysis rules
# ~~~~
import os
from os.path import join
from textwrap import dedent
from scripts.common import allocated, mk_dir_if_not_exist
from scripts.peakcall import outputIDR, zip_peak_files, \
    calc_effective_genome_fraction, get_manorm_sizes, getMacChip
from scripts.grouping import test_for_block
from itertools import combinations


# ~~ workflow configuration
workpath                        = config["project"]["workpath"]
bin_path                        = config["project"]["binpath"]
genome                          = config["options"]["genome"]
blocks                          = config["project"]["blocks"]
groupdata                       = config["project"]["groups"]
contrast                        = config["project"]["contrast"]
uropaver                        = config["tools"]["UROPAVER"]
gtf                             = config["references"][genome]["GTFFILE"]
chips                           = config['project']['peaks']['chips']

# ~~ directories
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
for (g1, g2) in list(combinations(config['project']['groups'].keys(), 2)):
    group_combos.append(f"{g1}_vs_{g2}")

# ~~ workflow switches
blocking = False if set(blocks.values()) in ({None}, {""}) else True
if reps == "yes": otherDirs.append(diffbind_dir)
mk_dir_if_not_exist(PeakTools + otherDirs)


# ~~ peak calling configuration and outputs
PeakToolsNG = [tool for tool in PeakTools if tool != "gem"]
PeakExtensions = {
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
}

FileTypesDiffBind = {
    "macsNarrow": "narrowPeak",
    "macsBroad": "narrowPeak",
    "sicer": "bed",
    "gem": "narrowPeak",
    "Genrich": "narrowPeak",
}

PeakExtensionsIDR = {
    "macsNarrow": "_peaks.narrowPeak",
    "macsBroad": "",
    "sicer": "_sicer.broadPeak",
}

FileTypesIDR = {
    "macsNarrow": "narrowPeak",
    "macsBroad": "broadPeak",
    "sicer": "broadPeak",
}

RankColIDR = {"macsNarrow": "q.value", "macsBroad": "q.value", "sicer": "q.value"}
IDRgroup, IDRsample1, IDRsample2, IDRpeaktool = outputIDR(
    groupswreps, groupdata, chip2input, PeakToolsNG
)
zipSample, zipTool, zipExt = zip_peak_files(chips, PeakTools, PeakExtensions)
contrastBlock = test_for_block(groupdata, contrast, blocks)
zipGroup1B, zipGroup2B, zipToolCB, contrastsB = zip_contrasts(contrastBlock, PeakTools)

localrules: UROPA_prep_in_db, UROPA_prep_in_macsN, UROPA_prep_in_macsB


# ~~ rules
rule diffbind_csv_macsB:
    input:
        beds                            = expand(join(macsB_dir, "{name}", "{name}_peaks.broadPeak"), name=chips),
    output:
        csvfile                         = expand(join(
                                            diffbind_dir,
                                            "{contrast}-macsNarrow",
                                            "{contrast}-macsNarrow_Diffbind_prep.csv",
                                          ), contrast=group_combos),
    params:
        rname                           = "diffbind_csv_macsB",
        this_peaktool                   = "macsBroad",
        peakcaller                      = "narrowPeak",
        this_peakextension              = "_peaks.broadPeak",
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        bam_dir                         = bam_dir,
        workpath                        = workpath,
    threads: 6
    run:
        for i, contrast in enumerate(group_combos):
            shell(dedent(
                """
                python {params.pythonscript} \\
                    --con """ + contrast + """ \\
                    --wp {params.workpath} \\
                    --pt {params.this_peaktool} \\
                    --pe {params.this_peakextension} \\
                    --bd {params.bam_dir} \\
                    --pc {params.peakcaller} \\
                    --csv {output.csvfile}
                """
            ))


rule diffbind_csv_macsN:
    input:
        beds                            = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips),
    output:
        csvfile                         = expand(join(
                                            diffbind_dir,
                                            "{contrast}-macsBroad",
                                            "{contrast}-macsBroad_Diffbind_prep.csv",
                                          ), contrast=group_combos),
                                          
    params:
        rname                           = "diffbind_csv_macsN",
        this_peaktool                   = "macsBroad",
        peakcaller                      = "narrowPeak",
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
                    --con """ + contrast + """ \\
                    --wp {params.workpath} \\
                    --pt {params.this_peaktool} \\
                    --pe {params.this_peakextension} \\
                    --bd {params.bam_dir} \\
                    --pc {params.peakcaller} \\
                    --csv {output.csvfile}
                """
            ))


rule diffbind_csv_deseq2:
    input:
        join(
            diffbind_dir,
            "{group1}_vs_{group2}-{PeakTool}",
            "{group1}_vs_{group2}-{PeakTool}_DiffbindDeseq2_full_list.txt",
        ),
    output:
        csvfile                         = expand(join(
                                            diffbind_dir,
                                            "{contrast}-DiffbindDeseq2",
                                            "{contrast}-DiffbindDeseq2_Diffbind_prep.csv",
                                          ), contrast=group_combos),                                
    params:
        rname                           = "diffbind_csv_deseq2",
        this_peaktool                   = "DiffbindDeseq2",
        peakcaller                      = "DiffBind",
        this_peakextension              = "_DiffbindDeseq2_fullList.bed",
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        bam_dir                         = bam_dir,
        workpath                        = workpath,
    threads: 32
    run:
        for i, contrast in enumerate(group_combos):
            shell(dedent(
                """
                python {params.pythonscript} \\
                    --con """ + contrast + """ \\
                    --wp {params.workpath} \\
                    --pt {params.this_peaktool} \\
                    --pe {params.this_peakextension} \\
                    --bd {params.bam_dir} \\
                    --pc {params.peakcaller} \\
                    --csv {output.csvfile}
                """
            ))


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


rule diffbind_edger:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        peak_counts                     = join(
                                            diffbind_dir, 
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"
                                          ),
    output:
        diffbind_report                 = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.html",
                                          ),
        up_file                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_up.bed",
                                          ),
        down_file                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_down.bed",
                                          ),
        full_list                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_full_list.bed",
                                          ),
    params:
        rname                           = "diffbind_edger",
        this_peaktool                   = "{PeakTool}",
        this_contrast                   = "{group1}_vs_{group2}",
        this_peakextension              = lambda w: PeakExtensions[w.PeakTool],
        peakcaller                      = lambda w: FileTypesDiffBind[w.PeakTool],
        rscript                         = join(bin_path, "DiffBind_v2_EdgeR.Rmd"),
        outdir                          = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}"),
    container:
        config["images"]["cfchip"]
    shell:
        dedent("""
        if [ ! -d \"""" + tmpdir + """\" ]; then mkdir -p \"""" + tmpdir + """\"; fi
        tmp=$(mktemp -d -p \"""" + tmpdir + """\")
        trap 'rm -rf "${{tmp}}"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}

        cat <<'EOF' > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.diffbind_report}", 
            params=list(csvfile="{input.csvfile}", peakcaller="{params.this_peaktool}", list_file="{output.full_list}", 
            up_file="{output.up_file}", down_file="{output.down_file}", contrasts="{params.this_contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        ls -al ${{tmp}}
        sh ${{tmp}}/rscript.sh
        """)


rule diffbind_deseq:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        peak_counts                     = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"
                                          ),
    output:
        diffbind_report                 = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_DeSeq2.html",
                                          ),
        up_file                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_up.bed",
                                          ),
        down_file                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_down.bed",
                                          ),
        full_list                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_full_list.txt",
                                          ),
    params:
        rname                           = "diffbind_deseq2",
        this_peaktool                   = "{PeakTool}",
        this_contrast                   = "{group1}_vs_{group2}",
        this_peakextension              = lambda w: PeakExtensions[w.PeakTool],
        peakcaller                      = lambda w: FileTypesDiffBind[w.PeakTool],
        rscript                         = join(bin_path, "DiffBind_v2_Deseq2.Rmd"),
        outdir                          = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}"),
    container:
        config["images"]["cfchip"]
    shell:
        dedent("""
        if [ ! -d \"""" + tmpdir + """\" ]; then mkdir -p \"""" + tmpdir + """\"; fi
        tmp=$(mktemp -d -p \"""" + tmpdir + """\")
        trap 'rm -rf "${{tmp}}"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}

        cat <<'EOF' > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.diffbind_report}", 
            params=list(csvfile="{input.csvfile}", peakcaller="{params.this_peaktool}", list_file="{output.full_list}", 
            up_file="{output.up_file}", down_file="{output.down_file}", contrasts="{params.this_contrast}", counts="{input.peak_counts}"))'
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
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        peak_counts                     = join(
                                            diffbind_dir, 
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"
                                          ),
    output:
        diffbind_block_report           = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_blocking_DeSeq2.html",
                                          ),
        full_list                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_block_Deseq2_full_list.txt",
                                          ),
    params:
        rname                           = "diffbind_deseq_block",
        blocking_rscript                = join(bin_path, "DiffBind_v2_Deseq2_block.Rmd"),
        outdir                          = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}"),
        this_peaktool                   = "{PeakTool}",
        this_contrast                   = "{group1}_vs_{group2}",
    container:
        config["images"]["cfchip"]
    shell:
        dedent("""
        if [ ! -d \"""" + tmpdir + """\" ]; then mkdir -p \"""" + tmpdir + """\"; fi
        tmp=$(mktemp -d -p \"""" + tmpdir + """\")
        trap 'rm -rf "${{tmp}}"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}
        cat <<'EOF' > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.blocking_rscript}", output_file="{output.diffbind_block_report}", 
            params=list(csvfile="{input.csvfile}", peakcaller="{params.this_peaktool}", list_file="{output.full_list}", 
            contrasts="{params.this_contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        sh ${{tmp}}/rscript.sh
        """)


rule diffbind_edger_blocking:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        peak_counts                     = join(
                                            diffbind_dir, 
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"),
    output:
        diffbind_block_report           = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_blocking_EdgeR.html",
                                          ),
        full_list                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_block_EdgeR_full_list.bed",
                                          ),
    params:
        rname                           = "diffbind_edger_block",
        blocking_rscript                = join(bin_path, "DiffBind_v2_EdgeR_block.Rmd"),
        outdir                          = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}"),
        this_peaktool                   = "{PeakTool}",
        this_contrast                   = "{group1}_vs_{group2}",
    container:
        config["images"]["cfchip"]
    shell:
        dedent("""
        if [ ! -d \"""" + tmpdir + """\" ]; then mkdir -p \"""" + tmpdir + """\"; fi
        tmp=$(mktemp -d -p \"""" + tmpdir + """\")
        trap 'rm -rf "${{tmp}}"' EXIT

        mkdir -p {params.outdir}
        cd {params.outdir}

        cat << EOF > ${{tmp}}/rscript.sh
        #!/bin/bash
        Rscript -e 'rmarkdown::render("{params.blocking_rscript}", output_file="{output.diffbind_block_report}",
            params=list(csvfile="{input.csvfile}", peakcaller="{params.this_peaktool}", list_file="{output.full_list}",
            contrasts="{params.this_contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        ls -al ${{tmp}}
        sh ${{tmp}}/rscript.sh
        """)



rule UROPA_prep_in_edger:
    input:
        join(diffbind_dir, "{group1}_vs_{group2}-DiffbindEdgeR", "{group1}_vs_{group2}-DiffbindEdgeR_Diffbind_fullList.bed"),
    params:
        rname                           = "UROPA_prep_in_edger",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [
                                            join(
                                                uropa_dir, 
                                                "DiffbindEdgeR", 
                                                "{group1}_vs_{group2}.DiffbindEdgeR."+pktype+".json"
                                            ) for pktype in peak_types
                                          ]
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.assay} \\
            -t {params.peak_types}
        """)


rule UROPA_prep_in_deseq2:
    input:
        join(diffbind_dir, "{group1}_vs_{group2}-DiffbindDeseq2", "{group1}_vs_{group2}-DiffbindDeseq2_Diffbind_fullList.bed"),
    params:
        rname                           = "UROPA_prep_in_deseq2",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [
                                            join(
                                                uropa_dir, 
                                                "DiffbindDeseq2", 
                                                "{group1}_vs_{group2}.DiffbindDeseq2."+pktype+".json"
                                            ) for pktype in peak_types
                                          ]
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.assay} \\
            -t {params.peak_types}
        """)


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
        this_json                       = [
                                            join(
                                                uropa_dir, 
                                                "macsBroad", 
                                                "{name}.macsBroad."+pktype+".json"
                                            ) for pktype in peak_types
                                          ],
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -t {params.peak_types}
        """)


rule UROPA_prep_in_macsN:
    input:
        join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"),
    params:
        rname                           = "UROPA_prep_in_macsN",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [
                                            join(
                                                uropa_dir, 
                                                "macsNarrow", 
                                                "{name}.macsNarrow."+pktype+".json"
                                            ) for pktype in peak_types
                                          ],
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -t {params.peak_types}
        """)


rule UROPA_macsNarrow:
    input:
        join(
            uropa_dir, 
            "macsNarrow", 
            "{name}.macsNarrow.{_type}.json"
        )
        