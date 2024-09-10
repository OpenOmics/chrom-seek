# Differential binding analysis rules
# ~~~~
import os
from os.path import join
from textwrap import dedent
from scripts.common import allocated, mk_dir_if_not_exist
from scripts.peakcall import outputIDR, zip_peak_files, \
    calc_effective_genome_fraction, get_manorm_sizes
from scripts.grouping import test_for_block


# ~~ workflow configuration
workpath                        = config["project"]["workpath"]
bin_path                        = config["project"]["binpath"]
genome                          = config["options"]["genome"]
blocks                          = config["project"]["blocks"]
groupdata                       = config["project"]["groups"]
contrast                        = config["project"]["contrast"]
uropaver                        = config["tools"]["UROPAVER"]
gtf                             = config["references"][genome]["GTFFILE"]


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
    "DiffbindEdgeR": "_Diffbind_EdgeR_full_list.txt",
    "DiffbindDeseq2": "_Diffbind_Deseq2_full_list.txt",
    "DiffbindEdgeRBlock": "_Diffbind_block_EdgeR_full_list.txt",
    "DiffbindDeseq2Block": "_Diffbind_block_Deseq2_full_list.txt",
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
    "macsBroad": "_peaks.broadPeak",
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

localrules: UROPA_prep_in


# ~~ rules
rule diffbind_csv:
    input:
        bed                             = lambda w: [
                                            join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool])
                                            for chip in chips
                                          ],
    output:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
    params:
        rname                           = "diffbind_csv",
        group1                          = "{group1}",
        group2                          = "{group2}",
        this_peaktool                   = "{PeakTool}",
        peakcaller                      = lambda w: FileTypesDiffBind[w.PeakTool],
        this_peakextension              = lambda w: PeakExtensions[w.PeakTool],
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        bam_dir                         = bam_dir,
        workpath                        = workpath,
    container:
        config["images"]["cfchip"]
    shell:
        """
        python {params.pythonscript} \
            --g1 {params.group1} \
            --g2 {params.group2} \
            --wp {params.workpath} \
            --pt {params.this_peaktool} \
            --pe {params.this_peakextension} \
            --bd {params.bam_dir} \
            --pc {params.peakcaller} \
            --csv {output.csvfile}
        """


rule diffbind_count:
    input:
        bed                             = lambda w: [
                                            join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool])
                                            for chip in chips
                                          ],
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
    output:
        peak_counts                     = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"),
        peak_list                       = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}_Diffbind_fullList.bed"),
    params:
        rname                           = "diffbind_count",
        this_contrast                   = "{group1}_vs_{group2}",
        this_peaktool                   = "{PeakTool}",
        this_script                     = join(bin_path, "DiffBind_v2_load.R"),
    threads: 32
    container:
        config["images"]["cfchip"]
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
        peak_counts                     = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"),
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
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_full_list.txt",
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
        peak_counts                     = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"),
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
        peak_counts                     = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"),
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
        peak_counts                     = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}_Diffbind_counts.rds"),
    output:
        diffbind_block_report           = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_blocking_EdgeR.html",
                                          ),
        up_file                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_block_EdgeR_up.bed",
                                          ),
        down_file                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_block_EdgeR_down.bed",
                                          ),
        full_list                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_block_EdgeR_full_list.txt",
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
            up_file="{output.up_file}", down_file="{output.down_file}", contrasts="{params.this_contrast}", counts="{input.peak_counts}"))'
        EOF

        chmod +x ${{tmp}}/rscript.sh
        echo "--"
        cat ${{tmp}}/rscript.sh
        echo "--"
        ls -al ${{tmp}}
        sh ${{tmp}}/rscript.sh
        """)


rule UROPA_prep_in:
    input:
        lambda w: join(workpath, w.PeakTool1, w.name, w.name + PeakExtensions[w.PeakTool2])
    params:
        fldr                            = join(uropa_dir, "{PeakTool1}"),
    log:
        join(uropa_dir, "{PeakTool1}", "{name}.{PeakTool2}.log")
    output:
        this_json                       = [
                                            join(
                                                uropa_dir, 
                                                "{PeakTool1}", 
                                                "{name}.{PeakTool2}."+pktype+".json"
                                            ) for pktype in peak_types
                                          ],
        reformat_bed                    = [
                                            join(
                                                uropa_dir, 
                                                "{PeakTool1}", 
                                                "{name}_{PeakTool2}_uropa_"+pktype+"_input.bed"
                                            ) for pktype in peak_types
                                          ],
    run:
        import json, csv
        if not os.path.exists(params.fldr): 
            os.mkdir(params.fldr, mode=0o775)
        csv_map = {
            'seqnames': 'chr',
            'start': 'start',
            'end': 'stop',
            'strand': 'strand'
        }
        bed_map = {
            0: 'chr',
            1: 'start',
            2: 'stop',
        }
        base_query = {
            "feature": "gene",
            "filter.attribute": "gene_type",
            "attribute.value": "protein_coding", 
            "feature.anchor": "start" 
        }

        assert len(input) == len(peak_types), 'Unequal # of inputs and peak types, something wrong!'
        for i, peak_type in enumerate(peak_types):
            json_construct = dict()
            json_construct['queries'] = []
            json_construct['show_attributes'] = ["gene_id", "gene_name", "gene_type"]
            json_construct["priority"] = "Yes"
            json_construct['gtf'] = gtf
            
            # reformat to standard bed
            sniff = open(input[i]).read(len('seqnames'))
            rdr = csv.DictReader(open(input[i]), delimiter='\t')
            newbed = []
            for row in rdr:
                # two different formats for bed file input 
                # need to back up somewhere see if theres an upstream method for 
                # fixing the discrepancy of formatting
                if sniff == 'seqnames':
                    # diff bind format
                    row = {csv_map[k.lower()]: v for k, v in row.items() if k.lower() in csv_map}
                else:
                    # ?? format
                    row = {bed_map[i]: v for i, (k, v) in enumerate(row.items()) if i < 3}
            with open(output.reformat_bed[i], 'w') as fo:
                wrt = csv.DictWriter(fo, fieldnames=csv_map.values())
                # wrt.writeheader() # bed file wo header row
                wrt.writerows(newbed)
            json_construct['bed'] = output.reformat_bed[i]

            if assay == 'cfchip':
                if peak_type == 'protTSS':
                    for _d in (3000, 10000, 100000):
                        this_q = base_query.copy()
                        this_q['distance'] = _d
                        json_construct['queries'].append(this_q)
            else:
                if peak_type == 'prot':
                    for _d in (5000, 100000):
                        this_q = base_query.copy()
                        del this_q["feature.anchor"]
                        this_q['distance'] = _d
                        json_construct['queries'].append(this_q)
                elif peak_type == 'genes':
                    this_query = {}
                    this_query['feature'] = 'gene'
                    for _d in (5000, 100000):
                        this_q = base_query.copy()
                        del this_q["feature.anchor"]
                        del this_q["filter.attribute"]
                        del this_q["attribute.value"]
                        this_q['distance'] = _d
                        json_construct['queries'].append(this_q)
                elif peak_type == 'protSEC':
                    # distance, feature.anchor
                    query_values = (
                        ([3000, 1000], "start"), 
                        (3000,         "end"), 
                        (100000,       "center"), 
                        (100000,       None)
                    )
                    for _distance, feature_anchor in query_values:
                        this_q = base_query.copy()
                        del this_q["feature.anchor"]
                        if feature_anchor: 
                            this_q["feature.anchor"] = feature_anchor
                        this_q['distance'] = _distance
                        json_construct['queries'].append(this_q)
                elif peak_type == 'protTSS':
                    for _d in ([3000, 1000], 10000, 100000):
                        this_q = base_query.copy()
                        this_q['distance'] = _d
                        json_construct['queries'].append(this_q)

            with open(output.this_json[i], 'w') as jo:
                json.dump(json_construct, jo, indent=4)
                jo.close()


rule UROPA:
    input:
        this_json                       = join(uropa_dir, "{PeakTool1}", "{name}.{PeakTool2}.{_type}.json"),
    output:
        txt                             = join(
                                            uropa_dir, 
                                            "{PeakTool1}", 
                                            "{name}_{PeakTool2}_uropa_{_type}_allhits.txt"
                                          ),
        bed1                            = temp(join(
                                            uropa_dir, 
                                            "{PeakTool1}", 
                                            "{name}_{PeakTool2}_uropa_{_type}_allhits.bed"
                                          )),
        bed2                            = temp(join(
                                            uropa_dir,
                                            "{PeakTool1}",
                                            "{name}_{PeakTool2}_uropa_{_type}_finalhits.bed",
                                          )),
        log                             = join(uropa_dir, "{PeakTool1}", "{name}_{PeakTool2}_uropa_{_type}.log")
    params:
        rname                           = "uropa",
        outroot                         = join(uropa_dir, "{PeakTool1}", "{name}_{PeakTool2}_uropa_{_type}"),
    threads: int(allocated("threads", "UROPA", cluster))
    container:
        config["images"]["uropa"]
    shell:
        """
        uropa -i {input.this_json} -p {params.outroot} -l {output.log} -t {threads} -s
        """