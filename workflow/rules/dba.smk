# Differential binding analysis rules
# ~~~~
import os
import json
from os.path import join
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
    "DiffbindEdgeR": "_Diffbind_EdgeR.bed",
    "DiffbindDeseq2": "_Diffbind_Deseq2.bed",
    "DiffbindEdgeRBlock": "_Diffbind_EdgeR_block.bed",
    "DiffbindDeseq2Block": "_Diffbind_Deseq2_block.bed",
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


# ~~ rules
rule init_diffbind:
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
        rname                           = "initialize_diffbind",
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


rule diffbind:
    input:
        csvfile                         = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        bed                             = lambda w: [
                                            join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool])
                                            for chip in chips
                                          ],
    output:
        diffbind_report                 = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind.html",
                                          ),
        Deseq2                          = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2.bed",
                                          ),
        EdgeR                           = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.bed",
                                          ),
        EdgeR_txt                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.txt",
                                          ),
        Deseq2_txt                      = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2.txt",
                                          ),
        EdgeR_ftxt                      = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_fullList.txt",
                                          ),
        Deseq2_ftxt                     = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_fullList.txt",
                                          ),
        # consensus_pks                   =  join(
        #                                     diffbind_dir, 
        #                                     "{group1}_vs_{group2}-{PeakTool}",
        #                                     "{group1}_vs_{group2}-{PeakTool}_Diffbind_consensusPeaks.bed"
        #                                   ),
    params:
        rname                           = "diffbind",
        this_peaktool                   = "{PeakTool}",
        this_contrast                   = "{group1}_vs_{group2}",
        this_peakextension              = lambda w: PeakExtensions[w.PeakTool],
        peakcaller                      = lambda w: FileTypesDiffBind[w.PeakTool],
        rscript                         = join(bin_path, "DiffBind_v2_ChIPseq.Rmd"),
        outdir                          = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}"),
    container:
        config["images"]["cfchip"]
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", \
                                        output_file="{output.diffbind_report}", \
                                        params=list( \
                                            csvfile="{input.csvfile}", \
                                            contrasts="{params.this_contrast}", \
                                            peakcaller="{params.this_peaktool}" \
                                        ) \
                                     )'
        """


rule diffbind_blocking:
    input:
        csvfile                        = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv",
                                          ),
        bed                             = lambda w: [
                                            join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool])
                                            for chip in chips
                                          ],
    output:
        diffbind_block_report           = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_blocking.html",
                                          ),
        Deseq2                          = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_block.bed",
                                          ),
        EdgeR                           = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_block.bed",
                                          ),
        EdgeR_txt                       = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_block.txt",
                                          ),
        Deseq2_txt                      = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_block.txt",
                                          ),
        EdgeR_ftxt                      = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_fullList_block.txt",
                                          ),
        Deseq2_ftxt                     = join(
                                            diffbind_dir,
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_fullList_block.txt",
                                          ),
        # consensus_pks                   =  join(
        #                                     diffbind_dir, 
        #                                     "{group1}_vs_{group2}-{PeakTool}",
        #                                     "{group1}_vs_{group2}-{PeakTool}_Diffbind_consensusPeaks_block.bed"
        #                                   ),
    params:
        rname                           = "diffbind_block",
        blocking_rscript                = join(bin_path, "DiffBind_v2_ChIPseq_block.Rmd"),
        outdir                          = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}"),
        this_peaktool                   = "{PeakTool}",
        this_contrast                   = "{group1}_vs_{group2}",
    container:
        config["images"]["cfchip"]
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.blocking_rscript}", \
                                        output_file="{output.diffbind_block_report}", \
                                        params=list( \
                                            csvfile="{input.csvfile}", \
                                            contrasts="{params.this_contrast}", \
                                            peakcaller="{params.this_peaktool}" \
                                        ) \
                                     )'
        """


localrules: UROPA_prep_in


rule UROPA_prep_in:
    input:
        lambda w: join(workpath, w.PeakTool1, w.name, w.name + PeakExtensions[w.PeakTool2]),
    params:
        fldr                            = join(uropa_dir, "{PeakTool1}"),
    output:
        json                            = join(uropa_dir, "{PeakTool1}", "{name}.{PeakTool2}.{type}.json"),
    run:
        # Dynamically creates UROPA config file
        if not os.path.exists(params.fldr): 
            os.mkdir(params.fldr, mode=0o775)

        json_construct = dict()
        json_construct['queries'] = []
        json_construct['show_attributes'] = ["gene_id", "gene_name", "gene_type"]
        json_construct["priority"] = "Yes"
        json_construct['gtf'] = gtf
        json_construct['bed'] = input[0]

        base_query = {
            "feature": "gene",
            "filter.attribute": "gene_type",
            "attribute.value": "protein_coding", 
            "feature.anchor": "start" 
        }

        if assay == 'cfchip':
            if wildcards.type == 'protTSS':
                for _d in (3000, 10000, 100000):
                    this_q = base_query.copy()
                    this_q['distance'] = _d
                    json_construct['queries'].append(this_q)
        else:
            if wildcards.type == 'prot':
                for _d in (5000, 100000):
                    this_q = base_query.copy()
                    del this_q["feature.anchor"]
                    this_q['distance'] = _d
                    json_construct['queries'].append(this_q)
            elif wildcards.type == 'genes':
                this_query = {}
                this_query['feature'] = 'gene'
                for _d in (5000, 100000):
                    this_q = base_query.copy()
                    del this_q["feature.anchor"]
                    del this_q["filter.attribute"]
                    del this_q["attribute.value"]
                    this_q['distance'] = _d
                    json_construct['queries'].append(this_q)
            elif wildcards.type == 'protSEC':
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
            elif wildcards.type == 'protTSS':
                for _d in ([3000, 1000], 10000, 100000):
                    this_q = base_query.copy()
                    this_q['distance'] = _d
                    json_construct['queries'].append(this_q)

        with open(output.json, 'w') as jo:
            json.dump(json_construct, jo, indent=4)
            jo.close()

        if not os.path.exists(output.json):
            raise FileNotFoundError(output.json + " does not exist!")



rule UROPA:
    input:
        json                            = join(uropa_dir, "{PeakTool1}", "{name}.{PeakTool2}.{type}.json"),
    output:
        txt                             = join(
                                            uropa_dir, 
                                            "{PeakTool1}", 
                                            "{name}_{PeakTool2}_uropa_{type}_allhits.txt"
                                          ),
        bed1                            = temp(join(
                                            uropa_dir, 
                                            "{PeakTool1}", 
                                            "{name}_{PeakTool2}_uropa_{type}_allhits.bed"
                                          )),
        bed2                            = temp(join(
                                            uropa_dir,
                                            "{PeakTool1}",
                                            "{name}_{PeakTool2}_uropa_{type}_finalhits.bed",
                                          )),
        log                             = join(uropa_dir, "{PeakTool1}", "{name}_{PeakTool2}_uropa_{type}.log")
    params:
        rname="uropa",
        outroot=join(uropa_dir, "{PeakTool1}", "{name}_{PeakTool2}_uropa_{type}"),
    threads: 4
    container:
        config["images"]["uropa"]
    shell:
        """
        uropa -i {input.json} -p {params.outroot} -l {output.log} -t {threads} -s
        """