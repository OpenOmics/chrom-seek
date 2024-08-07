# Differential binding analysis rules
# ~~~~
import os
import json
from os.path import join
from scripts.common import allocated, mk_dir_if_not_exist
from scripts.peakcall import outputIDR, zip_peak_files, calc_effective_genome_fraction, get_manorm_sizes
from scripts.grouping import test_for_block


# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
blocks                          = config['project']['blocks']
groupdata                       = config['project']['groups']
contrast                        = config['project']['contrast']
uropaver                        = config['tools']['UROPAVER']
gtf                             = config['references'][genome]['GTFFILE']

# ~~ directoriesxw
diffbind_dir_block              = join(workpath, "DiffBindBlock")
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
blocking                        = False if None in list(blocks.values()) else True
if reps == "yes": otherDirs.append(diffbind_dir)
mk_dir_if_not_exist(PeakTools + otherDirs)


# ~~ peak calling configuration and outputs
PeakToolsNG = [ tool for tool in PeakTools if tool != "gem" ]
PeakExtensions = {
    'macsNarrow': '_peaks.narrowPeak',
    'macsBroad': '_peaks.broadPeak',
    'sicer': '_broadpeaks.bed',
    'gem': '.GEM_events.narrowPeak' ,
    'MANorm': '_all_MA.bed',
    'DiffbindEdgeR': '_Diffbind_EdgeR.bed',
    'DiffbindDeseq2': '_Diffbind_Deseq2.bed', 
    'DiffbindEdgeRBlock': '_Diffbind_EdgeR_block.bed',
    'DiffbindDeseq2Block': '_Diffbind_Deseq2_block.bed',
    'Genrich': '.narrowPeak',
    'DiffBindQC': '_DiffBindQC_TMMcounts.bed'
}

FileTypesDiffBind = { 
    'macsNarrow': 'narrowPeak',
    'macsBroad': 'narrowPeak',
    'sicer': 'bed', 
    'gem': 'narrowPeak',
    'Genrich': 'narrowPeak'
}

PeakExtensionsIDR = { 
    'macsNarrow': '_peaks.narrowPeak',
    'macsBroad': '_peaks.broadPeak',
    'sicer': '_sicer.broadPeak'
}

FileTypesIDR = { 
    'macsNarrow': 'narrowPeak',
    'macsBroad': 'broadPeak',
    'sicer': 'broadPeak'
}

RankColIDR = { 
    'macsNarrow': 'q.value',
    'macsBroad': 'q.value',
    'sicer': 'q.value'
}
IDRgroup, IDRsample1, IDRsample2, IDRpeaktool =	outputIDR(groupswreps, groupdata, chip2input, PeakToolsNG)
zipSample, zipTool, zipExt = zip_peak_files(chips, PeakTools, PeakExtensions)
contrastBlock = test_for_block(groupdata, contrast, blocks)
zipGroup1B, zipGroup2B, zipToolCB, contrastsB = zip_contrasts(contrastBlock, PeakTools)

# ~~ rules 
rule diffbind:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ]
    output:
        html                            = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind.html"),
        Deseq2                          = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2.bed"),
        EdgeR                           = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.bed"),
        EdgeR_txt                       = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.txt"),
        Deseq2_txt                      = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2.txt"),
        EdgeR_ftxt                      = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_fullList.txt"),
        Deseq2_ftxt                     = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_fullList.txt"),
        html_block                      = provided(join(uropa_diffbind_dir, "{group1}_vs_{group2}-{PeakTool}", "{group1}_vs_{group2}-{PeakTool}_Diffbind_blocking.html"), blocking)
    params:
        # variables and wildcards used in the shell directive
        rname                           = "diffbind",
        group1                          = "{group1}",
        group2                          = "{group2}",
        this_peaktool                   = "{PeakTool}",
        this_contrast                   = "{group1}_vs_{group2}",
        this_peakextension              = lambda w: PeakExtensions[w.PeakTool],
        peakcaller                      = lambda w: FileTypesDiffBind[w.PeakTool],
        # scripts in the bin directory used in the shell directive
        rscript                         = join(bin_path, "DiffBind_v2_ChIPseq.Rmd"),
        pythonscript                    = join(bin_path, "prep_diffbind.py"),
        blocking_rscript                = join(bin_path, "DiffBind_v2_ChIPseq_block.Rmd"),
        # output base directories or full file locations
        outdir                          = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}"),
        csvfile                         = join(
                                            diffbind_dir, 
                                            "{group1}_vs_{group2}-{PeakTool}",
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv"
                                        ),
        Deseq2_block                    = provided(join(
                                            diffbind_dir_block, 
                                            "{group1}_vs_{group2}-{PeakTool}", 
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_block.bed"
                                        ), blocking),
        EdgeR_block                     = provided(join(
                                            diffbind_dir_block, 
                                            "{group1}_vs_{group2}-{PeakTool}", 
                                            "{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_block.bed"
                                        ), blocking),
        outdir_block                    = join(diffbind_dir_block, "{group1}_vs_{group2}-{PeakTool}"),
    container:
        config['images']['cfchip']
    shell: 
        """
        python {params.pythonscript} --g1 {params.group1} --g2 {params.group2} --wp {workpath} \
            --pt {params.this_peaktool} --pe {params.this_peakextension} --bd {bam_dir} \
            --pc {params.peakcaller} --csv {params.csvfile}
        cp {params.rscript} {params.outdir}
        cd {params.outdir}
        Rscript -e 'rmarkdown::render("DiffBind_v2_ChIPseq.Rmd", output_file= "{output.html}", \
            params=list(csvfile="{params.csvfile}", contrasts="{params.this_contrast}", peakcaller="{params.this_peaktool}"))'
        if [ ! -f {output.Deseq2} ]; then touch {output.Deseq2}; fi
        if [ ! -f {output.EdgeR} ]; then touch {output.EdgeR}; fi

        if [ '"""+str(blocking)+"""' == True ]; then
            echo "DiffBind with Blocking"
            Rscript -e 'rmarkdown::render("{params.blocking_rscript}", output_file= "{output.html_block}", \
                params=list(csvfile= "{params.csvfile}", contrasts="{params.this_contrast}", peakcaller="{params.this_peaktool}", dir="{params.outdir_block}"))'
            if [ ! -f {params.Deseq2_block} ]; then touch {params.Deseq2_block}; fi
            if [ ! -f {params.EdgeR_block} ]; then touch {params.EdgeR_block}; fi
        fi
        """


rule UROPA:
    input:
        lambda w: join(workpath, w.PeakTool1, w.name, w.name + PeakExtensions[w.PeakTool2]),
    output:
        txt                             = join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_allhits.txt'),
        bed1                            = temp(join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_allhits.bed')),
        bed2                            = temp(join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_finalhits.bed')),
        json                            = join(uropa_dir, '{PeakTool1}', '{name}.{PeakTool2}.{type}.json'),
    params:
        rname                           = "uropa",
        fldr                            = join(uropa_dir, '{PeakTool1}'),
        outroot                         = join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}'),
    threads: 4,
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
        shell.prefix(f"module load {uropaver};")
        shell("uropa -i " + output.json + " -p " + params.outroot + " -t " + str(threads) + " -s")


rule manorm:
    input:
        bam1                            = lambda w: join(bam_dir, groupdata[w.group1][0] + ".Q5DD.bam"),
        bam2                            = lambda w: join(bam_dir, groupdata[w.group2][0] + ".Q5DD.bam"),
        # ppqt                            = join(ppqt_dir, "Q5DD.ppqt.txt"), # ppqt input into manorm TODO
        peak1                           = lambda w: join(workpath, w.tool, groupdata[w.group1][0], groupdata[w.group1][0] + PeakExtensions[w.tool]),
        peak2                           = lambda w: join(workpath, w.tool, groupdata[w.group2][0], groupdata[w.group2][0] + PeakExtensions[w.tool]),
    output:
        xls                             = join(manorm_dir, "{group1}_vs_{group2}-{tool}", "{group1}_vs_{group2}-{tool}_all_MAvalues.xls"),
        bed                             = temp(join(manorm_dir, "{group1}_vs_{group2}-{tool}", "{group1}_vs_{group2}-{tool}_all_MA.bed")),
        wigA                            = join(manorm_dir, "{group1}_vs_{group2}-{tool}", "output_tracks", "{group1}_vs_{group2}_A_values.wig.gz"),
        wigM                            = join(manorm_dir, "{group1}_vs_{group2}-{tool}", "output_tracks", "{group1}_vs_{group2}_M_values.wig.gz"),
        wigP                            = join(manorm_dir, "{group1}_vs_{group2}-{tool}", "output_tracks", "{group1}_vs_{group2}_P_values.wig.gz"),
    params:
        rname                           = 'manorm',
        fldr                            = join(manorm_dir, "{group1}_vs_{group2}-{tool}"),
        bedtoolsver                     = config['tools']['BEDTOOLSVER'],
        manormver                       = "manorm/1.1.4",
        extsizes                        = lambda w, input: get_manorm_sizes(w.group1, w.group2, groupdata, "")
    shell:
        """
        if [ ! -e /lscratch/$SLURM_JOBID ]; then 
            mkdir /lscratch/$SLURM_JOBID
        fi
        cd /lscratch/$SLURM_JOBID
        module load {params.manormver}
        module load {params.bedtoolsver}
        bamToBed -i {input.bam1} > bam1.bed
        bamToBed -i {input.bam2} > bam2.bed
        cut -f 1,2,3 {input.peak1} > peak1.bed
        cut -f 1,2,3 {input.peak2} > peak2.bed
        manorm \
            --p1 peak1.bed \
            --p2 peak2.bed \
            --r1 bam1.bed \
            --r2 bam2.bed \
            {params.extsizes} \
            -o {params.fldr} \
            --name1 {wildcards.group1} \
            --name2 {wildcards.group2}
        gzip {params.fldr}/output_tracks/*wig
        mv {params.fldr}/{wildcards.group1}_vs_{wildcards.group2}_all_MAvalues.xls {output.xls}
        tail -n +2 {output.xls} | nl -w2 | awk -v OFS='\t' '{{print $2,$3,$4,$9$1,$6}}' > {output.bed}
        """