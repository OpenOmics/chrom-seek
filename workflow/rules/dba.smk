# Differential binding analysis rules
# ~~~~
from os.path import join
import os
from scripts.common import allocated, mk_dir_if_not_exist
from scripts.peakcall import outputIDR, zip_peak_files, calc_effective_genome_fraction
from scripts.blocking import test_for_block


# ~~ workflow configuration
workpath                        = config['project']['workpath']
genome                          = config['options']['genome']
blocks                          = config['project']['blocks']
groupdata                       = config['project']['groups']


# ~~ directories
bin_path                        = join(workpath, "workflow", "bin")
diffbind_dir_block              = join(workpath, "DiffBindBlock")
diffbind_dir2                   = join(workpath, "DiffBind_block")
diffbind_dir                    = join(workpath, "DiffBind")
bam_dir                         = join(workpath, "bam")
qc_dir                          = join(workpath, "PeakQC")
idr_dir                         = join(workpath, "IDR")
memechip_dir                    = join(workpath, "MEME")
homer_dir                       = join(workpath, "HOMER_motifs")
uropa_dir                       = join(workpath, "UROPA_annotations")
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
        html = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind.html"),
        Deseq2 = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2.bed"),
        EdgeR = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.bed"),
        EdgeR_txt = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.txt"),
        Deseq2_txt = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2.txt"),
        EdgeR_ftxt = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_fullList.txt"),
        Deseq2_ftxt = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_fullList.txt"),
        html_block = provided(join(diffbind_dir_block, "{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_blocking.html"), blocking)
    params:
        rname = "diffbind",
        rscript = join(workpath, "workflow", "scripts","DiffBind_v2_ChIPseq.Rmd"),
        outdir    = join(diffbind_dir, "{group1}_vs_{group2}-{PeakTool}"),
        contrast  = "{group1}_vs_{group2}",
        csvfile   = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_prep.csv"),
        pythonscript = join(workpath,"workflow","scripts","prep_diffbind.py"),
        PeakExtension= lambda w: PeakExtensions[w.PeakTool],
        peakcaller= lambda w: FileTypesDiffBind[w.PeakTool],
        group1="{group1}",
        group2="{group2}",
        PeakTool="{PeakTool}",
        blocking=blocking,
        blocking_rscript = join(workpath,"workflow","scripts","DiffBind_v2_ChIPseq_block.Rmd"),
        outdir_block= join(workpath,diffbind_dir_block,"{group1}_vs_{group2}-{PeakTool}"),
        Deseq2_block = provided(join(workpath, diffbind_dir_block,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_block.bed"), blocking),
        EdgeR_block = provided(join(workpath, diffbind_dir_block,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_block.bed"), blocking),
    container:
        config['images']['cfchip']
    shell: """
    python {params.pythonscript} --g1 {params.group1} --g2 {params.group2} --wp {workpath} \
         --pt {params.PeakTool} --pe {params.PeakExtension} --bd {bam_dir} \
         --pc {params.peakcaller} --csv {params.csvfile}
    cp {params.rscript} {params.outdir}
    cd {params.outdir}
    Rscript -e 'rmarkdown::render("DiffBind_v2_ChIPseq.Rmd", output_file= "{output.html}", 
    params=list(csvfile= "{params.csvfile}", contrasts= "{params.contrast}", peakcaller= "{params.PeakTool}"))'
    if [ ! -f {output.Deseq2} ]; then touch {output.Deseq2}; fi
    if [ ! -f {output.EdgeR} ]; then touch {output.EdgeR}; fi

    if [ '{params.blocking}' == True ]; then
        echo "DiffBind with Blocking"
        Rscript -e 'rmarkdown::render("{params.blocking_rscript}", output_file= "{output.html_block}", 
        params=list(csvfile= "{params.csvfile}", contrasts= "{params.contrast}", peakcaller= "{params.PeakTool}", dir= "{params.outdir_block}"))'
        if [ ! -f {params.Deseq2_block} ]; then touch {params.Deseq2_block}; fi
        if [ ! -f {params.EdgeR_block} ]; then touch {params.EdgeR_block}; fi
    fi
    """


if assay == "cfchip":
    rule UROPA:
        input:
            lambda w: [ join(workpath, w.PeakTool1, w.name, w.name + PeakExtensions[w.PeakTool2]) ]
        output:
            txt=join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_allhits.txt'),
            bed1=temp(join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_allhits.bed')),
            bed2=temp(join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_finalhits.bed')),
        params:
            rname="uropa",
            uropaver = config['tools']['UROPAVER'],
            fldr = join(uropa_dir, '{PeakTool1}'),
            json = join(uropa_dir, '{PeakTool1}','{name}.{PeakTool2}.{type}.json'),
            outroot = join(uropa_dir, '{PeakTool1}','{name}_{PeakTool2}_uropa_{type}'),
            gtf = config['references'][genome]['GTFFILE'],
            threads = 4,
        shell: """
        module load {params.uropaver};
        # Dynamically creates UROPA config file
        if [ ! -e {params.fldr} ]; then mkdir {params.fldr}; fi
        echo '{{"queries":[ ' > {params.json}
        if [ '{wildcards.type}' == 'protTSS' ]; then
             echo '      {{ "feature":"gene","distance":3000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start" }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":10000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start" }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":100000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start" }}],' >> {params.json}
        fi
        echo '"show_attributes":["gene_id", "gene_name","gene_type"],' >> {params.json}
        echo '"priority":"Yes",' >> {params.json}
        echo '"gtf":"{params.gtf}",' >> {params.json}
        echo '"bed": "{input}" }}' >> {params.json}
        uropa -i {params.json} -p {params.outroot} -t {params.threads} -s
        """
else:
    rule UROPA:
        input:
            lambda w: [ join(workpath, w.PeakTool1, w.name, w.name + PeakExtensions[w.PeakTool2]) ]
        output:
            txt=join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_allhits.txt'),
            bed1=temp(join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_allhits.bed')),
            bed2=temp(join(uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_finalhits.bed')),
        params:
            rname="uropa",
            uropaver = config['tools']['UROPAVER'],
            fldr = join(uropa_dir, '{PeakTool1}'),
            json = join(uropa_dir, '{PeakTool1}','{name}.{PeakTool2}.{type}.json'),
            outroot = join(uropa_dir, '{PeakTool1}','{name}_{PeakTool2}_uropa_{type}'),
            gtf = config['references'][genome]['GTFFILE'],
            threads = 4,
        shell: """
        module load {params.uropaver};
        # Dynamically creates UROPA config file
        if [ ! -e {params.fldr} ]; then mkdir {params.fldr}; fi
        echo '{{"queries":[ ' > {params.json}
        if [ '{wildcards.type}' == 'prot' ]; then
             echo '      {{ "feature":"gene","distance":5000,"filter.attribute":"gene_type","attribute.value":"protein_coding" }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":100000,"filter.attribute":"gene_type","attribute.value":"protein_coding" }}],' >> {params.json}
        elif [ '{wildcards.type}' == 'genes' ]; then
             echo '      {{ "feature":"gene","distance":5000 }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":100000 }}],' >> {params.json}
        elif [ '{wildcards.type}' == 'protSEC' ]; then
             echo '      {{ "feature":"gene","distance":[3000,1000],"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start" }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":3000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"end" }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":100000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"center" }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":100000,"filter.attribute":"gene_type","attribute.value":"protein_coding" }}],' >> {params.json}
        elif [ '{wildcards.type}' == 'protTSS' ]; then
             echo '      {{ "feature":"gene","distance":[3000,1000],"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start" }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":10000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start" }},' >> {params.json}
             echo '      {{ "feature":"gene","distance":100000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start" }}],' >> {params.json}

        fi
        echo '"show_attributes":["gene_id", "gene_name","gene_type"],' >> {params.json}
        echo '"priority":"Yes",' >> {params.json}
        echo '"gtf":"{params.gtf}",' >> {params.json}
        echo '"bed": "{input}" }}' >> {params.json}
        uropa -i {params.json} -p {params.outroot} -t {params.threads} -s
        """

rule manorm:
    input: 
        bam1 = lambda w: join(workpath,bam_dir, groupdata[w.group1][0] + ".Q5DD.bam"),
        bam2 = lambda w: join(workpath,bam_dir, groupdata[w.group2][0] + ".Q5DD.bam"),
        ppqt = join(workpath,bam_dir, "Q5DD.ppqt.txt"),
        peak1 = lambda w: join(workpath, w.tool, groupdata[w.group1][0], groupdata[w.group1][0] + PeakExtensions[w.tool]),
        peak2 = lambda w: join(workpath, w.tool, groupdata[w.group2][0], groupdata[w.group2][0] + PeakExtensions[w.tool]),
    output:
        xls = join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}","{group1}_vs_{group2}-{tool}_all_MAvalues.xls"),
        bed = temp(join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}","{group1}_vs_{group2}-{tool}_all_MA.bed")),
        wigA = join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}","output_tracks","{group1}_vs_{group2}_A_values.wig.gz"),
        wigM = join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}","output_tracks","{group1}_vs_{group2}_M_values.wig.gz"),
        wigP = join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}","output_tracks","{group1}_vs_{group2}_P_values.wig.gz"),
    params:
        rname='manorm',
        fldr = join(workpath,manorm_dir,"{group1}_vs_{group2}-{tool}"),
        bedtoolsver=config['tools']['BEDTOOLSVER'],
        sample1= lambda w: groupdata[w.group1][0],
        sample2= lambda w: groupdata[w.group2][0],
        manormver="manorm/1.1.4"
    run:
        commoncmd1 = "if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi "
        commoncmd2 = "cd /lscratch/$SLURM_JOBID; "
        commoncmd3 = "module load {params.manormver}; module load {params.bedtoolsver}; "
        cmd1 = "bamToBed -i {input.bam1} > bam1.bed; "
        cmd2 = "bamToBed -i {input.bam2} > bam2.bed; "
        cmd3 = "cut -f 1,2,3 {input.peak1} > peak1.bed; "
        cmd4 = "cut -f 1,2,3 {input.peak2} > peak2.bed; "
        file=list(map(lambda z:z.strip().split(),open(input.ppqt,'r').readlines()))
        extsize1 = [ ppqt[1] for ppqt in file if ppqt[0] == params.sample1 ][0]
        extsize2 = [ ppqt[1] for ppqt in file if ppqt[0] == params.sample2 ][0]
        cmd5 = "manorm --p1 peak1.bed --p2 peak2.bed --r1 bam1.bed --r2 bam2.bed --s1 " + extsize1  + " --s2 " + extsize2 + " -o {params.fldr} --name1 '" + wildcards.group1 + "' --name2 '" + wildcards.group2 + "'; "
        cmd6 = "gzip {params.fldr}/output_tracks/*wig; "
        cmd7 = "mv {params.fldr}/" + wildcards.group1 + "_vs_" + wildcards.group2 + "_all_MAvalues.xls {output.xls}; "
        cmd8 = "tail -n +2 {output.xls} | nl -w2 | awk -v OFS='\t' '{{print $2,$3,$4,$9$1,$6}}' > {output.bed}"
        shell(commoncmd1)
        shell( commoncmd2 + commoncmd3 + cmd1 + cmd2 + cmd3 + cmd4 + cmd5 + cmd6 + cmd7 + cmd8 )