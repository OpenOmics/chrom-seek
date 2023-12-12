# TODO: This Snakefile needs to be completely refactored.
# Python standard library
from os.path import join
import os

# Local imports
from scripts.common import (
    allocated
)

def outputIDR(groupswreps, groupdata, chip2input, tools):
    """
    Produces the correct output files for IDR. All supposed replicates
    should be directly compared when possible using IDR. IDR malfunctions
    with bed files and GEM so it will not run with either of those.
    Because there is no q-value calculated for SICER when there is no 
    input file, those samples are also ignored.
    """
    IDRgroup, IDRsample1, IDRsample2, IDRpeaktool = [], [], [], []
    for group in groupswreps:
        nsamples = len(groupdata[group])
        for i in range(nsamples):
            ctrlTF = chip2input[groupdata[group][i]] != ""
            for j in range(i+1,nsamples):
                if ctrlTF == (chip2input[groupdata[group][j]] != ""):
                    if ctrlTF == False:
                        tooltmp = [ tool for tool in tools if tool != "sicer" ]
                    else:
                        tooltmp = tools			           
                    IDRgroup.extend([group] * len(tooltmp))
                    IDRsample1.extend([groupdata[group][i]] * len(tooltmp))
                    IDRsample2.extend([groupdata[group][j]] * len(tooltmp))
                    IDRpeaktool.extend(tooltmp)
    return( IDRgroup, IDRsample1, IDRsample2, IDRpeaktool )


def zip_peak_files(chips, PeakTools, PeakExtensions):
    """Making input file names for FRiP"""
    zipSample, zipTool, zipExt = [], [], []
    for chip in chips:
        for PeakTool in PeakTools:
            zipSample.append(chip)
            zipTool.append(PeakTool)
            zipExt.append(PeakExtensions[PeakTool])
    return(zipSample, zipTool, zipExt)


def calc_effective_genome_fraction(effectivesize, genomefile):
    """
    calculate the effective genome fraction by calculating the
    actual genome size from a .genome-like file and then dividing
    the effective genome size by that number
    """
    lines=list(map(lambda x:x.strip().split("\t"),open(genomefile).readlines()))
    genomelen=0
    for chrom,l in lines:
        if not "_" in chrom and chrom!="chrX" and chrom!="chrM" and chrom!="chrY":
            genomelen+=int(l)
    return(str(float(effectivesize)/ genomelen))



# PREPARING TO DEAL WITH A VARIED SET OF PEAKCALL TOOLS
gem_dir = "gem"
macsB_dir = "macsBroad"
sicer_dir = "sicer"

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
    'Genrich': '.narrowPeak'
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


# CREATING DIRECTORIES
bam_dir='bam'
qc_dir='PeakQC'
idr_dir = 'IDR'
memechip_dir = "MEME"
homer_dir = "HOMER_motifs"
manorm_dir = "MANorm"
downstream_dir = "Downstream"

otherDirs = [qc_dir, homer_dir, uropa_dir]
if reps == "yes":
    # otherDirs.append(idr_dir)
    otherDirs.append(diffbind_dir)

for d in PeakTools + otherDirs:
        if not os.path.exists(join(workpath,d)):
                os.mkdir(join(workpath,d))


# Blocking code
diffbind_dir2 = "DiffBind_block"
blocks=config['project']['blocks']

def test_for_block(contrast, blocks):
   """ only want to run blocking on contrasts where all
   individuals are on both sides of the contrast """
   contrastBlock = [ ]
   for con in contrast:
       group1 = con[0]
       group2 = con[1]
       block1 = [ blocks[sample] for sample in groupdata[group1] ]
       block2 = [ blocks[sample] for sample in groupdata[group2] ]
       if len(block1) == len(block2):
           if len(set(block1).intersection(block2)) == len(block1):
                contrastBlock.append(con)
   return contrastBlock  


contrastBlock = test_for_block(contrast,blocks)
zipGroup1B, zipGroup2B, zipToolCB, contrastsB = zip_contrasts(contrastBlock, PeakTools)

rule diffbind:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ]
    output:
        html = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind.html"),
        Deseq2 = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2.bed"),
        EdgeR = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.bed"),
        EdgeR_txt = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR.txt"),
        Deseq2_txt = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2.txt"),
        EdgeR_ftxt = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_EdgeR_fullList.txt"),
        Deseq2_ftxt = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_Deseq2_fullList.txt"),
        html_block = provided(join(workpath,diffbind_dir_block,"{group1}_vs_{group2}-{PeakTool}","{group1}_vs_{group2}-{PeakTool}_Diffbind_blocking.html"), blocking)
    params:
        rname="diffbind",
        rscript = join(workpath,"workflow","scripts","DiffBind_v2_ChIPseq.Rmd"),
        outdir    = join(workpath,diffbind_dir,"{group1}_vs_{group2}-{PeakTool}"),
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
            join(workpath, uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_allhits.txt')
        params:
            rname="uropa",
            uropaver = config['tools']['UROPAVER'],
            fldr = join(workpath, uropa_dir, '{PeakTool1}'),
            json = join(workpath, uropa_dir, '{PeakTool1}','{name}.{PeakTool2}.{type}.json'),
            outroot = join(workpath, uropa_dir, '{PeakTool1}','{name}_{PeakTool2}_uropa_{type}'),
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
            join(workpath, uropa_dir, '{PeakTool1}', '{name}_{PeakTool2}_uropa_{type}_allhits.txt')
        params:
            rname="uropa",
            uropaver = config['tools']['UROPAVER'],
            fldr = join(workpath, uropa_dir, '{PeakTool1}'),
            json = join(workpath, uropa_dir, '{PeakTool1}','{name}.{PeakTool2}.{type}.json'),
            outroot = join(workpath, uropa_dir, '{PeakTool1}','{name}_{PeakTool2}_uropa_{type}'),
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