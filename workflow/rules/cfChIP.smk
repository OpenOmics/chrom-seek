# cell-free ChIP-seq
# ~~~~
# rules: picard_dedup, cfChIPtool, cfChIPcompile


# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
blocks                          = config['project']['blocks']
groupdata                       = config['project']['groups']

# Directory end points
bam_dir                         = join(workpath, "bam")
cfTool_dir                      = join(workpath, "cfChIPtool")
cfTool_subdir2                  = join(cfTool_dir, "BED", "H3K4me3")
qc_dir                          = join(workpath, "QC")


rule cfChIPtool:
    input: 
        out5                    = join(bam_dir, "{name}.Q5DD.bam.idxstat"),
    output:
        out1                    = join(cfTool_subdir2, "{name}.Q5DD.tagAlign.gz"),
        out2                    = join(cfTool_dir, "Output", "H3K4me3", "Signatures", "{name}.Q5DD.csv"),
        out3                    = join(cfTool_dir, "Samples", "H3K4me3", "{name}.Q5DD.rdata"),
    params:
        rname                   = 'cfChiP',
        rver                    = "R/4.1.0",
        toolkit                 = config['references'][genome]['cfChIP_TOOLS_SRC'],
        tmpfile                 = lambda w: join(cfTool_subdir2, w.name + ".Q5DD.tagAlign"),
        tag                     = lambda w: temp(join(bam_dir, w.name + ".Q5DD_tagAlign"))
    container:
        config['images']['cfchip']
    shell: 
        """
        cp {params.tag} {params.tmpfile}
        gzip {params.tmpfile}

        Rscript {params.toolkit}/ProcessBEDFiles.R \\
            -a {params.toolkit}/SetupFiles/H3K4me3 \\
            -r {cfTool_dir} \\
            -p {cfTool_dir} \\
            -m H3K4me3 \\
            -S {output.out1}
        """


rule cfChIPcompile:
    input:
        expand(join(cfTool_dir, "Output", "H3K4me3", "Signatures", "{name}.Q5DD.csv"), name=chips)
    output:
        txt                     = join(qc_dir, "H3K4me3_cfChIP_signature.txt"),
        pdf                     = join(qc_dir, "H3K4me3_cfChIP_signature.pdf")
    params:
        rname                   = "cfChIP2",
        script                  = join(bin_path, "cfChIP_signatures.R"),
        infolder                = join(cfTool_dir, "Output", "H3K4me3", "Signatures"),
    container:
        config['images']['cfchip']
    shell: 
        """
        Rscript -e "source('{params.script}'); mergeSignatures( '{params.infolder}', '{output.txt}' )";
        Rscript -e "source('{params.script}'); plotSignatures( '{output.txt}', '{output.pdf}' )";
        """


rule promoterTable1:
    input:
        expand(join(workpath,uropa_dir,'{PeakTool}','{name}_{PeakTool}_uropa_protTSS_allhits.txt'), PeakTool=PeakTools, name=chips),
    output:
        txt                     = join(uropa_dir, "promoterTable1", "{PeakTool}_promoter_overlap_summaryTable.txt")
    params:
        rname                   = "promoter1",
        script                  = join(bin_path, "promoterAnnotation_by_Gene.R"),
        infolder                = join(uropa_dir, '{PeakTool}')
    container:
        config['images']['cfchip']
    shell: 
        """
        Rscript -e "source('{params.script}'); peakcallVersion('{params.infolder}','{output.txt}')";
        """


rule promoterTable2:
    input:
        expand(join(diffbind_dir, '{name}_DiffbindDeseq2_uropa_protTSS_allhits.txt'), name=contrasts),
    output:
        txt                     = join(workpath,uropa_dir,"promoterTable2",'DiffbindDeseq2_{PeakTool}_promoter_overlap_summaryTable.txt'),
    params:
        rname                   = "promoter2",
        script1                 = join(bin_path, "promoterAnnotation_by_Gene.R"),
        script2                 = join(bin_path, "significantPathways.R"),
        infolder                = workpath,
        gtf                     = config['references'][genome]['GTFFILE'],
    container:
        config['images']['cfchip']
    shell: 
        """
        Rscript -e "source('{params.script1}'); diffbindVersion('{params.infolder}','{output.txt}')";
        Rscript -e "source('{params.script2}'); promoterAnnotationWrapper('{output.txt}','{params.gtf}','KEGG')";
        Rscript -e "source('{params.script2}'); promoterAnnotationWrapper('{output.txt}','{params.gtf}','Reactome')";
        """

rule diffbindQC:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ]
    output:
        html                    = join(qc_dir, "AllSamples-{PeakTool}", "AllSamples-{PeakTool}_DiffBindQC.html"),
        bed                     = join(workpath, "QC", "AllSamples-{PeakTool}", "AllSamples-{PeakTool}_DiffBindQC_TMMcounts.bed"),
    params:
       rname                    = "diffbindQC",
       contrast                 = "AllSamples",
       PeakTool                 = "{PeakTool}",
       rscript                  = join(bin_path, "DiffBind_v2_cfChIP_QC.Rmd"),
       outdir                   = join(qc_dir, "AllSamples-{PeakTool}"),
       csvfile                  = join(qc_dir, "AllSamples-{PeakTool}", "AllSamples-{PeakTool}_DiffBind_prep.csv"),
       pythonscript             = join(bin_path, "prep_diffbindQC.py"),
       PeakExtension            = lambda w: PeakExtensions[w.PeakTool],
       peakcaller               = lambda w: FileTypesDiffBind[w.PeakTool],
    container:
       config['images']['cfchip']
    shell: 
        """
        python {params.pythonscript} --wp {workpath} \
            --pt {params.PeakTool} --pe {params.PeakExtension} --bd {bam_dir} \
            --pc {params.peakcaller} --csv {params.csvfile}
        cp {params.rscript} {params.outdir}
        cd {params.outdir}
        Rscript -e 'rmarkdown::render("DiffBind_v2_cfChIP_QC.Rmd", output_file= "{output.html}", 
            params=list(csvfile= "{params.csvfile}", contrasts= "{params.contrast}", peakcaller= "{params.PeakTool}"))'
        """