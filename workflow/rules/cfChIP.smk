# cell-free DNA ChIP-seq rules:
#   - picard_dedup 
#   - cfChIPtool
#   - cfChIPcompile


rule cfChIPtool:
    input: 
        out5=join(workpath,bam_dir,"{name}.Q5DD.bam.idxstat"),
    output:
        out1=join(workpath,cfTool_subdir2,"{name}.Q5DD.tagAlign.gz"),
        out2=join(workpath,cfTool_dir,"Output","H3K4me3","Signatures","{name}.Q5DD.csv"),
        out3=join(workpath,cfTool_dir,"Samples","H3K4me3","{name}.Q5DD.rdata"),
    params:
        rname='cfChiP',
        rver="R/4.1.0",
        toolkit = config['references'][genome]['cfChIP_TOOLS_SRC'],
        tmpfile = lambda w: join(workpath,cfTool_subdir2, w.name + ".Q5DD.tagAlign"),
        tag=lambda w: temp(join(workpath,bam_dir, w.name+".Q5DD.tagAlign"))
    container:
        config['images']['cfchip']
    shell: """
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
        expand(join(workpath,cfTool_dir,"Output","H3K4me3","Signatures","{name}.Q5DD.csv"),name=chip)
    output:
        txt=join(workpath,"QC","H3K4me3_cfChIP_signature.txt"),
        pdf=join(workpath,"QC","H3K4me3_cfChIP_signature.pdf")
    params:
        rname="cfChIP2",
        script=join(workpath,"workflow","scripts","cfChIP_signatures.R"),
        infolder=join(workpath,cfTool_dir,"Output","H3K4me3","Signatures"),
    container:
        config['images']['cfchip']
    shell: """
    Rscript -e "source('{params.script}'); mergeSignatures( '{params.infolder}', '{output.txt}' )";
    Rscript -e "source('{params.script}'); plotSignatures( '{output.txt}', '{output.pdf}' )";
    """

rule promoterTable1:
    input:
        expand(join(workpath,uropa_dir,'{PeakTool}','{name}_{PeakTool}_uropa_protTSS_allhits.txt'),PeakTool=PeakTools,name=chips),
    output:
        txt=join(workpath,uropa_dir,downstream_dir,'{PeakTool}_promoter_overlap_summaryTable.txt'),
    params:
        rname="promoter1",
        script=join(workpath,"workflow","scripts","promoterAnnotation_by_Gene.R"),
        infolder= lambda w: join(workpath,uropa_dir, w.PeakTool)
    container:
        config['images']['cfchip']
    shell: """
    Rscript -e "source('{params.script}'); peakcallVersion('{params.infolder}','{output.txt}')";
    """

rule promoterTable2:
    input:
        expand(join(workpath,uropa_dir,diffbind_dir,'{name}_{PeakTool}_uropa_protTSS_allhits.txt'),PeakTool='DiffbindDeseq2',name=contrasts),
    output:
        txt=join(workpath,uropa_dir,downstream_dir,'{PeakTool}_promoter_overlap_summary_table.txt'),
    params:
        rname="promoter2",
        script1=join(workpath,"workflow","scripts","promoterAnnotation_by_Gene.R"),
        script2=join(workpath,"workflow","scripts","significantPathways.R"),
        infolder= workpath,
        gtf = config['references'][genome]['GTFFILE'],
    container:
        config['images']['cfchip']
    shell: """
    Rscript -e "source('{params.script1}'); diffbindVersion('{params.infolder}','{output.txt}')";
    Rscript -e "source('{params.script2}'); promoterAnnotationWrapper('{output.txt}','{params.gtf}','KEGG')";
    Rscript -e "source('{params.script2}'); promoterAnnotationWrapper('{output.txt}','{params.gtf}','Reactome')";
    """
