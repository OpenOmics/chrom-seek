# cell-free ChIP-seq
# ~~~~
# rules: picard_dedup, cfChIPtool, cfChIPcompile


# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
blocks                          = config['project']['blocks']
groupdata                       = config['project']['groups']
chips                           = config['project']['peaks']['chips']

# Directory end points
bam_dir                         = join(workpath, "bam")
cfTool_dir                      = join(workpath, "cfChIPtool")
cfTool_subdir2                  = join(cfTool_dir, "BED", "H3K4me3")
qc_dir                          = join(workpath, "QC")
diffbind_dir                    = join(workpath, "DiffBind")


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


rule promoterTable1_macsN:
    input:
        expand(join(uropa_dir, "macsNarrow", "{name}_uropa_protTSS_allhits.txt"), name=chips)
    output:
        txt                     = join(uropa_dir, "promoterTable1", "macsNarrow_promoter_overlap_summaryTable.txt")
    params:
        rname                   = "promoterTable1",
        script                  = join(bin_path, "promoterAnnotation_by_Gene.R"),
        infolder                = join(uropa_dir, "macsNarrow")
    container:
        config['images']['cfchip']
    shell: 
        """
        Rscript -e "source('{params.script}'); peakcallVersion('{params.infolder}','{output.txt}')";
        """


rule promoterTable2:
    input:
        full_list               = join(
                                    diffbind_dir, 
                                    "{group1}_vs_{group2}-{PeakTool}", 
                                    "{group1}_vs_{group2}-{PeakTool}_Diffbind_fullList.bed"
                                  ),
        peak_annos              = [join(uropa_diffbind_dir, 
                                        "{group1}_vs_{group2}-{PeakTool}-DeSeq2",
                                        "{group1}_vs_{group2}_{PeakTool}_DeSeq2_protTSS_uropa_allhits.txt") 
                                    for pk_type in peak_types],
    output:
        join(
            uropa_diffbind_dir,
            "{group1}_vs_{group2}-{PeakTool}-DeSeq2",
            "{group1}_vs_{group2}-{PeakTool}_Diffbind_DeSeq2_promoter_overlap_summaryTable.txt"
        ),
    params:
        rname                   = "promoterTable2",
        script1                 = join(bin_path, "promoterAnnotation_by_Gene.R"),
        script2                 = join(bin_path, "significantPathways.R"),
        infolder                = workpath,
        gtf                     = config['references'][genome]['GTFFILE'],
    container:
        config['images']['cfchip']
    shell: 
        """
        Rscript -e "source('{params.script1}'); diffbindVersion('{params.infolder}', '{output}')";
        Rscript -e "source('{params.script2}'); promoterAnnotationWrapper('{output}', '{params.gtf}', 'KEGG')";
        Rscript -e "source('{params.script2}'); promoterAnnotationWrapper('{output}', '{params.gtf}', 'Reactome')";
        """


rule diffbindQC_macsN:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_TMMcounts.csv"),
        umap                        = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_DiffBindQC_UMAP.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_macsN",
        peak_tool                   = "macsNarrow",
        peak_type                   = "narrow",
        rscript                     = join(bin_path, "DiffBind_v2_QC_cfchip.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-macsNarrow"),
        pythonscript                = join(bin_path, "prep_diffbindQC.py"),
    container:
       config['images']['diffbind']
    shell:
        """
        python {params.pythonscript} \\
            -s {input.sample_bams} \\
            -c {input.control_bams} \\
            -p {input.samples_peaks} \\
            -t {params.peak_type} \\
            -o {output.csvfile}
        cp {params.rscript} {params.outdir}
        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", umapfile="{output.umap}", 
            counts_bed="{output.countsbed}", counts_csv="{output.countscsv}",
            peakcaller="{params.peak_tool}"))'
        """


rule diffbindQC_macsB:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(macsB_dir, "{name}", "{name}_peaks.broadPeak"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_TMMcounts.csv"),
        umap                        = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_DiffBindQC_UMAP.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_macsB",
        peak_tool                   = "macsBroad",
        peak_type                   = "narrow",
        rscript                     = join(bin_path, "DiffBind_v2_QC_cfchip.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-macsBroad"),
        pythonscript                = join(bin_path, "prep_diffbindQC.py"),
    container:
       config['images']['diffbind']
    shell:
        """
        python {params.pythonscript} \\
            -s {input.sample_bams} \\
            -c {input.control_bams} \\
            -p {input.samples_peaks} \\
            -t {params.peak_type} \\
            -o {output.csvfile}
        cp {params.rscript} {params.outdir}
        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", umapfile="{output.umap}", 
            counts_bed="{output.countsbed}", counts_csv="{output.countscsv}",
            peakcaller="{params.peak_tool}"))'
        """


rule diffbindQC_genrich:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(genrich_dir, "{name}", "{name}.narrowPeak"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_TMMcounts.csv"),
        umap                        = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_DiffBindQC_UMAP.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_genrich",
        peak_tool                   = "Genrich",
        peak_type                   = "narrow",
        rscript                     = join(bin_path, "DiffBind_v2_QC_cfchip.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-Genrich"),
        pythonscript                = join(bin_path, "prep_diffbindQC.py"),
    container:
       config['images']['diffbind']
    shell:
        """
        python {params.pythonscript} \\
            -s {input.sample_bams} \\
            -c {input.control_bams} \\
            -p {input.samples_peaks} \\
            -t {params.peak_type} \\
            -o {output.csvfile}
        cp {params.rscript} {params.outdir}
        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", umapfile="{output.umap}", 
            counts_bed="{output.countsbed}", counts_csv="{output.countscsv}",
            peakcaller="{params.peak_tool}"))'
        """


rule diffbindQC_SEACR:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(seacr_dir, "{name}.stringent.bed"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_TMMcounts.csv"),
        umap                        = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_DiffBindQC_UMAP.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_SEACR",
        peak_tool                   = "SEACR",
        peak_type                   = "raw",
        rscript                     = join(bin_path, "DiffBind_v2_QC_cfchip.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-SEACR"),
        pythonscript                = join(bin_path, "prep_diffbindQC.py"),
    container:
       config['images']['diffbind']
    shell:
        """
        python {params.pythonscript} \\
            -s {input.sample_bams} \\
            -c {input.control_bams} \\
            -p {input.samples_peaks} \\
            -t {params.peak_type} \\
            -o {output.csvfile}
        cp {params.rscript} {params.outdir}
        cd {params.outdir}
        Rscript -e 'rmarkdown::render("{params.rscript}", output_file="{output.html}",
            params=list(csvfile="{output.csvfile}", umapfile="{output.umap}", 
            counts_bed="{output.countsbed}", counts_csv="{output.countscsv}",
            peakcaller="{params.peak_tool}"))'
        """