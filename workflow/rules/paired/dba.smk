rule diffbindQC_macsN:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC.html"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBind_prep.csv"),
        countspca                   = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_PCA_peaks_and_counts.txt"),
        peakspca                    = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_PCA_peaks_only.txt"),
        countstmm                   = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_TMMcounts.csv"),
        countsrpkm                  = join(diffbind_qc_dir, "AllSamples-macsNarrow", "AllSamples-macsNarrow_DiffBindQC_RPKMcounts.csv"),
    params:
        rname                       = "diffbindQC_macsN",
        peak_tool                   = "macsNarrow",
        peak_type                   = "narrow",
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
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
            peakcaller="{params.peak_tool}", pca_counts_file="{output.countspca}", pca_peaks_file="{output.peakspca}",
            tmm_counts_file="{output.countstmm}", rpkm_counts_file="{output.countsrpkm}"))'
        """


rule diffbindQC_macsB:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(macsB_dir, "{name}", "{name}_peaks.broadPeak"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC.html"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBind_prep.csv"),
        countspca                   = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_PCA_peaks_and_counts.txt"),
        peakspca                    = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_PCA_peaks_only.txt"),
        countstmm                   = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_TMMcounts.csv"),
        countsrpkm                  = join(diffbind_qc_dir, "AllSamples-macsBroad", "AllSamples-macsBroad_DiffBindQC_RPKMcounts.csv"),
    params:
        rname                       = "diffbindQC_macsB",
        peak_tool                   = "macsBroad",
        peak_type                   = "narrow",
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
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
            peakcaller="{params.peak_tool}", pca_counts_file="{output.countspca}", pca_peaks_file="{output.peakspca}",
            tmm_counts_file="{output.countstmm}", rpkm_counts_file="{output.countsrpkm}"))'
        """


rule diffbindQC_genrich:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(genrich_dir, "{name}", "{name}.narrowPeak"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC.html"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBind_prep.csv"),
        countspca                   = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_PCA_peaks_and_counts.txt"),
        peakspca                    = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_PCA_peaks_only.txt"),
        countstmm                   = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_TMMcounts.csv"),
        countsrpkm                  = join(diffbind_qc_dir, "AllSamples-Genrich", "AllSamples-Genrich_DiffBindQC_RPKMcounts.csv"),
    params:
        rname                       = "diffbindQC_genrich",
        peak_tool                   = "Genrich",
        peak_type                   = "narrow",
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
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
            peakcaller="{params.peak_tool}", pca_counts_file="{output.countspca}", pca_peaks_file="{output.peakspca}",
            tmm_counts_file="{output.countstmm}", rpkm_counts_file="{output.countsrpkm}"))'
        """


rule diffbindQC_SEACR:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(seacr_dir, "{name}.stringent.bed"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC.html"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBind_prep.csv"),
        countspca                   = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_PCA_peaks_and_counts.txt"),
        peakspca                    = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_PCA_peaks_only.txt"),
        countstmm                   = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_TMMcounts.csv"),
        countsrpkm                  = join(diffbind_qc_dir, "AllSamples-SEACR", "AllSamples-SEACR_DiffBindQC_RPKMcounts.csv"),
    params:
        rname                       = "diffbindQC_SEACR",
        peak_tool                   = "SEACR",
        peak_type                   = "raw",
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
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
            peakcaller="{params.peak_tool}", pca_counts_file="{output.countspca}", pca_peaks_file="{output.peakspca}",
            tmm_counts_file="{output.countstmm}", rpkm_counts_file="{output.countsrpkm}"))'
        """