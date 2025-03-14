# Global workflow variables
searc_dir                       = join(workpath, "searc")
bg_dir                          = join(workpath, "bedgraph")
tmpdir                          = config['options']['tmp_dir']
genome                          = config['options']['genome']
seacr_dir                       = join(workpath, "SEACR")


rule bam2bedgraph:
    input: join(bam_dir, "{name}.Q5DD.bam")
    output: temp(join(bg_dir, "{name}.bedgraph"))
    container: config['images']['seacr']
    params:
        rname                           = 'bam2bedgraph',
        tmpdir                          = tmpdir,
        reference                       = config['references'][genome]['REFLEN'],
    threads:
        int(cluster['bam2bedgraph'].get('threads', cluster['__default__']['threads']))
    shell:
        """
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT
        export TMPDIR="${{tmp}}" # used by sort

        # bam to bedgraph
        samtools sort -@{threads} -n {input} > ${{tmp}}/name_sorted_sample.bam
        bedtools bamtobed -bedpe -i ${{tmp}}/name_sorted_sample.bam > ${{tmp}}/sample.bed
        awk '$1==$4 && $6-$2 < 1000 {{print $0}}' ${{tmp}}/sample.bed > ${{tmp}}/sample.clean.bed
        cut -f 1,2,6 ${{tmp}}/sample.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${{tmp}}/sample.fragments.bed
        bedtools genomecov -bg -i ${{tmp}}/sample.fragments.bed -g {params.reference} > {output}
        """



rule SEACR:
    input:
        exp                             = join(bg_dir, "{name}.bedgraph"),
        control                         = lambda w: join(bg_dir, f"{chip2input[w.name]}.bedgraph")
                                            if chip2input[w.name] else [],
    output:
        peaks                           = join(seacr_dir, "{name}.stringent.bed")
    params:
        rname                           = 'SEARC',
        out_dir                         = seacr_dir,
        control_flag                    = lambda w, input: input.control if input.control else "0.01",
    container:
        config['images']['seacr']
    shell:
        """
        cd {params.out_dir}
        SEACR_1.3.sh {input.exp} {params.control_flag} non stringent {wildcards.name}
        """


rule diffbindQC_SEARC:
    input:
        sample_bams                 = expand(join(bam_dir, "{name}.Q5DD.bam"), name=chips),
        control_bams                = [join(bam_dir, f"{chip2input[sample_name]}.Q5DD.bam") for sample_name in chips],
        samples_peaks               = expand(join(seacr_dir, "{name}.stringent.bed"), name=chips)
    output:
        html                        = join(diffbind_qc_dir, "AllSamples-SEARC", "AllSamples-SEARC_DiffBindQC.html"),
        countsbed                   = join(diffbind_qc_dir, "AllSamples-SEARC", "AllSamples-SEARC_DiffBindQC_TMMcounts.bed"),
        countscsv                   = join(diffbind_qc_dir, "AllSamples-SEARC", "AllSamples-SEARC_DiffBindQC_TMMcounts.csv"),
        umap                        = join(diffbind_qc_dir, "AllSamples-SEARC", "AllSamples-SEARC_DiffBindQC_DiffBindQC_UMAP.csv"),
        csvfile                     = join(diffbind_qc_dir, "AllSamples-SEARC", "AllSamples-SEARC_DiffBind_prep.csv"),
    params:
        rname                       = "diffbindQC_SEARC",
        peak_tool                   = "SEARC",
        peak_type                   = "raw",
        rscript                     = join(bin_path, "DiffBind_v2_QC.Rmd"),
        outdir                      = join(diffbind_qc_dir, "AllSamples-SEARC"),
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
