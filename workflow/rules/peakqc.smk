# Quality control metrics and reports
# ~~~~
# Generally applicable quality control rules
from scripts.common import allocated
from textwrap import dedent
from os.path import join


# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
paired_end                      = False if config['project']['nends'] == 1 else True
samples                         = config['samples']
ends                            = [1] if not paired_end else [1, 2]
assay                           = config['options']['assay']
PeakTools                       = get_peaktools(assay)
chips                           = config['project']['peaks']['chips']

# ~~ directories
peakqc_dir                      = join(workpath, "PeakQC")
tmpdir                          = config['options']['tmp_dir']



rule FRiP_macsN:
    input:
        peaks                   = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        tbl                     = join(peakqc_dir, "FRiP", "macsNarrow", "macsNarrow.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_macsN",
        tblscript               = join(bin_path, "frip.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir,
        this_config             = join(workpath, 'config.json')
    container: 
        config['images']['python']
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.tblscript} \\
            -p {input.peaks} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {output.tbl} \\
            -x 16
        """


rule FRiP_Genrich:
    input:
        peaks                   = join(genrich_dir, "{name}", "{name}.narrowPeak"),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        tbl                     = join(peakqc_dir, "FRiP", "Genrich", "Genrich.{name}.Q5DD.FRiP_table.txt")
    params:
        rname                   = "FRiP_Genrich",
        tblscript               = join(bin_path, "frip.py"),
        plotscript              = join(bin_path, "FRiP_plot.R"),
        this_config             = join(workpath, 'config.json'),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir,
    container: 
        config['images']['python']
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.tblscript} \\
            -p {input.peaks} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {output.tbl} \\
            -x 16
        """


rule FRiP_macsB:
    input:
        peaks                   = expand(join(macsB_dir, "{name}", "{name}_peaks.broadPeak"), name=chips),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        tbl                     = join(peakqc_dir, "FRiP", "macsBroad", "macsBroad.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_macsB",
        tblscript               = join(bin_path, "frip.py"),
        plotscript              = join(bin_path, "FRiP_plot.R"),
        this_config             = join(workpath, 'config.json'),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir,
    container: 
        config['images']['python']
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.tblscript} \\
            -p {input.peaks} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {output.tbl} \\
            -x 16
        """


rule FRiP_SEACR:
    input:
        peaks                   = join(seacr_dir, "{name}", "{name}.stringent.bed"),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        tbl                     = join(peakqc_dir, "FRiP", "SEACR", "SEACR.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_SEACR",
        tblscript               = join(bin_path, "frip.py"),
        plotscript              = join(bin_path, "FRiP_plot.R"),
        this_config             = join(workpath, 'config.json'),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir,
    container: 
        config['images']['python']
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.tblscript} \\
            -p {input.peaks} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {output.tbl} \\
            -x 16
        """


rule frip_summary:
    input:
        tbl                     = expand(join(peakqc_dir, "FRiP", "{pktool}", "{pktool}.{name}.Q5DD.FRiP_table.txt"), pktool=PeakTools, name=samples)
    output:
        heatmap                 = join(peakqc_dir, "FRiP", "summary.Q5DD.FRiP_heatmap.pdf"),
        bar_plot                = join(peakqc_dir, "FRiP", "summary.Q5DD.FRiP_barplot.pdf"),
        # scatter_plot            = join(peakqc_dir, "FRiP", "summary.Q5DD.FRiP_scatter.pdf")
    params:
        rname                   = "frip_summary",
        plotscript              = join(bin_path, "frip_summary_plot.py"),
        this_config             = join(workpath, 'config.json')
    container: 
        config['images']['python']
    shell:
        dedent("""
        {params.plotscript} \\
            --hm {output.heatmap} \\
            -t {input.tbl} \\
            -b {output.bar_plot} \\
            -c {params.this_config} \\
            -v
        """)


rule jaccard_genrich:
    input:
        expand(join(genrich_dir, "{SID}", "{SID}.narrowPeak"), SID=chips),
    output:
        table                   = join(peakqc_dir, "jaccard", 'Genrich_jaccard.txt'),
        pcaplot                 = join(peakqc_dir, "jaccard", 'Genrich_jaccard_pca.pdf'),
        pcatab                  = join(peakqc_dir, "jaccard", 'Genrich_jaccard_pca.tsv'),
        heatmap                 = join(peakqc_dir, "jaccard", 'Genrich_jaccard_heatmap.pdf'),
        heatmaptab              = join(peakqc_dir, "jaccard", 'Genrich_jaccard_heatmap.tsv'),
    params:
        rname                   = "jaccard_genrich",
        outroot                 = join(peakqc_dir, "jaccard"),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \\
            -i "{input}" \\
            --caller Genrich \\
            --pcatab {output.pcatab} \\
            --outtable {output.table} \\
            --pcaplot {output.pcaplot} \\
            --outheatmap {output.heatmap} \\
            --tabheatmap {output.heatmaptab} \\
            -g {params.genome}
        """)


rule jaccard_macsbroad:
    input:
        expand(join(macsB_dir, "{SID}", "{SID}_peaks.broadPeak"), SID=chips),
    output:
        table                   = join(peakqc_dir, "jaccard", 'macsBroad_jaccard.txt'),
        pcaplot                 = join(peakqc_dir, "jaccard", 'macsBroad_jaccard_pca.pdf'),
        pcatab                  = join(peakqc_dir, "jaccard", 'macsBroad_jaccard_pca.tsv'),
        heatmap                 = join(peakqc_dir, "jaccard", 'macsBroad_jaccard_heatmap.pdf'),
        heatmaptab              = join(peakqc_dir, "jaccard", 'macsBroad_jaccard_heatmap.tsv'),
    params:
        rname                   = "jaccard_macsbroad",
        outroot                 = join(peakqc_dir, "jaccard"),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \\
            -i "{input}" \\
            --caller macsBroad \\
            --pcatab {output.pcatab} \\
            --outtable {output.table} \\
            --pcaplot {output.pcaplot} \\
            --outheatmap {output.heatmap} \\
            --tabheatmap {output.heatmaptab} \\
            -g {params.genome}
        """)


rule jaccard_macsnarrow:
    input:
        expand(join(macsN_dir, "{SID}", "{SID}_peaks.narrowPeak"), SID=chips),
    output:
        table                   = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard.txt'),
        pcaplot                 = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard_pca.pdf'),
        pcatab                  = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard_pca.tsv'),
        heatmap                 = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard_heatmap.pdf'),
        heatmaptab              = join(peakqc_dir, "jaccard", 'macsNarrow_jaccard_heatmap.tsv'),
    params:
        rname                   = "jaccard_macsnarrow",
        outroot                 = join(peakqc_dir, "jaccard"),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \\
            -i "{input}" \\
            --caller macsNarrow \\
            --pcatab {output.pcatab} \\
            --outtable {output.table} \\
            --pcaplot {output.pcaplot} \\
            --outheatmap {output.heatmap} \\
            --tabheatmap {output.heatmaptab} \\
            -g {params.genome}
        """)


rule jaccard_seacr:
    input:
        expand(join(seacr_dir, "{SID}", "{SID}.stringent.bed"), SID=chips),
    output:
        table                   = join(peakqc_dir, "jaccard", 'SEACR_jaccard.txt'),
        pcaplot                 = join(peakqc_dir, "jaccard", 'SEACR_jaccard_pca.pdf'),
        pcatab                  = join(peakqc_dir, "jaccard", 'SEACR_jaccard_pca.tsv'),
        heatmap                 = join(peakqc_dir, "jaccard", 'SEACR_jaccard_heatmap.pdf'),
        heatmaptab              = join(peakqc_dir, "jaccard", 'SEACR_jaccard_heatmap.tsv')
    params:
        rname                   = "jaccard_seacr",
        outroot                 = join(peakqc_dir, "jaccard"),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN'],
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \\
            -i "{input}" \\
            --caller SEARC \\
            --pcatab {output.pcatab} \\
            --outtable {output.table} \\
            --pcaplot {output.pcaplot} \\
            --outheatmap {output.heatmap} \\
            --tabheatmap {output.heatmaptab} \\
            -g {params.genome}
        """)


rule jaccard_summary:
    input:
        pca_cords               = expand(join(peakqc_dir, "jaccard", '{pkcaller}_jaccard_pca.tsv'), pkcaller=PeakTools),
        hm_cords                = expand(join(peakqc_dir, "jaccard", '{pkcaller}_jaccard_heatmap.tsv'), pkcaller=PeakTools),
    output:
        pcaplot                 = join(peakqc_dir, "jaccard", 'jaccard_summary_pca.pdf'),
        hmplot                 = join(peakqc_dir, "jaccard", 'jaccard_summary_heatmap.pdf'),
    params:
        rname                   = "jaccard_summary",
        genome                  = config['references'][genome]['REFLEN'],
        script                  = join(bin_path, "jaccard_summary.py"),
        tmpdir                  = tmpdir
    container: 
        config['images']['python']
    shell: 
        dedent("""
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        export TMPDIR="${{tmp}}"
        trap 'rm -rf "${{tmp}}"' EXIT
        python {params.script} \
            --pca {input.pca_cords} \
            --hm {input.hm_cords}
        """)
