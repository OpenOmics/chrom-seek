# Quality control rules
# ~~~~
# Common quality-control rules: preseq, NRF, rawfastqc,
#   fastqc, fastq_screen, multiQC
from os.path import join
from scripts.common import get_bam_ext, get_fqscreen_outputs


# ~~ workflow configuration
workpath                        = config['project']['workpath']
bin_path                        = config['project']['binpath']
genome                          = config['options']['genome']
paired_end                      = False if config['project']['nends'] == 1 else True
samples                         = config['samples']
ends                            = [1] if not paired_end else [1, 2]
 

# ~~ directories
qc_dir                          = join(workpath, "QC")
kraken_dir                      = join(workpath, 'kraken')
deeptools_dir                   = join(workpath, 'deeptools')
peakqc_dir                      = join(workpath, "PeakQC")
extra_fingerprint_dir           = join(deeptools_dir, 'sorted_fingerprint')


rule preseq:
    """
    Quality step to estimate library complexity. Low library complexity may indicate
    an issue with library preparation or sample storage (FFPE samples) where very
    little input RNA was over-amplified or the sample may be highly degraded.
    @Input:
        Sorted, duplicate marked genomic BAM file (scatter)
    @Output:
        Logfile containing library complexity information
    """
    input:
        bam                     = join(bam_dir, "{name}.sorted.bam"),
    output:
        ccurve                  = join(qc_dir, "{name}.ccurve"),
    params:
        rname                   = "preseq",
        preseqver               = config['tools']['PRESEQVER'],
    shell: 
        """
        module load {params.preseqver};
        preseq c_curve \\
            -B \\
            -o {output.ccurve} \\
            {input.bam}            
        """


rule NRF:
    """
    Quality step computes the expected yield for theoretical 
    larger experiments and the associated confidence intervals
    @Input:
        Sorted BAM file (scatter)
    @Output:
        Output is a text file with four columns. The total number 
        of reads, average expected number of distinct reads, and
        the lower and upper limits of the confidence interval.
        And pyhton code produces NRF = distinct_reads/tot_reads,
        PBC1 = one_pair/distinct_reads, and PBC2 = one_pair/two_pair.
    """
    input:
        bam                     = join(bam_dir, "{name}.sorted.bam"),
    output:
        preseq                  = join(qc_dir, "{name}.preseq.dat"),
        preseqlog               = join(qc_dir, "{name}.preseq.log"),
        nrf                     = temp(join(qc_dir, "{name}.nrf")),
    params:
        rname                   = 'NRF',
        samtoolsver             = config['tools']['SAMTOOLSVER'],
        rver                    = config['tools']['RVER'],
        preseqver               = config['tools']['PRESEQVER'],
        nrfscript               = join(bin_path, "atac_nrf.py"),
    threads: 16
    shell: 
        """
        module load {params.preseqver};
        preseq lc_extrap \\
            -P \\
            -B \\
            -D \\
            -o {output.preseq} \\
            {input.bam} \\
            -seed 12345 \\
            -v \\
            -l 100000000000 \\
        2> {output.preseqlog}
        python {params.nrfscript} \\
            {output.preseqlog} \\
        > {output.nrf}
        """


rule rawfastqc:
    """
    Quality-control step to assess sequencing quality of the raw data prior removing
    adapter sequences. FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Raw FastQ files (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """
    input:
        expand(join(workpath, "{name}.R{rn}.fastq.gz"), name=samples, rn=list(map(str, ends)))
    output:
        expand(join(qc_dir, "rawfastQC", "{name}.R{rn}_fastqc.html"), name=samples, rn=ends),
    params:
        rname                   = 'rawfastqc',
        outdir                  = join(qc_dir, "rawfastQC"),
        tmpdir                  = tmpdir,
    envmodules: 
        config['tools']['FASTQCVER']
    threads:
        int(allocated("threads", "rawfastqc", cluster))
    shell: 
        """
        # fastqc storage on lscratch b/c nfs bug
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        fastqc \\
            {input} \\
            -t {threads} \\
            -o "${{tmp}}"
            
        find "${{tmp}}" \\
            -type f \\
            \\( -name '*.html' -o -name '*.zip' \\) \\
            -exec cp {{}} {params.outdir} \\; 
        """


rule fastqc:
    """
    Quality-control step to assess sequencing quality of the raw data after removing
    adapter sequences. This step is run after trim_pe rule. FastQC is run after adapter
    trimming to evalute if the adapter sequences were properly removed.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Trimmed FastQC reports and zip file containing data quality information
    """
    input:
        expand(join(trim_dir, "{name}.R{rn}.trim.fastq.gz"), name=samples, rn=ends)
    output:
        expand(join(qc_dir, 'fastQC', "{name}.R{rn}.trim_fastqc.html"), name=samples, rn=ends),
    params:
        rname                   = 'fastqc',
        outdir                  = join(qc_dir, "fastQC"),
        tmpdir                  = tmpdir,
    envmodules: 
        config['tools']['FASTQCVER']
    threads:
        int(allocated("threads", "fastqc", cluster))
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        # Running fastqc with local
        # disk or a tmpdir, fastqc
        # has been observed to lock
        # up gpfs filesystems, adding
        # this on request by HPC staff
        fastqc \\
            {input} \\
            -t {threads} \\
            -o "${{tmp}}"
        
        # Copy output files from tmpdir
        # to output directory
        find "${{tmp}}" \\
            -type f \\
            \\( -name '*.html' -o -name '*.zip' \\) \\
            -exec cp {{}} {params.outdir} \\;
        """


rule fastq_screen:
    """
    Quality-control step to screen for different sources of contamination.
    FastQ Screen compares your sequencing data to a set of different reference
    genomes to determine if there is contamination. It allows a user to see if
    the composition of your library matches what you expect.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        FastQ Screen report and logfiles
    """
    input:
        expand(join(trim_dir, "{name}.R{rn}.trim.fastq.gz"), name=samples, rn=ends)
    output:
        get_fqscreen_outputs(paired_end, samples, qc_dir)
    params:
        rname                   = 'fqscreen',
        outdir                  = join(qc_dir, "FQscreen"),
        outdir2                 = join(qc_dir, "FQscreen2"),
        fastq_screen            = config['bin']['FASTQ_SCREEN'],
        fastq_screen_config1    = config['shared_resources']['FASTQ_SCREEN_CONFIG_P1'],
        fastq_screen_config2    = config['shared_resources']['FASTQ_SCREEN_CONFIG_P2'],
    envmodules:
        config['tools']['BOWTIE2VER'],
        config['tools']['PERLVER'],
    threads: 
        int(allocated("threads", "fastq_screen", cluster))
    shell: 
        """
        # First pass of contamination screening
        {params.fastq_screen} \\
            --conf {params.fastq_screen_config1} \\
            --outdir {params.outdir} \\
            --threads {threads} \\
            --subset 1000000 \\
            --aligner bowtie2 \\
            --force \\
            {input}
        # Second pass of contamination screening
        {params.fastq_screen} \\
            --conf {params.fastq_screen_config2} \\
            --outdir {params.outdir2} \\
            --threads {threads} \\
            --subset 1000000 \\
            --aligner bowtie2 \\
            --force \\
            {input}
        """


rule kraken:
    """
    Quality-control step to assess for potential sources of microbial contamination.
    If there are high levels of microbial contamination, Kraken will provide an
    estimation of the taxonomic composition. Kraken is used in conjunction with
    Krona to produce an interactive reports.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Kraken logfile and interative krona report
    """
    input:
        fq1                     = join(trim_dir, "{name}.R1.trim.fastq.gz"),
        fq2                     = provided(join(trim_dir,"{name}.R2.trim.fastq.gz"), paired_end)
    output:
        krakenout               = join(kraken_dir, "{name}.trim.kraken_bacteria.out.txt"),
        krakentaxa              = join(kraken_dir, "{name}.trim.kraken_bacteria.taxa.txt"),
        kronahtml               = join(kraken_dir, "{name}.trim.kraken_bacteria.krona.html"),
    params:
        rname                   = 'kraken',
        outdir                  = kraken_dir,
        bacdb                   = config['shared_resources']['KRAKENBACDB'],
        tmpdir                  = tmpdir,
        paired_end              = paired_end
    threads: 
        int(allocated("threads", "kraken_pe", cluster)),
    envmodules:
        config['tools']['KRAKENVER'],
        config['tools']['KRONATOOLSVER'],
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT
        
        # Copy kraken2 db to /lscratch or temp 
        # location to reduce filesystem strain
        cp -rv {params.bacdb} ${{tmp}}/;
        kdb_base=$(basename {params.bacdb})
        if [ '{params.paired_end}' == True ]; then
            kraken2 --db ${{tmp}}/${{kdb_base}} \\
                --threads {threads} --report {output.krakentaxa} \\
                --output {output.krakenout} \\
                --gzip-compressed \\
                --paired {input.fq1} {input.fq2}
        else
            kraken2 --db ${{tmp}}/${{kdb_base}} \\
                --threads {threads} --report {output.krakentaxa} \\
                --output {output.krakenout} \\
                --gzip-compressed \\
                {input.fq1}
        fi
        
        # Generate Krona Report
        cut -f2,3 {output.krakenout} | \\
            ktImportTaxonomy - -o {output.kronahtml}
        """


rule multiqc:
    """
    Reporting step to aggregate sample statistics and quality-control information
    across all samples. This will be one of the last steps of the pipeline. The inputs
    listed here are to ensure that this step runs last. During runtime, MultiQC will
    recurively crawl through the working directory and parse files that it supports.
    @Input:
        List of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report and a QC metadata table
    """
    input: 
        expand(join(qc_dir, "FQscreen", "{name}.R1.trim_screen.txt"), name=samples),
        expand(join(qc_dir, "FQscreen2", "{name}.R1.trim_screen.txt"), name=samples),
        expand(join(qc_dir, "{name}.ccurve"), name=samples),
        expand(join(qc_dir, "rawfastQC", "{name}.R1_fastqc.html"), name=samples),
        expand(join(qc_dir, "fastQC", "{name}.R1.trim_fastqc.html"), name=samples),
        expand(join(kraken_dir, "{name}.trim.kraken_bacteria.krona.html"), name=samples),
        expand(join(bam_dir, "{name}.Q5DD.bam.flagstat"), name=samples),
        expand(join(bam_dir, "{name}.Q5.bam.flagstat"), name=samples),
        join(deeptools_dir, "spearman_readcounts.Q5DD.tab"),
        join(deeptools_dir, "fingerprint.raw.Q5DD.tab"),
	join(deeptools_dir,"TSS_profile.Q5DD.tab")
    output:
        join(workpath, "multiqc_report.html")
    params:
        rname                   = "multiqc",
        multiqc                 = config['tools']['MULTIQCVER'],
	    qcconfig                = join(workpath, config['shared_resources']['MULTIQC_CONFIG']),
	    excludedir              = join(workpath, extra_fingerprint_dir),
    shell: 
        """
        module load {params.multiqc}
        multiqc \\
            -f \\
            -c {params.qcconfig} \\
            --interactive \\
            -e cutadapt \\
            --ignore {params.excludedir} \\
            -d """ + workpath + """
        """


rule insert_size:
    """
    Quality step calculates number of reads per insert size.
    @Input:
        Sorted only bam file, and also bam files that were sorted, 
        filtered by mapQ a value, and deduplicated (extensions: sorted and Q5DD),
        for all samples.
    @Output:
        Number of reads per insert size and their histogram
    """
    input:
        bam                     = lambda w : join(bam_dir, w.name + "." + w.ext + "." + get_bam_ext(w.ext, paired_end))
    output:
        txt                     = join(qc_dir, "{name}.{ext}.insert_size_metrics.txt"),
        pdf                     = join(qc_dir, "{name}.{ext}.insert_size_histogram.pdf"),
    params:
        rname                   = "insert_size",
        picardver               = config['tools']['PICARDVER'],
        rver                    = config['tools']['RVER'],
        javaram                 = '16g',
    shell: 
        """
        module load {params.picardver} {params.rver};
        java -Xmx{params.javaram} -jar ${{PICARDJARPATH}}/picard.jar CollectInsertSizeMetrics \\
            -I {input.bam} \\
            -O {output.txt} \\
            -H {output.pdf}
        """

rule deeptools_QC:
    input:
        [ join(bw_dir, name + ".Q5DD.RPGC.bw") for name in samples ] 
    output:
        heatmap                 = join(deeptools_dir, "spearman_heatmap.Q5DD.pdf"),
        pca                     = join(deeptools_dir, "pca.Q5DD.pdf"),
        npz                     = temp(join(deeptools_dir, "Q5DD.npz")),
	mqc                     = join(deeptools_dir, "spearman_readcounts.Q5DD.tab"),
	png                     = join(deeptools_dir, "spearman_heatmap.Q5DD_mqc.png")
    params:
        rname                   = "deeptools_QC",
        parent_dir              = deeptools_dir,
        deeptoolsver            = config['tools']['DEEPTOOLSVER'],
        # this should be the sample names to match the bigwigs in the same order
        labels                  = samples
    threads: 24
    shell: 
        """    
        module load {params.deeptoolsver}
        if [ ! -d "{params.parent_dir}" ]; then mkdir "{params.parent_dir}"; fi
        multiBigwigSummary bins -b {input} -p {threads} -l {params.labels} -out {output.npz}
        plotCorrelation -in {output.npz} -o {output.heatmap} -c 'spearman' -p 'heatmap' \\
               --skipZeros --removeOutliers --outFileCorMatrix {output.mqc}
        plotCorrelation -in {output.npz} -o {output.png} -c 'spearman' -p 'heatmap' --skipZeros --removeOutliers
        plotPCA -in {output.npz} -o {output.pca}
        """

rule deeptools_fingerprint:
    input:
        [ join(bam_dir, name + ".Q5DD.bam") for name in samples ] 
    output:
        image=join(deeptools_dir, "fingerprint.Q5DD.pdf"),
        raw=temp(join(deeptools_dir, "fingerprint.raw.Q5DD.tab")),
        metrics=join(deeptools_dir, "fingerprint.metrics.Q5DD.tsv"),
    params:
        rname                   = "deeptools_fingerprint",
        parent_dir              = deeptools_dir,
        deeptoolsver            = config['tools']['DEEPTOOLSVER'],
        # this should be the sample names to match the bigwigs in the same order
        labels                  = samples
    threads: int(allocated("threads", "deeptools_fingerprint", cluster)),
    shell: 
        """    
        module load {params.deeptoolsver}
        if [ ! -d "{params.parent_dir}" ]; then mkdir "{params.parent_dir}"; fi
        if [ \"""" + str(paired_end) + """\" == False ]; then
           extension_option="-e 200"
        else
           extension_option=""
        fi
        
        plotFingerprint -b {input} --labels {params.labels} -p {threads} --skipZeros \\
                        --outQualityMetrics {output.metrics} --plotFile {output.image} --outRawCounts {output.raw} \\
                        ${{extension_option}}
        """


rule deeptools_gene_all:
    input:
        [ join(bw_dir, name + ".Q5DD.RPGC.bw") for name in samples ] 
    output:
        TSSline=join(deeptools_dir,"TSS_profile.Q5DD.pdf"),
        TSSmat=temp(join(deeptools_dir,"TSS.Q5DD.mat.gz")),
        bed=temp(join(deeptools_dir,"geneinfo.Q5DD.bed")),
	mqc=join(deeptools_dir,"TSS_profile.Q5DD.tab")
    params:
        rname                   = "deeptools_gene_all",
        parent_dir              = deeptools_dir,
        deeptoolsver            = config['tools']['DEEPTOOLSVER'],
        # this should be the sample names to match the bigwigs in the same order
        labels                  = samples,
        prebed                  = config['references'][genome]['GENEINFO'],
    threads: 4
    # eventually threads should be 16
    shell: 
        """    
        module load {params.deeptoolsver}
        if [ ! -d "{params.parent_dir}" ]; then mkdir "{params.parent_dir}"; fi
        grep --line-buffered 'protein_coding' {params.prebed} | awk -v OFS='\t' -F'\t' '{{print $1, $2, $3, $5, ".", $4}}' > {output.bed}
       computeMatrix reference-point -S {input} -R {output.bed} -p {threads} \\
                     --referencePoint TSS --upstream 3000 --downstream 3000 --skipZeros \\
                     -o {output.TSSmat} --samplesLabel {params.labels}
       plotProfile -m {output.TSSmat} -out {output.TSSline} \\
                   --yAxisLabel 'average RPGC' --plotType 'se' --legendLocation upper-left \\
                   --numPlotsPerRow 5 --outFileNameData {output.mqc}
        """


rule FRiP_macsN:
    input:
        bed                     = expand(join(macsN_dir, "{name}", "{name}_peaks.narrowPeak"), name=chips),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        join(workpath, "PeakQC", "macsNarrow.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_macsN",
        outroot                 = join(peakqc_dir, "macsNarrow"),
        script                  = join(bin_path, "frip.py"),
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
        export TMPDIR="{params.tmpdir}"
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.script} \\
            -p {input.bed} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {params.outroot} \\
            -x {threads}
        """


rule FRiP_Genrich:
    input:
        bed                     = expand(join(genrich_dir, "{name}", "{name}.narrowPeak"), name=chips),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        join(workpath, "PeakQC", "Genrich.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_macsN",
        outroot                 = join(peakqc_dir, "macsNarrow"),
        script                  = join(bin_path, "frip.py"),
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

        python {params.script} \\
            -p {input.bed} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {params.outroot} \\
            -x {threads}
        """


rule FRiP_macsB:
    input:
        bed                     = expand(join(macsB_dir, "{name}", "{name}_peaks.broadPeak"), name=chips),
        bam                     = join(bam_dir, "{name}.Q5DD.bam"),
    output:
        join(workpath, "PeakQC", "macsBroad.{name}.Q5DD.FRiP_table.txt"),
    params:
        rname                   = "FRiP_macsB",
        outroot                 = join(peakqc_dir, "macsBroad"),
        script                  = join(bin_path, "frip.py"),
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
        export TMPDIR="{params.tmpdir}"
        trap 'rm -rf "${{tmp}}"' EXIT

        python {params.script} \\
            -p {input.bed} \\
            -b {input.bam} \\
            -g {params.genome} \\
            -o {params.outroot} \\
            -x {threads}
        """


rule jaccard:
    input:
        lambda w: [ join(workpath, w.PeakTool, chip, chip + PeakExtensions[w.PeakTool]) for chip in chips ],
    output:
        join(qc_dir, '{PeakTool}_jaccard.txt'),
    params:
        rname                   = "jaccard",
        outroot                 = lambda w: join(qc_dir, w.PeakTool),
        script                  = join(bin_path, "jaccard_score.py"),
        genome                  = config['references'][genome]['REFLEN']
    envmodules:
        config['tools']['BEDTOOLSVER']
    shell: 
        """
        python {params.script} \\
            -i "{input}" \\
            -o "{params.outroot}" \\
            -g {params.genome}
        """