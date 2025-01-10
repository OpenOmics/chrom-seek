# Quality control metrics and reports
# ~~~~
# Paired-end peak qc rules
from scripts.common import allocated

# ~~ workflow configuration
ends                            = [1, 2]

# ~~ directories
trim_dir                        = join(workpath, 'trim')
kraken_dir                      = join(workpath, 'kraken')


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
        first_txt               = expand(join(
                                    qc_dir, "FQscreen", "{name}.R{rn}.trim_screen.txt"), name=samples, rn=[1, 2]
                                  ),
        first_png               = expand(join(
                                    qc_dir, "FQscreen", "{name}.R{rn}.trim_screen.png"), name=samples, rn=[1, 2]
                                  ),
        second_txt              = expand(join(
                                    qc_dir, "FQscreen2", "{name}.R{rn}.trim_screen.txt"), name=samples, rn=[1, 2]
                                  ),
        second_png              = expand(join(
                                    qc_dir, "FQscreen2", "{name}.R{rn}.trim_screen.png"), name=samples, rn=[1, 2]
                                  )
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
        fq2                     = join(trim_dir, "{name}.R2.trim.fastq.gz")
    output:
        krakenout               = join(kraken_dir, "{name}.trim.kraken_bacteria.out.txt"),
        krakentaxa              = join(kraken_dir, "{name}.trim.kraken_bacteria.taxa.txt"),
        kronahtml               = join(kraken_dir, "{name}.trim.kraken_bacteria.krona.html"),
    params:
        rname                   = 'kraken',
        outdir                  = kraken_dir,
        bacdb                   = config['shared_resources']['KRAKENBACDB'],
        tmpdir                  = tmpdir
    threads: 
        int(allocated("threads", "kraken", cluster)),
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

        kraken2 --db ${{tmp}}/${{kdb_base}} \\
            --threads {threads} --report {output.krakentaxa} \\
            --output {output.krakenout} \\
            --gzip-compressed \\
            --paired {input.fq1} {input.fq2}
        
        # Generate Krona Report
        cut -f2,3 {output.krakenout} | \\
            ktImportTaxonomy - -o {output.kronahtml}
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
        bam                     = join(bam_dir, "{name}.{ext}.bam")
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