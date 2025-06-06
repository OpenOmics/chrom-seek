# Quality control metrics and reports
# ~~~~
# Single-ended quality control rules
from scripts.common import allocated

# ~~ workflow configuration
ends                            = [1]

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
        first_txt               = expand(join(qc_dir, "FQscreen", "{name}.R1.trim_screen.txt"), name=samples),
        first_png               = expand(join(qc_dir, "FQscreen", "{name}.R1.trim_screen.png"), name=samples),
        second_txt              = expand(join(qc_dir, "FQscreen2", "{name}.R1.trim_screen.txt"), name=samples),
        second_png              = expand(join(qc_dir, "FQscreen2", "{name}.R1.trim_screen.png"), name=samples),
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
    output:
        krakentaxa              = join(kraken_dir, "{name}.trim.kraken_bacteria.taxa.txt"),
        krakenout               = temp(join(kraken_dir, "{name}.trim.kraken_bacteria.out.txt")),
    params:
        rname                   = 'kraken',
        outdir                  = kraken_dir,
        bacdb                   = config['shared_resources']['KRAKENBACDB'],
        tmpdir                  = tmpdir,
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

        kraken2 --db ${{tmp}}/${{kdb_base}} \\
            --threads {threads} --report {output.krakentaxa} \\
            --output {output.krakenout} \\
            --gzip-compressed \\
            {input.fq1}
        """
