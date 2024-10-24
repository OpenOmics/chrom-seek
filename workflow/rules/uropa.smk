# ~~ Differential binding analysis rules ~~
from os.path import join, sep
from textwrap import dedent

localrules: UROPA_prep_in_macsB, UROPA_prep_in_macsN, \
    UROPA_prep_diffbind_edgeR, UROPA_prep_diffbind_DeSeq2

rule UROPA_prep_diffbind_edgeR:
    input: 
        join(
            diffbind_dir, 
            "{group1}_vs_{group2}-{PeakTool}", 
            "{group1}_vs_{group2}-{PeakTool}_Diffbind_fullList.bed"
        ),
    params:
        rname                           = "UROPA_prep_in_diffbind",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [join(
                                                uropa_diffbind_dir, 
                                                "{group1}_vs_{group2}.{PeakTool}.DiffBind.EdgeR." + pktype + ".json"
                                          ) for pktype in peak_types],               
    log: join(local_log_dir, "UROPA_prep_in_diffbind", "{group1}_vs_{group2}-{PeakTool}_diffbind_prep.log")
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -b {input} \\
            -t {params.peak_types}
        """)


rule UROPA_prep_diffbind_DeSeq2:
    input: 
        join(
            diffbind_dir, 
            "{group1}_vs_{group2}-{PeakTool}", 
            "{group1}_vs_{group2}-{PeakTool}_Diffbind_fullList.bed"
        ),
    params:
        rname                           = "UROPA_prep_in_diffbind",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [join(
                                                uropa_diffbind_dir, 
                                                "{group1}_vs_{group2}.{PeakTool}.DiffBind.DeSeq2." + pktype + ".json"
                                          ) for pktype in peak_types],
                                          
    log: join(local_log_dir, "UROPA_prep_in_diffbind", "{group1}_vs_{group2}-{PeakTool}_diffbind_prep.log")
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -b {input} \\
            -t {params.peak_types}
        """)


rule UROPA_diffbind:
    input: 
        json                            = join(
                                            uropa_dir, 
                                            "DiffBind", 
                                            "{group1}_vs_{group2}.{PeakTool}.DiffBind.{differential_app}.{_type}.json"
                                          ),
    params:
        rname                           = "UROPA_diffbind",
        outroot                         = join(uropa_diffbind_dir, "{group1}_vs_{group2}-{PeakTool}-{differential_app}"),
        uropa_prefix                    = "{group1}_vs_{group2}_{PeakTool}_{differential_app}_{_type}_uropa",
    output: 
        join(
            uropa_diffbind_dir, 
            "{group1}_vs_{group2}-{PeakTool}-{differential_app}",
            "{group1}_vs_{group2}_{PeakTool}_{differential_app}_{_type}_uropa_allhits.txt"
        ),
    threads: int(allocated("threads", "UROPA_diffbind", cluster)),
    container: config["images"]["uropa"]
    shell: 
        dedent("""
            uropa \\
                -i {input} \\
                -o {params.outroot} \\
                -p {params.uropa_prefix} \\
                -t {threads} -s
            """)

# ~~ macs Narrow peak annotation  ~~ #
rule UROPA_prep_in_macsN:
    input:
        join(macsN_dir, "{name}", "{name}_peaks.narrowPeak")
    params:
        rname                           = "UROPA_prep_in_macsN",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [join(
                                            uropa_dir, 
                                            "macsNarrow", 
                                            "{name}.macsNarrow." + pktype + ".json"
                                          ) for pktype in peak_types],
    log: join(local_log_dir, "UROPA_prep_in_macsN", "{name}_uropa_prep.log")
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -b {input} \\
            -t {params.peak_types}
        """)


rule UROPA_macsN:
    input: join(uropa_dir, "macsNarrow", "{name}.macsNarrow.{_type}.json")
    params: 
        rname                           = "UROPA_macsN",
        outroot                         = join(uropa_dir, "macsNarrow"),
    output: join(uropa_dir, "macsNarrow", "{name}_uropa_{_type}_allhits.txt"),
    threads: int(allocated("threads", "UROPA_macsN", cluster)),
    container: config["images"]["uropa"]
    shell: "uropa -i {input} -p {wildcards.name}_uropa_{wildcards._type} -o {params.outroot} -t {threads} -s"


# ~~ macs Broad peak annotation ~~ #
rule UROPA_prep_in_macsB:
    input: 
        join(macsB_dir, "{name}", "{name}_peaks.broadPeak"),
    params:
        rname                           = "UROPA_prep_in_macsB",
        this_script                     = join(bin_path, "uropa_input.py"),
        this_gtf                        = gtf,
        this_assay                      = assay,
        peak_types                      = ' '.join(peak_types),
    output:
        this_json                       = [join(
                                            uropa_dir, 
                                            "macsBroad", 
                                            "{name}.macsBroad." + pktype + ".json"
                                          ) for pktype in peak_types],
    shell:
        dedent("""
        {params.this_script} \\
            -g {params.this_gtf} \\
            -o {output.this_json} \\
            -a {params.this_assay} \\
            -b {input} \\
            -t {params.peak_types}
        """)


rule UROPA_macsBroad:
    input: join(uropa_dir, "macsBroad", "{name}.macsBroad.{_type}.json")
    params:
        rname = "UROPA_macsB",
        outroot = join(uropa_dir, "macsBroad"),
    output: join(uropa_dir, "macsBroad", "{name}_uropa_{_type}_allhits.txt"),
    threads: int(allocated("threads", "UROPA_macsBroad", cluster)),
    container: config["images"]["uropa"],
    shell: "uropa -i {input} -p {wildcards.name}_uropa_{wildcards._type} -o {params.outroot} -t {threads} -s"