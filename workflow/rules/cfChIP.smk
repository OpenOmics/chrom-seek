# cell-free DNA ChIP-seq rules:
#   - picard_dedup 
#   - cfChIPtool
#   - cfChIPcompile

rule picard_dedup:
    input: 
        bam2=join(workpath,bam_dir,"{name}.Q5.bam")
    output:
        out5=join(workpath,bam_dir,"{name}.Q5DD.bam"),
        out5f=join(workpath,bam_dir,"{name}.Q5DD.bam.flagstat"),
        out5i=join(workpath,bam_dir,"{name}.Q5DD.bam.idxstat"),
        out6=join(workpath,bam_dir,"{name}.bwa.Q5.duplic"),
        out7=temp(join(workpath,bam_dir,"{name}.Q5DD.tagAlign"))
    params:
        rname='dedup',
        picardver=config['tools']['PICARDVER'],
        samtoolsver=config['tools']['SAMTOOLSVER'],
        rver=config['tools']['RVER'],
        javaram='16g',
        tmpBam="{name}.Q5DD.withXY.bam",
        rscript=join(config['references'][genome]['cfChIP_TOOLS_SRC'], "bam2fragment.R")
    shell: """
    module load {params.samtoolsver};
    module load {params.picardver};
    module load {params.rver}; 
    if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
    cd /lscratch/$SLURM_JOBID
    
    java -Xmx{params.javaram} \\
      -jar $PICARDJARPATH/picard.jar MarkDuplicates \\
      INPUT={input.bam2} \\
      OUTPUT={params.tmpBam} \\
      TMP_DIR=/lscratch/$SLURM_JOBID \\
      VALIDATION_STRINGENCY=SILENT \\
      REMOVE_DUPLICATES=true \\
      METRICS_FILE={output.out6}
    
    samtools index {params.tmpBam}
    samtools view -b {params.tmpBam} chr{{1..22}} > {output.out5}
    samtools index {output.out5}
    samtools flagstat {output.out5} > {output.out5f}
    samtools idxstats {output.out5} > {output.out5i}
    
    Rscript {params.rscript} {params.tmpBam} {output.out7}
    """


rule cfChIPtool:
    input: 
        join(workpath,bam_dir,"{name}.Q5DD.tagAlign")
    output:
        out1=join(workpath,cfTool_subdir2,"{name}.Q5DD.tagAlign.gz"),
        out2=join(workpath,cfTool_dir,"Output","H3K4me3","Signatures","{name}.Q5DD.csv"),
        out3=join(workpath,cfTool_dir,"Samples","H3K4me3","{name}.Q5DD.rdata"),
    params:
        rname='cfChiP',
        rver="R/4.1.0",
        toolkit = config['references']['cfChIP_TOOLS_SRC'],
        tmpfile = lambda w: join(workpath,cfTool_subdir2, w.name + ".Q5DD.tagAlign"),
    container:
        config['images']['cfchip']
    shell: """
    cp {input} {params.tmpfile}
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
        expand(join(workpath,cfTool_dir,"Output","H3K4me3","Signatures","{name}.Q5DD.csv"),name=H3K4me3samples)
    output:
        txt=join(workpath,qc_dir,"H3K4me3_cfChIP_signature.txt"),
        pdf=join(workpath,qc_dir,"H3K4me3_cfChIP_signature.pdf")
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
