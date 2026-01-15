# <code>chrom-seek <b>run</b></code>

## 1. About 
The `chrom-seek` executable is composed of several inter-related sub commands. Please see `chrom-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>chrom-seek <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running chrom-seek with one of its available data-processing pipelines. 

Setting up the chrom-seek pipeline is fast and easy! In its most basic form, <code>chrom-seek <b>run</b></code> only has *five required inputs*. To run an available pipeline with your data raw data, please provide a space seperated list of FastQ (globbing is supported), an output directory to store results, a reference genome for alignment and annotation, an assay type to invoke a specific data-processing pipeline, and a peak call file to set sample metadata. 

## 2. Synopsis
```text
$ chrom-seek run [--help] \
      [--mode {slurm,local}] [--job-name JOB_NAME] [--batch-id BATCH_ID] \
      [--tmp-dir TMP_DIR] [--silent] [--sif-cache SIF_CACHE] \ 
      [--singularity-cache SINGULARITY_CACHE] \
      [--dry-run] [--threads THREADS] \
      [--contrasts CONTRASTS] \
      --assay {cfChIP,ChIP,ATAC,cutnrun} \
      --genome GENOME \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --peakcall PEAKCALL
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument and an output directory to store results via `--output` argument, define an assay type to select an appropriate data-processing pipeline via `--assay` argument, select a reference genome to be used for alignment and annotation via `--genome` argument, and a peakcall file to define groups/inputs/blocking factors for each sample.

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--assay {cfChIP,ChIP,ATAC,cutnrun}`  
> **Assay type or data-processing pipeline.**  
> *type: string*  
> 
> This option defines which pipeline will be run. chrom-seek supports the processing of bulk ChIP-seq (ChIP), cell-free DNA ChIP-seq (cfChIP), ATAC-seq (ATAC) samples, and CUT & RUN or CUT & TAG (cutnrun). Please select from one of the following data-processing pipelines: `cfChIP`, `ChIP`, `ATAC`, `cutnrun`.
> 
> ***Example:*** `--assay ChIP`

---
  `--genome {hg19,hg38,mm10,mm39}`  
> **Reference genome.**  
> *type: string*  
> 
> This option defines the reference genome of the samples for alignment and annotation. There are prebuilt reference files for human (hg19 and hg38), mouse (mm10 and mm39), and rhesus data. Please select one of the following options: `hg19`, `hg38`, `mm10`, `mm39`,`rheMac10`.
> 
> ***Example:*** `--genome hg19`

---
  `--input INPUT [INPUT ...]`  
> **Input FastQ.**  
> *type: file(s)*  
> 
> One or more FastQ files can be provided. From the command-line, each input file should seperated by a space. Globbing is supported! This makes selecting FastQ files easy. FastQ files should always be gzipp-ed. Only list files you want processed as all files in the list will be run through the initial pipeline steps. All file merging must be done before running the pipeline.
> 
> ***Example:*** `--input .tests/*.R?.fastq.gz`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/chrom-seek_out`

---  
  `--peakcall PEAKCALL`
> **Peakcall file.**   
> *type: file*
> 
> Path to a sample sheet in TSV format used to map each ChIP/ATAC sample to its group label(s) for downstream comparisons and, for *ChIP-seq only*, optionally pair a sample to an *input control*. The header must include `Sample` and `Group`, and may optionally include `InputControl` and/or `Block` (column names are case-insensitive).
>
> #### Required columns
> * **Sample**: Basename of the sample (derived from the sample’s R1 FastQ by removing the read/extension suffix), e.g. `WT_S4.R1.fastq.gz` becomes `WT_S4` and `WT_S4_R1_001.fastq.gz` becomes `WT_S4`.
> * **Group**: Group label(s) for the sample. Multiple groups may be provided as a comma-separated list. Each sample must be assigned to at least one group. Group names currently cannot include `.`, `-`, or `_`, and group names should **not** be substrings of other group names (e.g. avoid `WT` and `WT_Treated` together).
>
> #### Optional columns
> * **InputControl**: Basename of the corresponding input control sample for the given `Sample` (derived the same way as above). This column is used to pair each ChIP sample to its matched input control for correction during peak calling. *ATAC-seq samples should never provide `InputControl`.*
>* **Block**: Blocking factor used to avoid duplicate correlations between repeated observations (commonly biological replicate ID or subject/individual ID, e.g., multiple samples from the same individual).
>
> 
> 
> **Contents of example peakcalls file:** 
> ```
> Sample	InputControl	Group
> WT_S1	  IN_S1	        G1,G3
> WT_S2	  IN_S2       	G1,G3
> WT_S3	  IN_S3	        G1
> WT_S4	  IN_S4	        G2,G4
> WT_S5	  IN_S5	        G2,G4
> WT_S6	  IN_S6	        G2
> ```
> ***Example:*** `--peakcall /data/$USER/peakcall.tsv`

### 2.2 Analysis options

Each of the following arguments are optional, and do not need to be provided. 

#### 2.2.1 Differential Binding/Accessibility

  `--contrasts CONTRASTS`
> **Contrasts file.**   
> *type: file*
>   
> This tab delimited (TSV) file is used to setup comparisons within different groups of samples. Please see the `--peakcall` option above for more information about how to define groups within a set of samples. This file consists of two columns containing the names of two groups to compare. The names defined in this file must also exist in the peakcall file.  
> 
> *Please note:* the ordering of groups is preserved when creating contrasts. This is important because it dicates how to interpret the direction of the fold-change for your comparison. In the example below, the first comparison can be interpreted as *G2 vs. G1*. This would result in the following contrast: `G2-G1`. Within the context of differential binding analysis, a positive fold-change would indicate that the G2 group has higher levels of binding (cfChIP/ChIP) or accessibility (ATAC) at X region.  
> 
> **Contents of example contrasts file:**  
> ```
> G2 	G1
> G4 	G1
> G4 	G3
> ```
> ***Example:*** `--contrasts /data/$USER/contrasts.tsv` 

### 2.3 Orchestration options

Each of the following arguments are optional, and do not need to be provided. 

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
> 
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,local}`  
> **Execution Method.**  
> *type: string*  
> *default: slurm*
> 
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm or local. 
> 
> ***slurm***    
> The slurm execution method will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running chrom-seek in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution.
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a local execution mode. 
> 
> ***Example:*** `--mode slurm`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*
> *default: pl:chrom-seek*
> 
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "pl:chrom-seek".
> 
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `--output OUTPUT/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
> 
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `chrom-seek cache` subcommand can be used to create a local SIF cache. Please see `chrom-seek cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running chrom-seek with this option when ever possible.
> 
> ***Example:*** `--sif-cache /data/$USER/SIFs`

---  
  `--threads THREADS`   
> **Max number of threads for each process.**  
> *type: int*  
> *default: 2*
> 
> Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
> 
> ***Example:*** `--threads 12`


---  
  `--tmp-dir TMP_DIR`   
> **Max number of threads for each process.**  
> *type: path*  
> *default: `/lscratch/$SLURM_JOB_ID`*
> 
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOB_ID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
> 
> ***Example:*** `--tmp-dir /scratch/$USER/`

### 2.4 Miscellaneous options  
Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example usage on Biowulf
```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
sinteractive -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./chrom-seek run --assay ChIP \
                  --genome hg19 \
                  --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --peakcall .tests/peakcall.tsv \
                  --sif-cache /data/OpenOmics/SIFs/ \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the chrom-seek pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode.
./chrom-seek run --assay ChIP \
                  --genome hg19 \
                  --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --peakcall .tests/peakcall.tsv \
                  --sif-cache /data/OpenOmics/SIFs/ \
                  --mode slurm
```
