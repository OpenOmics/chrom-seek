<div align="center">
   
  <h1>chrom-seek 🔬</h1>
  
  **_An awesome set of epigenetic pipelines_**

   [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.16921008-blue)](https://doi.org/10.5281/zenodo.16921008) [![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/OpenOmics/chrom-seek?color=blue&include_prereleases)](https://github.com/OpenOmics/chrom-seek/releases) [![Docker Pulls](https://img.shields.io/docker/pulls/skchronicles/cfchip_toolkit)](https://hub.docker.com/r/skchronicles/cfchip_toolkit) <br> [![tests](https://github.com/OpenOmics/chrom-seek/workflows/tests/badge.svg)](https://github.com/OpenOmics/chrom-seek/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/chrom-seek/workflows/docs/badge.svg)](https://github.com/OpenOmics/chrom-seek/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/chrom-seek?color=brightgreen)](https://github.com/OpenOmics/chrom-seek/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/chrom-seek)](https://github.com/OpenOmics/chrom-seek/blob/main/LICENSE) 
  
  <i>
    This is the home of the pipeline, chrom-seek. Its long-term goals: to accurately call and annotate peaks, to infer cell types in cell-free samples, and to boldly quantify diferential binding or accessibility like no pipeline before!
  </i>
</div>

## Overview

Welcome to chrom-seek! Before getting started, we highly recommend reading through [chrom-seek's documentation](https://openomics.github.io/chrom-seek/).

The **`./chrom-seek`** pipeline is composed of several interrelated sub-commands to set up and run the pipeline across different systems. Each of the available sub-commands performs different functions: 

 * [<code>chrom-seek <b>run</b></code>](https://openomics.github.io/chrom-seek/usage/run/): Run the chrom-seek pipeline with your input files.
 * [<code>chrom-seek <b>unlock</b></code>](https://openomics.github.io/chrom-seek/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>chrom-seek <b>install</b></code>](https://openomics.github.io/chrom-seek/usage/install/): Download reference files locally.
 * [<code>chrom-seek <b>cache</b></code>](https://openomics.github.io/chrom-seek/usage/cache/): Cache remote resources locally, coming soon!

**chrom-seek** is an awesome set of pipelines designed specifically for cell-free ChIP-seq, bulk ChIP-seq, and bulk ATAC-seq sequencing data. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/chrom-seek/usage/run/) section of each available sub-command.

For more information about issues or troubleshooting a problem, please check out our [FAQ](https://openomics.github.io/chrom-seek/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/chrom-seek/issues).

## Dependencies

**Requires:** `singularity>=3.5`  `snakemake<=7.32.4`

At the current moment, the pipeline uses a mixture of environment modules and docker images; however, this will be changing soon! In the very near future, the pipeline will only use docker images. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline will rely on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation

Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/chrom-seek.git
# Change your working directory
cd chrom-seek/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./chrom-seek -h
```

## Contribute 

This site is a living document, created for and by members like you. chrom-seek is maintained by the members of OpenOmics and is improved by continuous feedback! We encourage you to contribute new content and make improvements to existing content via pull requests to our [GitHub repository](https://github.com/OpenOmics/chrom-seek).


## Cite

Please note that more citation styles and releases can be found on the chrom-seek [zenodo page](https://doi.org/10.5281/zenodo.16921008).

If you use this pipeline, please cite it as below:

<details>
  <summary><b><i>@BibText</i></b></summary>
 
```text
@software{routsong_2026_16921008,
  author       = {Markowitz, Tovah and
                  Routsong, Ryan and
                  Khleborodova, Asya and
                  Kuhn, Skyler},
  title        = {OpenOmics/chrom-seek},
  month        = aug,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.16921008},
  url          = {https://doi.org/10.5281/zenodo.16921008},
}
```

</details>

<details>
  <summary><b><i>@APA</i></b></summary>

```text
Markowitz, T., Routsong, R., Khleborodova, A., & Kuhn, S. (2026). OpenOmics/chrom-seek [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.16921008
```

</details>


## References

<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
