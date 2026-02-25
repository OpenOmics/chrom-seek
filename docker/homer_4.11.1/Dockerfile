# https://github.com/mattgalbraith/homer-docker-singularity
################## BASE IMAGE ######################
FROM ubuntu:22.04 AS build

################## 1st STAGE: SAMTOOLS INSTALLATION ######################
ARG ENV_NAME="samtools"
ARG VERSION="1.16.1"

ENV DEBIAN_FRONTEND=noninteractive
ENV PACKAGES="ca-certificates gcc mono-mcs libncurses5-dev libncursesw5-dev zlib1g zlib1g-dev bzip2 libbz2-dev liblzma-dev \
    libhtscodecs2 build-essential wget libcurl4 libcurl4-gnutls-dev"
# need libcurl4 and libcurl4-gnutls-dev for GCS access - see https://github.com/samtools/samtools/issues/862

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Get, configure, make, and install samtools 
RUN wget https://github.com/samtools/samtools/releases/download/${VERSION}/samtools-${VERSION}.tar.bz2 && \
    bzip2 -d samtools-${VERSION}.tar.bz2 && \
    tar xvf samtools-${VERSION}.tar && \
    cd samtools-${VERSION} && \
    ./configure && \
    make && \
    make install

################## 2ND STAGE: HOMER INSTALLATION ######################
FROM ubuntu:22.04

################## METADATA ######################
LABEL base_image="ubuntu:22.04"
LABEL version="1.1.0"
LABEL software="Homer"
LABEL software.version.homer="4.11.1"
LABEL software.version.samtools="1.9.0"
LABEL about.summary="HOMER (Hypergeometric Optimization of Motif EnRichment) is a suite of tools for Motif Discovery and next-gen sequencing analysis."
LABEL about.home="http://homer.ucsd.edu/homer/"
LABEL about.documentation="http://homer.ucsd.edu/homer/"
LABEL about.license_file=""
LABEL about.license.homer="GPLv3"
LABEL about.license.samtools="MIT/Expat"

################## MAINTAINER ######################
LABEL maintainer="Matthew Galbraith <matthew.galbraith@cuanschutz.edu>"

ENV DEBIAN_FRONTEND=noninteractive
ENV PACKAGES="libnss-sss r-base r-base-dev tabix wget ca-certificates libssl-dev libcurl4-openssl-dev libxml2-dev \
    mono-mcs libncurses5-dev libncursesw5-dev zlib1g zlib1g-dev bzip2 libbz2-dev liblzma-dev libhtscodecs2"
# dependency samtools installed and copied from 1st stage to control version and save space

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /usr/local/bin/samtools/ /usr/local/bin

# Install required Bioconductor package(s) # SLOW and MAY NOT NEED THESE for tagdir and ucsc file creation + may want to control R version? (looks like 3.6.3 is default)
RUN R -e 'install.packages(c("BiocManager", "edgeR", "DESeq2"))'

# Install and configure Homer
# RUN mkdir /opt/homer/ && cd /opt/homer && wget http://homer.ucsd.edu/homer/configureHomer.pl && /usr/bin/perl configureHomer.pl -install # defaults to latest version
RUN mkdir /opt/homer/ && cd /opt/homer && wget http://homer.ucsd.edu/homer/configureHomer.pl && /usr/bin/perl configureHomer.pl -install homer -version v4.11.1 # should install specific version

# NOT USED: Genome & Annotation data is too large
# RUN cd /opt/homer && /usr/bin/perl configureHomer.pl -install hg38 && /usr/bin/perl configureHomer.pl -install human-o

# Genome & Annotation data is too large to include in this image, so it
# requires a directory to be mounted at /opt/homerdata/ containing
#  - config.txt - homer configuration file with directories pointing to paths like "data/accession"
#  - data  - folder containing homer annotation data files

# Add softlinks for config file and data directory
RUN rm -rf /opt/homer/data && ln -s /opt/homerdata/data /opt/homer/data
RUN rm -f /opt/homer/config.txt && ln -s /opt/homerdata/config.txt /opt/homer/config.txt

# Add homer installation dir to PATH
ENV PATH=${PATH}:/opt/homer/bin/