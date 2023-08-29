# Python image for chrom-seek,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

LABEL maintainer="Asya Khleborodova <asya.khleborodova@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="chrom-seek" \
    about.summary="An awesome set of epigenetic pipelines" \
    about.home="https://github.com/OpenOmics/chrom-seek" \
    about.documentation="openomics.github.io/chrom-seek/" \
    about.tags="python for chip-seek"

# Overview of Dependencies

#  • python version >= 3:            python3 (3.10.4)
#    • Python Package: optparse, pybedtools, pysam, pandas,
#                      argparse, scipy, scikit-learn, matplotlib,
#                      seaborn                      


############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        autotools-dev \
        automake \
        autoconf \
        cmake \
        curl \
        gawk \
        git \
        grep \
        gzip \
        locales \
        libtool \
        make \
        unzip \
        wget \
        vcftools \
        bedtools \
        samtools \
        python3 \
        python3-pip \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

############### INSTALL ################
# Make python3 the default interpreter
# and install Python Packages 
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
	&& pip3 install argparse \
    && pip3 install pybedtools \
    && pip3 install pysam \
    && pip3 install pandas \
    && pip3 install scipy \
    && pip3 install scikit-learn \
    && pip3 install matplotlib \
    && pip3 install seaborn

################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
COPY Dockerfile /opt2
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2