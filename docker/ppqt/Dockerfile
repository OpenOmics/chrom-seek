# Base image for PPQT/v1.2.2-6984a713,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of PhantomPeakQualTools:
#  - https://github.com/kundajelab/phantompeakqualtools#dependencies
# Overview:
#  -+ samtools                  # apt-get, installs: 1.13-4
#  -+ gawk                      # apt-get, installs: 1:5.1.0
#  -+ libboost-all-dev          # apt-get, installs: 1.74.0
#  -+ R>=3.4                    # apt-get, installs: 4.1.2-1
#     @requires: spp            # from Github ppqt, installs: 1.14
#     @requires: Rsamtools      # from Bioconductor
#     @requires: biocmanager    # from CRAN
#     @requires: devtools       # from CRAN
#     @requires: ggplot2        # from CRAN
#     @requires: tidyr          # from CRAN
#     @requires: caTools        # from CRAN
#     @requires: snow           # from CRAN
#     @requires: snowfall       # from CRAN
#     @requires: bitops         # from CRAN
#  -+ python>3.7                # apt-get, installs: 3.10.12
#     @requires: numpy          # pip install from pypi
#     @requires: pysam          # pip install from pypi
#     @requires: pandas         # pip install from pypi
#     @requires: scipy          # pip install from pypi
LABEL maintainer="kundajelab" \
   base_image="ubuntu:22.04" \
   version="v0.2.0"   \
   software="ppqt/v1.2.2-6984a713" \
   about.summary="Computes informative enrichment and quality measures for ChIP-seq, DNase-seq, MNase-seq data." \
   about.home="https://github.com/kundajelab/phantompeakqualtools" \
   about.documentation="https://github.com/kundajelab/phantompeakqualtools/blob/master/README.md" \
   about.tags="Epigenetics"

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
       cmake \
       curl \
       gawk \
       git \
       gzip \
       g++ \
       libboost-dev \
       libcurl4-openssl-dev \
       libssl-dev \
       libxml2-dev \
       locales \
       # python/3.10.6
       python3 \
       python3-pip \
       # R/4.1.2-1
       r-base \
       r-base-dev \
       r-recommended \
       # R packages (CRAN)
       r-cran-biocmanager \
       r-cran-devtools \
       r-cran-catools \
       r-cran-snow \
       r-cran-snowfall \
       r-cran-bitops \
       r-cran-dplyr \
       r-cran-ggplot2 \
       r-cran-tidyr \
       # R packages (Bioconductor)
       r-bioc-rsamtools \
       # samtools/1.13-4
       samtools \
       wget \
       zlib1g-dev \
   && apt-get clean && apt-get purge \
   && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8


############### INSTALL ################
# Install any bioinformatics packages
# available with pypi or CRAN/BioC
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
   && pip3 install matplotlib \
   && pip3 install numpy \
   && pip3 install pandas \
   && pip3 install scipy


############# PPQT/v1.2.2 ##############
# Install PPQT from src manually,
# Installs ppqt/v1.2.2 via commit:
# 6984a713aba0218b76bacc63f2fb5087425fd6a3
RUN mkdir -p /opt2/ppqt/v1.2.2-6984a713 \
   && cd /opt2/ppqt/v1.2.2-6984a713/ \
   && git clone https://github.com/kundajelab/phantompeakqualtools.git \
   && cd /opt2/ppqt/v1.2.2-6984a713/phantompeakqualtools \
   && git checkout 6984a713aba0218b76bacc63f2fb5087425fd6a3 \
   # Fixes a bug in the current version of
   # ppqt due to undeclared variable name,
   # probably a typo from the authors or 
   # the name of the variable changed at
   # some point in time. It is pretty 
   # harmless and only occurs if it does 
   # not recongize the extension of the
   # input file name (*.bam* / *.tagAlign*) 
   && sed -i 's/align.fname/align.filename/' /opt2/ppqt/v1.2.2-6984a713/phantompeakqualtools/run_spp.R \
   && wget https://github.com/hms-dbmi/spp/archive/refs/tags/1.15.2.tar.gz -P /opt2/ppqt/v1.2.2-6984a713/phantompeakqualtools \
   && Rscript -e 'install.packages("/opt2/ppqt/v1.2.2-6984a713/phantompeakqualtools/1.15.2.tar.gz")' 
ENV PATH="${PATH}:/opt2/ppqt/v1.2.2-6984a713/phantompeakqualtools"
WORKDIR /opt2


################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/ppqt_1-2-2.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2