# Dockerfile for cfChIP Toolkit
# https://www.nature.com/articles/s41587-020-00775-6
# Using Ubuntu Focal (20.04 LTS) as base image
FROM ubuntu:20.04
LABEL maintainer=kuhnsa@nih.gov

#  • R/4.1
#    - Biobase            (Bioconductor)
#    - cba                (CRAN)
#    - circlize           (CRAN)
#    - clusterProfiler    (Bioconductor)
#    - ComplexHeatmap     (Bioconductor)
#    - cowplot            (CRAN)
#    - ctc                (Bioconductor)
#    - DiffBind           (Bioconductor/3.11 Archive, requires <3.X)
#    - enrichplot         (Bioconductor)
#    - enrichR            (CRAN)
#    - GenomicAlignments  (Bioconductor)
#    - GenomicRanges      (Bioconductor)
#    - ggforce            (CRAN)
#    - ggplot2            (CRAN)
#    - ggpubr             (CRAN)
#    - ggrepel            (CRAN)
#    - MASS               (CRAN)
#    - Matrix             (CRAN)
#    - NMF                (CRAN)
#    - optparse           (CRAN)
#    - preprocessCore     (Bioconductor)
#    - RColorBrewer       (CRAN)
#    - ReactomePA         (Bioconductor)
#    - reshape2           (CRAN)
#    - rtracklayer        (Bioconductor)
#    - tools              (r-core)

# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with host filesyetem 
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apt-transport-https \
        autotools-dev \
        automake \
        autoconf \
        build-essential \
        bzip2 \
        ca-certificates \
        cmake \
        curl \
        dirmngr \
        g++ \
        gcc \
        gfortran \
        git \
        gnupg \
        gzip \
        libcairo2-dev \
        libfontconfig1-dev \
        libgeos-dev \
        libhdf5-dev \
        libtool \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
        locales \
        make \
        pandoc \
        python3 \
        python3-pip \
        software-properties-common \
        samtools \
        unzip \
        wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install python packages
RUN ln -sf /usr/bin/python3.8 /usr/bin/python \
    && ln -sf /usr/bin/python3.8 /usr/bin/python3
RUN pip3 install --upgrade pip \
    && pip3 install argparse \
    && pip3 install pandas \
    && pip3 install pysam \
    && pip3 install scipy \
    && pip3 install matplotlib \
    && pip3 install XlsxWriter

# Install R (4.X) and required packages
# For more information, check out: 
# https://cran.r-project.org/bin/linux/ubuntu/
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && apt-get -y install r-base r-base-core r-recommended r-base-dev \
    && apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev
# Install required R packages
# CRAN packages
RUN Rscript -e 'install.packages(c("Bioconductor", "argparse", "cba", "circlize", "cowplot", "enrichR", "ggforce", "ggplot2", "ggpubr", "ggrepel", "MASS", "Matrix", "NMF", "optparse", "RColorBrewer", "reshape2"), repos="http://cran.r-project.org")'
# Bioconductor packages
RUN Rscript -e 'BiocManager::install(c("Biobase", "clusterProfiler", "ComplexHeatmap", "ctc", "enrichplot", "GenomicAlignments", "GenomicRanges", "preprocessCore", "ReactomePA", "rtracklayer"))'

# Install DiffBind from Bioconductor (3.11) Archive
# This version of Bioconductor is tied to R/4.0
# and it contains DiffBind version, 2.15.2 
RUN Rscript -e 'install.packages(c("locfit", "gplots", "NMF"), repos="http://cran.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("edgeR", "limma", "systemPipeR", "DESeq2"))'
RUN wget "https://bioconductor.org/packages/3.11/bioc/src/contrib/Archive/DiffBind/DiffBind_2.15.2.tar.gz" \
    && R CMD INSTALL DiffBind_2.15.2.tar.gz \
    && rm DiffBind_2.15.2.tar.gz

# Extra packages for new DiffBind Rmarkdown report
RUN Rscript -e 'install.packages(c("VennDiagram", "ggprism", "DT"), repos="http://cran.r-project.org")'
# Bioconductor packages
RUN Rscript -e 'BiocManager::install(c("org.Hs.eg.db"))'

# Add Dockerfile and export env variables
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2
