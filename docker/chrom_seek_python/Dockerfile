FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Skyler Kuhn <skyler.kuhn@nih.gov>, Ryan Routsong <routsongrm@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.1"   \
    software="chrom-seek" \
    about.summary="An awesome set of epigenetic pipelines" \
    about.home="https://github.com/OpenOmics/chrom-seek" \
    about.documentation="openomics.github.io/chrom-seek/" \
    about.tags="python for chip-seek"

############### INIT ################
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get install -y \
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
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

############### INSTALL ################
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip 
RUN pip3 install pybedtools \
    pysam \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    plotly \
    seaborn

################ POST #################
COPY Dockerfile /opt2
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"