# Dockerfile for cfChIP Toolkit
# https://www.nature.com/articles/s41587-020-00775-6
# Uses Ubuntu Focal (20.04 LTS) as base image
# https://github.com/OpenOmics/chrom-seek/commit/48c3903261fe1c5979b2f6b62beb03cbd62c7468
FROM skchronicles/cfchip_toolkit:v0.3.0
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="skchronicles/cfchip_toolkit:v0.3.0" \
    version="v0.4.0"   \
    software="chrom-seek_cfchip_toolkit" \
    about.summary="An awesome set of Epigenetic Pipeline" \
    about.home="https://github.com/OpenOmics/chrom-seek" \
    about.documentation="https://openomics.github.io/chrom-seek/" \
    about.tags="cell-free ChIP-seq"

############# OVERVIEW ##################
# Adds the following packages:
#  • R/4.X
#    - KEGG.db              (BioConductor/3.11)


############### INSTALL #################
# Extra packages for clusterProfiler,
# installing the KEGG database locally
# to avoid making requests to their 
# servers. This avoids the issues and
# any potential failures when we try
# to run KEGG and their servers are
# down or overwhelmed. clusterProfiler
# can use an local KEGG db install
# with: 'use_internal_data=T'. The
# old R package KEGG.db cannot be 
# installed anymore using newer
# versions of R/bioconductor; however,
# the older archive version tied to 
# Bioconductor/3.11 can be installed.
WORKDIR /opt2
RUN wget "https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz" \
    && R CMD INSTALL KEGG.db_3.2.4.tar.gz \
    && rm KEGG.db_3.2.4.tar.gz


################ POST ###################
# Add Dockerfile and export env variables
RUN mv /opt2/Dockerfile /opt2/cfchip_toolkit-v0.4.0.dockerfile
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2
