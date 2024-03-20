# Dockerfile for cfChIP Toolkit
# https://www.nature.com/articles/s41587-020-00775-6
# Uses Ubuntu Focal (20.04 LTS) as base image
# https://github.com/OpenOmics/chrom-seek/commit/48c3903261fe1c5979b2f6b62beb03cbd62c7468
FROM skchronicles/cfchip_toolkit:v0.4.0
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="skchronicles/cfchip_toolkit:v0.3.0" \
    version="v0.5.0"   \
    software="chrom-seek_cfchip_toolkit" \
    about.summary="An awesome set of Epigenetic Pipeline" \
    about.home="https://github.com/OpenOmics/chrom-seek" \
    about.documentation="https://openomics.github.io/chrom-seek/" \
    about.tags="cell-free ChIP-seq"

############# OVERVIEW ##################
# Adds the following packages:
#  • R/4.X
#    - umap              (CRAN)


############### INSTALL #################
# Extra packages for differential bind
# analysis in the cfChIP pipeline 
WORKDIR /opt2
RUN Rscript -e 'install.packages(c("umap"), repos="http://cran.r-project.org")'


################ POST ###################
# Add Dockerfile and export env variables
RUN mv /opt2/Dockerfile /opt2/cfchip_toolkit-v0.5.0.dockerfile
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2
