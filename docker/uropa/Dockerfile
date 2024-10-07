FROM conda/miniconda3
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict
RUN conda install bioconda::uropa
