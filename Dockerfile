#############################################################
# Dockerfile to build a running env container for starscope #
#############################################################

## From conda docker
FROM continuumio/miniconda3

## Maintainer
MAINTAINER oben <obennoname@gmail.com>

## Setup conda env
USER root
## install tini for init process
RUN apt update && apt install -y git tini
## use mamba instead of conda command
RUN conda install mamba -n base -c conda-forge

## Setup workdir
WORKDIR /app

#### download app source
##RUN git clone https://github.com/obenno/scATAC-seq ./scATAC-seq
##WORKDIR /app/scATAC-seq
## create conda env with app env file
##RUN mamba env create -f scRNAseq_env.yml

## All packages will be installed to base env
## and could be invoked directly
COPY conda/scATAC_packages.txt .
RUN mamba install -n base --file scATAC_packages.txt

## install SeuratDisk
RUN mamba install -c bioconda -c conda-forge r-hdf5r r-cli r-crayon r-matrix r-r6 r-rlang r-withr r-stringi r-remotes
RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk", lib = "/opt/conda/lib/R/library")'

## copy entrypoint.sh
##COPY entrypoint.sh .

## Follow Dockstore's guide
## switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -m -r -g ubuntu -u 1000 ubuntu
RUN chown -R ubuntu: /app
USER ubuntu

## Setup bashrc file for the ubuntu user
##RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
##RUN echo "conda activate starscope_env" >> ~/.bashrc

##SHELL ["/bin/bash", "--login", "-c"]
## setup default cmd (deprecated)
## no entrypoint and cmd will be set
##RUN chmod +x entrypoint.sh
##ENTRYPOINT ["./entrypoint.sh"]
