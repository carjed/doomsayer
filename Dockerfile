# FROM jupyter/all-spark-notebook:c7fb6660d096
# FROM rocker/tidyverse:3.4.2
# FROM jupyter/r-notebook

# Modified from https://github.com/jupyter/docker-stacks/blob/master/r-notebook/Dockerfile
# Distributed under the terms of the Modified BSD License.
FROM jupyter/scipy-notebook:c7fb6660d096
# FROM jupyter/minimal-notebook:033056e6d164

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

LABEL maintainer="Jupyter Project <jupyter@googlegroups.com>"

USER root

# R pre-requisites
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    fonts-dejavu \
    tzdata \
    gfortran \
    gcc && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# R packages including IRKernel which gets installed globally.
RUN conda config --system --append channels r && \
    conda install --quiet --yes \
    'r-base=3.4.2' \
    'r-irkernel=0.8*' \
    'r-devtools=1.12*' \
    'r-tidyverse=1.1*' \
    'r-shiny=1.0*' \
    'r-rmarkdown=1.6*' \
    'r-rcurl=1.95*' && \
    conda clean -tipsy

USER ${NB_USER}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

# run install.r script to load R package dependencies
# RUN R --quiet -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"
RUN if [ -f install.r ]; then R --quiet -f install.r; fi

# FROM jupyter/scipy-notebook:c7fb6660d096
# ADD pip_reqs.txt ./
# ADD env.yml ./

# RUN conda create -n doomsayer-environment python=3.6 anaconda
# RUN source activate doomsayer-environment
# RUN conda env export > environment.yml

# create conda environment and install Python library dependencies
RUN conda env create -n doomsayer -f env.yml

# activate the conda environment
ENV PATH /opt/conda/envs/doomsayer/bin:$PATH

# old version; installs via pip
# causes problems with cyvcf2--better to use bioconda
# RUN pip install -r pip_reqs.txt
