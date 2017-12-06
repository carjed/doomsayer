# FROM jupyter/all-spark-notebook:c7fb6660d096
# FROM rocker/tidyverse:3.4.2
# FROM jupyter/r-notebook

# Modified from https://github.com/jupyter/docker-stacks/blob/master/r-notebook/Dockerfile
# Distributed under the terms of the Modified BSD License.
FROM jupyter/minimal-notebook:c7fb6660d096

LABEL maintainer="Jedidiah Carlson <jed.e.carlson@gmail.com>"

# add contents of repo to ${HOME}
COPY . ${HOME}

# go to root to own dir and pre-requisite installs
USER root
RUN chown -R ${NB_UID}:${NB_UID} ${HOME}

# R pre-requisites
RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  fonts-dejavu \
  tzdata \
  gfortran \
  file \
  git \
  libapparmor1 \
  libcurl4-openssl-dev \
  libedit2 \
  libssl-dev \
  lsb-release \
  psmisc \
  python-setuptools \
  sudo \
  wget \
  libxml2-dev \
  libcairo2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadb-client-lgpl-dev \
  libpq-dev \
  libssh2-1-dev \
  gcc && apt-get clean && \
  rm -rf /var/lib/apt/lists/*

# switch to user for R install
USER ${NB_USER}

# R packages
RUN conda install --quiet --yes \
    'r-base=3.4.2' \
    'r-irkernel=0.7*' \
    'r-devtools=1.13*' \
    'r-tidyverse=1.2*' && \
    conda clean -tipsy && \
    fix-permissions $CONDA_DIR

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
